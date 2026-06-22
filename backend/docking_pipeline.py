#!/usr/bin/env python3
"""
Molecular Docking Pipeline
Orchestrates the complete docking workflow
"""

import base64
from pathlib import Path
from typing import Dict, Optional
import asyncio

from .protein_handler import ProteinHandler
from .ligand_handler import LigandHandler
from .vina_runner import VinaRunner


class DockingPipeline:
    """Main docking pipeline orchestrator"""
    
    def __init__(
        self,
        work_dir: Path,
        grid_padding: float = 4.0,
        exhaustiveness: int = 8,
        num_modes: int = 9
    ):
        """
        Initialize docking pipeline
        
        Args:
            work_dir: Working directory for temporary files
            grid_padding: Padding for docking box (Angstroms)
            exhaustiveness: Vina search exhaustiveness
            num_modes: Number of docking modes to generate
        """
        self.work_dir = Path(work_dir)
        self.work_dir.mkdir(parents=True, exist_ok=True)
        
        self.grid_padding = grid_padding
        self.exhaustiveness = exhaustiveness
        self.num_modes = num_modes
        
        # Initialize handlers
        self.protein_handler = ProteinHandler()
        self.ligand_handler = LigandHandler()
        self.vina_runner = VinaRunner()
    
    async def run_full_pipeline(
        self,
        pdb_id: str,
        ligand_input: str,
        ligand_type: str = "cid"
    ) -> Dict:
        """
        Execute complete docking pipeline
        
        Args:
            pdb_id: 4-character PDB identifier
            ligand_input: PubChem CID or SMILES string
            ligand_type: Type of ligand input ("cid" or "smiles")
            
        Returns:
            Dictionary with docking results
        """
        
        try:
            print(f"\n{'='*60}")
            print(f"Starting Docking Pipeline")
            print(f"PDB ID: {pdb_id}")
            print(f"Ligand: {ligand_input} ({ligand_type})")
            print(f"{'='*60}\n")
            
            # Step 1: Fetch and prepare protein
            print("[1/5] Fetching protein structure...")
            protein_pdb = await self.protein_handler.fetch_protein(pdb_id, self.work_dir)
            
            print("[2/5] Preparing protein...")
            protein_pdbqt = await self.protein_handler.prepare_protein(protein_pdb)
            
            # Step 2: Calculate docking box
            print("[3/5] Calculating docking box dimensions...")
            center, size = self.protein_handler.get_protein_center_and_dimensions(
                protein_pdb,
                padding=self.grid_padding
            )
            
            # Step 3: Prepare ligand
            print("[4/5] Preparing ligand...")
            ligand_pdbqt = await self.ligand_handler.prepare_ligand(
                ligand_input,
                ligand_type,
                self.work_dir
            )
            
            # Get ligand properties
            # Extract PDB file path for properties calculation
            ligand_pdb = self.work_dir / "ligand.pdb"
            ligand_props = self.ligand_handler.get_ligand_properties(ligand_pdb)
            
            # Step 4: Run docking
            print("[5/5] Running AutoDock Vina...")
            docking_results = await self.vina_runner.run_docking(
                protein_pdbqt=protein_pdbqt,
                ligand_pdbqt=ligand_pdbqt,
                center=center,
                size=size,
                output_dir=self.work_dir,
                exhaustiveness=self.exhaustiveness,
                num_modes=self.num_modes
            )
            
            # Step 5: Prepare results
            print("\nPreparing results...")
            
            # Read PDB files for visualization
            protein_pdb_content = protein_pdb.read_text()
            ligand_pdb_content = ligand_pdb.read_text()
            
            # Encode for transmission
            protein_pdb_b64 = base64.b64encode(protein_pdb_content.encode()).decode()
            ligand_pdb_b64 = base64.b64encode(ligand_pdb_content.encode()).decode()
            
            # Compile final results
            results = {
                "pdb_id": pdb_id,
                "ligand_input": ligand_input,
                "ligand_type": ligand_type,
                "binding_affinity": docking_results.get("best_affinity"),
                "poses": docking_results.get("poses", []),
                "protein_pdb": protein_pdb_b64,
                "ligand_pdb": ligand_pdb_b64,
                "grid_info": {
                    "center": center,
                    "size": size,
                    "padding": self.grid_padding
                },
                "ligand_properties": ligand_props,
                "status": "completed"
            }
            
            print(f"\n{'='*60}")
            print(f"Docking Completed Successfully!")
            print(f"Best Binding Affinity: {results['binding_affinity']:.2f} kcal/mol")
            print(f"Number of poses: {len(results['poses'])}")
            print(f"{'='*60}\n")
            
            return results
            
        except Exception as e:
            print(f"\nERROR in pipeline: {str(e)}")
            raise
