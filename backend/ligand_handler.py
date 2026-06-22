#!/usr/bin/env python3
"""
Ligand Handler Module
Handles ligand structure retrieval, preparation, and conversion to 3D
"""

import subprocess
from pathlib import Path
from typing import Optional
import requests
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors
from rdkit.Chem import rdMolDescriptors


class LigandHandler:
    """Handles ligand structure operations"""
    
    PUBCHEM_API_URL = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
    
    def __init__(self):
        pass
    
    async def validate_ligand(self, ligand_input: str, ligand_type: str = "cid") -> bool:
        """
        Validate ligand input
        
        Args:
            ligand_input: PubChem CID or SMILES string
            ligand_type: Type of input ("cid" or "smiles")
            
        Returns:
            True if valid, False otherwise
        """
        try:
            if ligand_type == "cid":
                cid = int(ligand_input)
                # Try to fetch
                url = f"{self.PUBCHEM_API_URL}/compound/cid/{cid}/property/MolecularFormula/json"
                response = requests.get(url, timeout=5)
                return response.status_code == 200
            
            elif ligand_type == "smiles":
                mol = Chem.MolFromSmiles(ligand_input)
                return mol is not None
            
            return False
            
        except Exception as e:
            print(f"Error validating ligand: {e}")
            return False
    
    async def fetch_ligand_from_pubchem(self, cid: str, work_dir: Path) -> Path:
        """
        Download ligand structure from PubChem
        
        Args:
            cid: PubChem Compound ID
            work_dir: Working directory
            
        Returns:
            Path to SDF file
        """
        try:
            cid = str(cid).strip()
            sdf_file = work_dir / f"ligand_{cid}.sdf"
            
            if sdf_file.exists():
                print(f"Using cached ligand: {sdf_file}")
                return sdf_file
            
            # Download SDF from PubChem
            url = f"{self.PUBCHEM_API_URL}/compound/cid/{cid}/record/SDF"
            print(f"Downloading ligand from: {url}")
            
            response = requests.get(url, timeout=30)
            response.raise_for_status()
            
            sdf_file.write_text(response.text)
            print(f"Downloaded ligand to: {sdf_file}")
            
            return sdf_file
            
        except Exception as e:
            raise RuntimeError(f"Failed to fetch ligand from PubChem: {str(e)}")
    
    async def smiles_to_structure(self, smiles: str, work_dir: Path) -> Path:
        """
        Convert SMILES string to 3D SDF structure
        
        Args:
            smiles: SMILES string
            work_dir: Working directory
            
        Returns:
            Path to SDF file
        """
        try:
            # Create molecule from SMILES
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                raise ValueError(f"Invalid SMILES: {smiles}")
            
            # Add hydrogens
            mol = Chem.AddHs(mol)
            
            # Generate 3D coordinates
            AllChem.EmbedMolecule(mol, randomSeed=42)
            AllChem.UMMFFOptimizeMolecule(mol)
            
            # Save as SDF
            sdf_file = work_dir / "ligand_from_smiles.sdf"
            writer = Chem.SDWriter(str(sdf_file))
            writer.write(mol)
            writer.close()
            
            print(f"Generated ligand from SMILES: {sdf_file}")
            return sdf_file
            
        except Exception as e:
            raise RuntimeError(f"Failed to convert SMILES to structure: {str(e)}")
    
    async def prepare_ligand(
        self,
        ligand_input: str,
        ligand_type: str,
        work_dir: Path
    ) -> Path:
        """
        Prepare ligand for docking (generate 3D, convert to PDBQT)
        
        Args:
            ligand_input: PubChem CID or SMILES string
            ligand_type: Type of input ("cid" or "smiles")
            work_dir: Working directory
            
        Returns:
            Path to PDBQT file
        """
        try:
            # Get SDF file
            if ligand_type == "cid":
                sdf_file = await self.fetch_ligand_from_pubchem(ligand_input, work_dir)
            elif ligand_type == "smiles":
                sdf_file = await self.smiles_to_structure(ligand_input, work_dir)
            else:
                raise ValueError(f"Unknown ligand type: {ligand_type}")
            
            # Convert SDF to PDB
            pdb_file = work_dir / "ligand.pdb"
            await self._sdf_to_pdb(sdf_file, pdb_file)
            
            # Convert PDB to PDBQT
            pdbqt_file = work_dir / "ligand.pdbqt"
            await self._pdb_to_pdbqt(pdb_file, pdbqt_file)
            
            print(f"Prepared ligand: {pdbqt_file}")
            return pdbqt_file
            
        except Exception as e:
            raise RuntimeError(f"Failed to prepare ligand: {str(e)}")
    
    async def _sdf_to_pdb(self, sdf_file: Path, pdb_file: Path):
        """
        Convert SDF to PDB format
        """
        try:
            # Read SDF
            suppl = Chem.SDMolSupplier(str(sdf_file))
            if len(suppl) == 0:
                raise ValueError("No molecules in SDF file")
            
            mol = suppl[0]
            if mol is None:
                raise ValueError("Failed to read molecule from SDF")
            
            # Ensure 3D coordinates
            if mol.GetNumConformers() == 0:
                mol = Chem.AddHs(mol)
                AllChem.EmbedMolecule(mol, randomSeed=42)
                AllChem.UMMFFOptimizeMolecule(mol)
            
            # Convert to PDB
            pdb_block = Chem.MolToPDBBlock(mol)
            pdb_file.write_text(pdb_block)
            
        except Exception as e:
            raise RuntimeError(f"Failed to convert SDF to PDB: {str(e)}")
    
    async def _pdb_to_pdbqt(self, pdb_file: Path, pdbqt_file: Path):
        """
        Convert PDB to PDBQT format
        """
        try:
            # Try using meeko
            try:
                cmd = f"mk_prepare_ligand.py -i {pdb_file} -o {pdbqt_file}"
                result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
                if result.returncode == 0:
                    return
            except:
                pass
            
            # Fallback: simple conversion
            lines = pdb_file.read_text().split('\n')
            pdbqt_lines = []
            
            atom_count = 0
            for line in lines:
                if line.startswith('ATOM') or line.startswith('HETATM'):
                    atom_count += 1
                    # Ensure proper formatting
                    if len(line) < 66:
                        line = line + (' ' * (66 - len(line)))
                    
                    # Add charge (simplified)
                    pdbqt_line = line[:66] + '  0.00  0.00'
                    pdbqt_lines.append(pdbqt_line)
                else:
                    pdbqt_lines.append(line)
            
            pdbqt_file.write_text('\n'.join(pdbqt_lines))
            
        except Exception as e:
            raise RuntimeError(f"Failed to convert PDB to PDBQT: {str(e)}")
    
    def get_ligand_properties(self, pdb_file: Path) -> dict:
        """
        Calculate ligand molecular properties
        
        Args:
            pdb_file: Path to PDB file
            
        Returns:
            Dictionary with molecular properties
        """
        try:
            # Read PDB and convert to mol
            mol = Chem.MolFromPDBFile(str(pdb_file))
            if mol is None:
                return {}
            
            properties = {
                "molecular_weight": round(Descriptors.MolWt(mol), 2),
                "logp": round(Descriptors.MolLogP(mol), 2),
                "hbd": Descriptors.NumHDonors(mol),  # H-bond donors
                "hba": Descriptors.NumHAcceptors(mol),  # H-bond acceptors
                "rotatable_bonds": Descriptors.NumRotatableBonds(mol),
                "molecular_formula": rdMolDescriptors.CalcMolFormula(mol)
            }
            
            return properties
            
        except Exception as e:
            print(f"Error calculating ligand properties: {e}")
            return {}
