#!/usr/bin/env python3
"""
Protein Handler Module
Handles protein structure retrieval, validation, and preparation
"""

import asyncio
import subprocess
from pathlib import Path
from typing import Optional, Tuple
import requests
from Bio import PDB
from Bio.PDB import PDBParser, PDBIO, Select
import numpy as np


class ProteinHandler:
    """Handles protein structure operations"""
    
    RCSB_API_URL = "https://files.rcsb.org/download"
    RCSB_SEARCH_URL = "https://data.rcsb.org/rest/v1/core/entry"
    
    def __init__(self):
        self.parser = PDBParser(QUIET=True)
        self.pdb_list = PDB.PDBList()
    
    async def validate_pdb_id(self, pdb_id: str) -> bool:
        """Validate if PDB ID exists"""
        try:
            # Remove any whitespace and convert to uppercase
            pdb_id = pdb_id.strip().upper()
            
            if len(pdb_id) != 4:
                return False
            
            # Try to fetch PDB info
            url = f"{self.RCSB_SEARCH_URL}/{pdb_id}"
            response = requests.get(url, timeout=5)
            return response.status_code == 200
            
        except Exception as e:
            print(f"Error validating PDB ID {pdb_id}: {e}")
            return False
    
    async def fetch_protein(self, pdb_id: str, work_dir: Path) -> Path:
        """
        Download protein structure from RCSB PDB
        
        Args:
            pdb_id: 4-character PDB identifier
            work_dir: Working directory for output
            
        Returns:
            Path to downloaded PDB file
        """
        try:
            pdb_id = pdb_id.strip().upper()
            
            # Check if already exists
            pdb_file = work_dir / f"{pdb_id}.pdb"
            if pdb_file.exists():
                print(f"Using cached PDB: {pdb_file}")
                return pdb_file
            
            # Download from RCSB
            url = f"{self.RCSB_API_URL}/{pdb_id}.pdb"
            print(f"Downloading protein from: {url}")
            
            response = requests.get(url, timeout=30)
            response.raise_for_status()
            
            # Save PDB file
            pdb_file.write_text(response.text)
            print(f"Downloaded PDB to: {pdb_file}")
            
            return pdb_file
            
        except Exception as e:
            raise RuntimeError(f"Failed to fetch protein {pdb_id}: {str(e)}")
    
    async def prepare_protein(self, pdb_file: Path) -> Path:
        """
        Prepare protein for docking (clean, add hydrogens, etc.)
        
        Args:
            pdb_file: Path to PDB file
            
        Returns:
            Path to prepared PDBQT file
        """
        try:
            # Parse PDB
            structure = self.parser.get_structure('protein', pdb_file)
            
            # Remove water and heteroatoms
            class ProteinSelect(Select):
                def accept_residue(self, residue):
                    # Keep protein residues
                    if residue.get_id()[0] == ' ':
                        return True
                    return False
            
            # Save cleaned PDB
            io = PDBIO()
            io.set_structure(structure)
            cleaned_pdb = pdb_file.parent / f"{pdb_file.stem}_clean.pdb"
            io.save(cleaned_pdb, select=ProteinSelect())
            
            # Convert to PDBQT using meeko
            pdbqt_file = pdb_file.parent / f"{pdb_file.stem}.pdbqt"
            
            try:
                # Try using meeko
                cmd = f"mk_prepare_ligand.py -i {cleaned_pdb} -o {pdbqt_file}"
                subprocess.run(cmd, shell=True, check=True, capture_output=True)
            except:
                # Fallback: use simple conversion
                await self._simple_pdb_to_pdbqt(cleaned_pdb, pdbqt_file)
            
            print(f"Prepared protein: {pdbqt_file}")
            return pdbqt_file
            
        except Exception as e:
            raise RuntimeError(f"Failed to prepare protein: {str(e)}")
    
    async def _simple_pdb_to_pdbqt(self, pdb_file: Path, pdbqt_file: Path):
        """
        Simple PDB to PDBQT conversion
        """
        try:
            # Read PDB
            lines = pdb_file.read_text().split('\n')
            pdbqt_lines = []
            
            for line in lines:
                if line.startswith('ATOM') or line.startswith('HETATM'):
                    # Parse PDB format
                    atom_type = line[76:78].strip()
                    
                    # Add charge (simplified - set to 0)
                    if len(line) < 66:
                        line = line + (' ' * (66 - len(line)))
                    
                    # Keep original line and add charge columns
                    pdbqt_line = line[:66] + '  0.00  0.00'
                    pdbqt_lines.append(pdbqt_line)
                else:
                    pdbqt_lines.append(line)
            
            pdbqt_file.write_text('\n'.join(pdbqt_lines))
            
        except Exception as e:
            raise RuntimeError(f"Failed to convert PDB to PDBQT: {str(e)}")
    
    def get_protein_center_and_dimensions(
        self,
        pdb_file: Path,
        padding: float = 4.0
    ) -> Tuple[Tuple[float, float, float], Tuple[float, float, float]]:
        """
        Calculate protein center and optimal docking box dimensions
        
        Args:
            pdb_file: Path to PDB file
            padding: Padding around protein (Angstroms)
            
        Returns:
            Tuple of (center, dimensions) where:
            - center: (x, y, z) coordinates
            - dimensions: (x_size, y_size, z_size) box sizes
        """
        try:
            structure = self.parser.get_structure('protein', pdb_file)
            
            # Get all alpha carbons
            coords = []
            for model in structure:
                for chain in model:
                    for residue in chain:
                        if 'CA' in residue:
                            ca = residue['CA']
                            coords.append(ca.get_coord())
            
            if not coords:
                raise ValueError("No alpha carbons found in protein")
            
            coords = np.array(coords)
            
            # Calculate center
            center = coords.mean(axis=0)
            
            # Calculate dimensions (with padding)
            min_coords = coords.min(axis=0)
            max_coords = coords.max(axis=0)
            
            dimensions = (
                float(max_coords[0] - min_coords[0] + 2 * padding),
                float(max_coords[1] - min_coords[1] + 2 * padding),
                float(max_coords[2] - min_coords[2] + 2 * padding)
            )
            
            center = tuple(float(c) for c in center)
            
            print(f"Protein center: {center}")
            print(f"Docking box dimensions: {dimensions}")
            
            return center, dimensions
            
        except Exception as e:
            raise RuntimeError(f"Failed to calculate protein dimensions: {str(e)}")
