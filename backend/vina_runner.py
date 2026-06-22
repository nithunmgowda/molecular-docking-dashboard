#!/usr/bin/env python3
"""
AutoDock Vina Runner Module
Executes molecular docking simulations
"""

import subprocess
import re
from pathlib import Path
from typing import List, Dict, Tuple, Optional
import sys


class VinaRunner:
    """Wrapper for AutoDock Vina docking engine"""
    
    def __init__(self, vina_path: str = "vina"):
        """
        Initialize Vina runner
        
        Args:
            vina_path: Path to Vina executable
        """
        self.vina_path = vina_path
        self._check_vina_available()
    
    def _check_vina_available(self) -> bool:
        """
        Check if Vina is installed and accessible
        """
        try:
            result = subprocess.run(
                [self.vina_path, "--help"],
                capture_output=True,
                text=True,
                timeout=5
            )
            return result.returncode == 0
        except (FileNotFoundError, subprocess.TimeoutExpired):
            print("Warning: AutoDock Vina not found. Trying alternative paths...")
            # Try alternative paths
            for alt_path in ["vina_1.2.5", "autodock_vina", "/usr/bin/vina"]:
                try:
                    result = subprocess.run(
                        [alt_path, "--help"],
                        capture_output=True,
                        timeout=5
                    )
                    if result.returncode == 0:
                        self.vina_path = alt_path
                        return True
                except:
                    pass
            return False
    
    async def run_docking(
        self,
        protein_pdbqt: Path,
        ligand_pdbqt: Path,
        center: Tuple[float, float, float],
        size: Tuple[float, float, float],
        output_dir: Path,
        exhaustiveness: int = 8,
        num_modes: int = 9,
        energy_range: float = 3.0
    ) -> Dict:
        """
        Run AutoDock Vina docking
        
        Args:
            protein_pdbqt: Path to protein PDBQT file
            ligand_pdbqt: Path to ligand PDBQT file
            center: Docking box center (x, y, z)
            size: Docking box size (x, y, z)
            output_dir: Directory for output files
            exhaustiveness: Search exhaustiveness (1-32)
            num_modes: Number of docking modes
            energy_range: Energy range for clustering (kcal/mol)
            
        Returns:
            Dictionary with docking results
        """
        try:
            # Output file
            output_file = output_dir / "docking_output.pdbqt"
            log_file = output_dir / "docking_log.txt"
            
            # Build Vina command
            cmd = [
                self.vina_path,
                "--receptor", str(protein_pdbqt),
                "--ligand", str(ligand_pdbqt),
                "--center_x", str(center[0]),
                "--center_y", str(center[1]),
                "--center_z", str(center[2]),
                "--size_x", str(size[0]),
                "--size_y", str(size[1]),
                "--size_z", str(size[2]),
                "--out", str(output_file),
                "--log", str(log_file),
                "--exhaustiveness", str(exhaustiveness),
                "--num_modes", str(num_modes),
                "--energy_range", str(energy_range),
                "--cpu", "1"  # Use single CPU to avoid overhead
            ]
            
            print(f"Running Vina with command: {' '.join(cmd)}")
            
            # Run Vina
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=300  # 5 minute timeout
            )
            
            if result.returncode != 0:
                error_msg = result.stderr or result.stdout
                raise RuntimeError(f"Vina failed: {error_msg}")
            
            # Parse results
            results = await self._parse_vina_output(log_file, output_file)
            
            print(f"Docking completed successfully")
            print(f"Best affinity: {results['best_affinity']:.2f} kcal/mol")
            
            return results
            
        except subprocess.TimeoutExpired:
            raise RuntimeError("Docking timed out (exceeded 5 minutes)")
        except Exception as e:
            raise RuntimeError(f"Failed to run docking: {str(e)}")
    
    async def _parse_vina_output(self, log_file: Path, pdbqt_file: Path) -> Dict:
        """
        Parse Vina output files to extract results
        
        Args:
            log_file: Path to Vina log file
            pdbqt_file: Path to output PDBQT file
            
        Returns:
            Dictionary with parsed results
        """
        try:
            results = {
                "best_affinity": None,
                "poses": [],
                "output_file": str(pdbqt_file)
            }
            
            if not log_file.exists():
                return results
            
            # Read log file
            log_content = log_file.read_text()
            
            # Extract docking results table
            # Pattern: "   1   -7.3      0.000      0.000"
            pattern = r'\s+(\d+)\s+(-\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)'
            
            for match in re.finditer(pattern, log_content):
                mode = int(match.group(1))
                affinity = float(match.group(2))
                rmsd_lb = float(match.group(3))
                rmsd_ub = float(match.group(4))
                
                pose = {
                    "mode": mode,
                    "affinity": affinity,
                    "rmsd_lb": rmsd_lb,
                    "rmsd_ub": rmsd_ub
                }
                results["poses"].append(pose)
            
            if results["poses"]:
                results["best_affinity"] = results["poses"][0]["affinity"]
            
            return results
            
        except Exception as e:
            print(f"Error parsing Vina output: {e}")
            return {"best_affinity": None, "poses": [], "output_file": str(pdbqt_file)}
