#!/usr/bin/env python3
"""
Molecular Docking Dashboard - FastAPI Backend
Automated pipeline for protein-ligand docking simulations
"""

import os
import sys
import json
import shutil
import tempfile
from pathlib import Path
from typing import Optional
from fastapi import FastAPI, HTTPException, BackgroundTasks
from fastapi.responses import FileResponse, HTMLResponse
from fastapi.staticfiles import StaticFiles
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel
import uvicorn

# Import backend modules
from backend.docking_pipeline import DockingPipeline
from backend.protein_handler import ProteinHandler
from backend.ligand_handler import LigandHandler

# ============================================================================
# Configuration
# ============================================================================

APP_NAME = "Molecular Docking Dashboard"
APP_VERSION = "1.0.0"
WORK_DIR = Path("uploads")
WORK_DIR.mkdir(exist_ok=True)

# Docking parameters
GRID_PADDING = 4.0  # Padding around protein
SEARCH_EXHAUSTIVENESS = 8  # Search thoroughness (1-32)
NUM_MODES = 9  # Number of docking modes to generate

# ============================================================================
# FastAPI Setup
# ============================================================================

app = FastAPI(
    title=APP_NAME,
    version=APP_VERSION,
    description="Automated molecular docking dashboard"
)

# Enable CORS
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# ============================================================================
# Pydantic Models
# ============================================================================

class DockingRequest(BaseModel):
    """Request model for docking simulation"""
    pdb_id: str
    ligand_input: str
    ligand_type: str = "cid"  # "cid" or "smiles"

class DockingResponse(BaseModel):
    """Response model for docking results"""
    status: str
    message: str
    binding_affinity: Optional[float] = None
    poses: Optional[list] = None
    protein_pdb: Optional[str] = None
    ligand_pdb: Optional[str] = None
    grid_info: Optional[dict] = None

class HealthResponse(BaseModel):
    """Health check response"""
    status: str
    version: str

# ============================================================================
# Routes - Health & Status
# ============================================================================

@app.get("/", response_class=HTMLResponse)
async def serve_dashboard():
    """Serve the main dashboard HTML"""
    try:
        with open("index.html", "r") as f:
            return f.read()
    except FileNotFoundError:
        return HTMLResponse(
            "<h1>Dashboard not found. Ensure index.html is in the root directory.</h1>",
            status_code=404
        )

@app.get("/health", response_model=HealthResponse)
async def health_check():
    """Health check endpoint"""
    return HealthResponse(
        status="healthy",
        version=APP_VERSION
    )

@app.get("/api/info")
async def get_info():
    """Get application information"""
    return {
        "app_name": APP_NAME,
        "version": APP_VERSION,
        "status": "running"
    }

# ============================================================================
# Routes - Docking Simulation
# ============================================================================

@app.post("/api/dock", response_model=DockingResponse)
async def run_docking(request: DockingRequest, background_tasks: BackgroundTasks):
    """
    Main docking endpoint
    
    Process:
    1. Validate inputs
    2. Download protein from RCSB PDB
    3. Prepare ligand from PubChem or SMILES
    4. Calculate docking box dimensions
    5. Run AutoDock Vina
    6. Return results and visualization data
    """
    
    work_session = WORK_DIR / f"session_{id(request)}"
    work_session.mkdir(exist_ok=True)
    
    try:
        # Validate PDB ID
        pdb_id = request.pdb_id.strip().upper()
        if len(pdb_id) != 4:
            raise ValueError("PDB ID must be exactly 4 characters")
        
        # Validate ligand input
        if request.ligand_type == "cid":
            try:
                cid = int(request.ligand_input)
            except ValueError:
                raise ValueError("PubChem CID must be a number")
        elif request.ligand_type == "smiles":
            if not request.ligand_input.strip():
                raise ValueError("SMILES string cannot be empty")
        else:
            raise ValueError("Ligand type must be 'cid' or 'smiles'")
        
        # Initialize pipeline
        pipeline = DockingPipeline(
            work_dir=work_session,
            grid_padding=GRID_PADDING,
            exhaustiveness=SEARCH_EXHAUSTIVENESS,
            num_modes=NUM_MODES
        )
        
        # Run docking pipeline
        results = await pipeline.run_full_pipeline(
            pdb_id=pdb_id,
            ligand_input=request.ligand_input,
            ligand_type=request.ligand_type
        )
        
        # Schedule cleanup
        background_tasks.add_task(cleanup_session, work_session)
        
        return DockingResponse(
            status="success",
            message="Docking completed successfully",
            binding_affinity=results.get("binding_affinity"),
            poses=results.get("poses"),
            protein_pdb=results.get("protein_pdb"),
            ligand_pdb=results.get("ligand_pdb"),
            grid_info=results.get("grid_info")
        )
        
    except Exception as e:
        print(f"Error in docking: {str(e)}", file=sys.stderr)
        raise HTTPException(status_code=400, detail=str(e))

@app.post("/api/validate-pdb")
async def validate_pdb(pdb_id: str):
    """
    Validate if PDB ID exists
    """
    try:
        handler = ProteinHandler()
        is_valid = await handler.validate_pdb_id(pdb_id.strip().upper())
        return {"valid": is_valid, "pdb_id": pdb_id.upper()}
    except Exception as e:
        raise HTTPException(status_code=400, detail=str(e))

@app.post("/api/validate-ligand")
async def validate_ligand(ligand_input: str, ligand_type: str = "cid"):
    """
    Validate ligand input (CID or SMILES)
    """
    try:
        handler = LigandHandler()
        is_valid = await handler.validate_ligand(ligand_input, ligand_type)
        return {"valid": is_valid, "ligand_input": ligand_input, "type": ligand_type}
    except Exception as e:
        raise HTTPException(status_code=400, detail=str(e))

# ============================================================================
# Utility Functions
# ============================================================================

def cleanup_session(session_path: Path):
    """Cleanup temporary session files"""
    try:
        if session_path.exists():
            shutil.rmtree(session_path)
            print(f"Cleaned up session: {session_path}")
    except Exception as e:
        print(f"Error cleaning up session: {e}", file=sys.stderr)

# ============================================================================
# Main
# ============================================================================

if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description="Molecular Docking Dashboard")
    parser.add_argument("--host", default="0.0.0.0", help="Host to bind to")
    parser.add_argument("--port", type=int, default=8000, help="Port to bind to")
    parser.add_argument("--reload", action="store_true", help="Enable auto-reload")
    
    args = parser.parse_args()
    
    print(f"\n{'='*60}")
    print(f"  {APP_NAME} v{APP_VERSION}")
    print(f"{'='*60}")
    print(f"Starting server at http://{args.host}:{args.port}")
    print(f"Dashboard: http://localhost:{args.port}")
    print(f"API Docs: http://localhost:{args.port}/docs")
    print(f"{'='*60}\n")
    
    uvicorn.run(
        "app:app",
        host=args.host,
        port=args.port,
        reload=args.reload
    )
