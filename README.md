# Molecular Docking Dashboard

A modern, fully-automated web dashboard for molecular docking simulations. Simply input a PDB ID (protein) and a PubChem CID or SMILES string (ligand), and the system automatically handles all steps from structure preparation to visualization.

## Features

✨ **User-Friendly Interface**
- Clean, modern UI built with Tailwind CSS
- Real-time progress tracking
- Interactive 3D visualization of docking results

⚙️ **Fully Automated Pipeline**
- Automatic protein structure retrieval from RCSB PDB
- Ligand preparation from PubChem or SMILES strings
- Automatic docking box dimension calculation
- AutoDock Vina integration for molecular docking

📊 **Results & Visualization**
- Binding affinity table with detailed metrics
- Interactive 3D molecular visualization
- 2D structure images
- Download results as CSV/PDB files

## System Requirements

- Python 3.8+
- pip (Python package manager)
- ~2GB free disk space
- Linux/macOS/Windows

## Installation

### 1. Clone the Repository
```bash
git clone https://github.com/yourusername/molecular-docking-dashboard.git
cd molecular-docking-dashboard
```

### 2. Create a Virtual Environment (Recommended)
```bash
python -m venv venv

# On Windows:
venv\Scripts\activate

# On macOS/Linux:
source venv/bin/activate
```

### 3. Install Dependencies
```bash
pip install -r requirements.txt
```

### 4. Install System Dependencies

#### For Ubuntu/Debian:
```bash
sudo apt-get update
sudo apt-get install -y build-essential python3-dev autodock-vina
```

#### For macOS:
```bash
brew install autodock-vina
```

#### For Windows:
Download and install AutoDock Vina from: http://vina.scripps.edu/download.html

### 5. Download and Setup Vina
The dashboard includes automatic Vina setup. On first run, it will download the appropriate Vina binary for your system.

## Running the Application

### Start the Backend Server
```bash
python app.py
```

The server will start at `http://localhost:8000`

### Access the Dashboard
Open your browser and navigate to:
```
http://localhost:8000
```

## Usage

1. **Enter Protein Target**: Input a 4-character PDB ID (e.g., `1BVN`)
2. **Enter Ligand**: Either:
   - PubChem CID (e.g., `5280863`)
   - SMILES String (e.g., `CC(=O)Oc1ccccc1C(=O)O`)
3. **Click "Run Docking"**: The system will:
   - Download and prepare the protein
   - Prepare the ligand from PubChem or SMILES
   - Automatically calculate docking box dimensions
   - Run AutoDock Vina simulation
   - Display results and 3D visualization

4. **View Results**:
   - Binding affinity scores (kcal/mol)
   - Interactive 3D structure
   - Docking poses ranked by energy

## Project Structure

```
molecular-docking-dashboard/
├── app.py                      # FastAPI backend application
├── requirements.txt            # Python dependencies
├── index.html                  # Frontend dashboard UI
├── static/
│   ├── css/
│   │   └── styles.css         # Tailwind CSS styles
│   └── js/
│       └── main.js            # Frontend JavaScript
├── backend/
│   ├── __init__.py
│   ├── docking_pipeline.py    # Core docking pipeline
│   ├── protein_handler.py     # Protein structure handling
│   ├── ligand_handler.py      # Ligand preparation
│   └── vina_runner.py         # AutoDock Vina wrapper
├── uploads/                    # Temporary working directory
└── README.md                   # This file
```

## Key Packages

| Package | Purpose |
|---------|----------|
| `fastapi` | Web framework backend |
| `uvicorn` | ASGI server |
| `biopython` | Protein structure parsing |
| `rdkit` | Chemistry/ligand handling |
| `requests` | API calls to RCSB/PubChem |
| `numpy` | Numerical calculations |
| `pandas` | Data management |
| `meeko` | AutoDock PDBQT preparation |

## API Endpoints

### POST `/api/dock`
Runs a complete docking simulation.

**Request:**
```json
{
  "pdb_id": "1BVN",
  "ligand_input": "5280863",
  "ligand_type": "cid"
}
```

**Response:**
```json
{
  "status": "success",
  "binding_affinity": -7.3,
  "poses": [
    {
      "mode": 1,
      "affinity": -7.3,
      "rmsd_lb": 0.0,
      "rmsd_ub": 0.0
    }
  ],
  "protein_pdb": "base64_encoded_pdb",
  "ligand_pdb": "base64_encoded_pdb",
  "visualization_data": {}
}
```

## Troubleshooting

### Issue: "Vina not found"
**Solution:** 
```bash
# Manually install AutoDock Vina
pip install vina
```

### Issue: "PDB not found"
**Solution:** Ensure you're using a valid 4-character PDB ID from http://www.rcsb.org

### Issue: "Failed to prepare ligand"
**Solution:** Check that:
- PubChem CID is valid (https://pubchem.ncbi.nlm.nih.gov/)
- SMILES string is correctly formatted

### Issue: "Port 8000 already in use"
**Solution:**
```bash
python app.py --port 8001
```

## Performance Notes

- **First run**: May take 1-2 minutes for initial downloads and setup
- **Typical docking**: 30 seconds to 2 minutes depending on protein size
- **Network**: Requires internet for initial structure downloads

## Advanced Configuration

Edit `app.py` to modify:
- **Docking grid parameters**: `GRID_PADDING`, `GRID_POINTS`
- **Search exhaustiveness**: `EXHAUSTIVENESS` (default: 8)
- **CPU threads**: `NUM_MODES` (default: 9)
- **Server port**: Change in `uvicorn.run()`

## Contributing

Found a bug? Have suggestions? Open an issue or submit a pull request!

## License

MIT License - See LICENSE file for details

## Citing This Work

If you use this dashboard in your research, please cite:

```bibtex
@software{docking_dashboard2024,
  title={Molecular Docking Dashboard},
  author={Your Name},
  year={2024},
  url={https://github.com/yourusername/molecular-docking-dashboard}
}
```

Also cite the tools used:
- **AutoDock Vina**: Morris et al., 2009
- **RDKit**: Landrum et al.
- **BioPython**: Cock et al., 2009

## References

- RCSB PDB: https://www.rcsb.org/
- PubChem: https://pubchem.ncbi.nlm.nih.gov/
- AutoDock Vina: http://vina.scripps.edu/
- RDKit: https://www.rdkit.org/
- BioPython: https://biopython.org/

## Support

For issues or questions:
1. Check the Troubleshooting section
2. Review the logs in the console
3. Open a GitHub issue with:
   - Your PDB ID and ligand input
   - Error message
   - Your system (OS, Python version)

---

**Happy docking! 🧬**
