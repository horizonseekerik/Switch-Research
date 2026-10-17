#!/usr/bin/env bash
# ==============================================================================
# PROJECT JANUS: AZURE 4-CPU SUB-2 um^2 OPTICAL SWITCH SEARCH & 3D SIMULATION
# Target: Footprint < 2.00 um^2 | IL <= 0.15 dB (Am) | IL <= 0.20 dB (Cr) | XT <= -20.0 dB
# ==============================================================================
set -e

AZURE_CPUS=4
TIMESTAMP=$(date -u +"%Y%m%d_%H%M%S")
RESULTS_DIR="sub2um_results_${TIMESTAMP}"
mkdir -p "${RESULTS_DIR}"

echo "======================================================================"
echo "  PROJECT JANUS: AZURE 4-CPU SUB-2 um^2 OPTICAL SWITCH SEARCH"
echo "  Target: Footprint < 2.00 um^2 | IL <= 0.15 dB (Am) | IL <= 0.20 dB (Cr)"
echo "  Timestamp: $(date -u +"%Y-%m-%d %H:%M:%SZ")"
echo "======================================================================"

# Auto-activate Miniforge Conda environment if available
if [ -f "$HOME/miniforge3/etc/profile.d/conda.sh" ]; then
    source "$HOME/miniforge3/etc/profile.d/conda.sh"
    conda activate base || true
fi

echo "[INFO] Azure VM CPU count configured: ${AZURE_CPUS} cores"
python -c "import sys; print(f'Python runtime: {sys.version}')"
python -c "try:
    import meep as mp
    print(f'MEEP FDTD Status: Active (v{mp.__version__}) - Full-Wave 3D Maxwell Mode')
except ImportError:
    print('MEEP FDTD Status: Analytical CMT Mode (Run ./install_meep_azure.sh for full wave)')"


# Phase 1: Primary 3D Maxwell & Modal CMT Simulation
echo ""
echo "======================================================================"
echo "  [PRIMARY] EXECUTING 3D VECTOR MAXWELL SWITCH SIMULATION (SUB-2 um^2)"
echo "  Target: 1.620 um^2 Cell (w=240nm, gap=18nm, L_c=2.559um, L_bend=0.38um)"
echo "======================================================================"
python sub2um_switch_3d_meep.py --state both --output-dir "${RESULTS_DIR}"
cp "${RESULTS_DIR}/sub2um_3d_switch_results.json" ./sub2um_3d_switch_results.json

# Phase 2: Combinatorial Parameter Grid Search (7,200 designs)
echo ""
echo "======================================================================"
echo "  [SECONDARY] RUNNING 4-CPU SUB-2 um^2 PARAMETER GRID SEARCH"
echo "======================================================================"
python sub2um_switch_search.py --cpus "${AZURE_CPUS}" --output-dir "${RESULTS_DIR}"
cp "${RESULTS_DIR}/sub2um_top_candidates.json" ./sub2um_top_candidates.json
cp "${RESULTS_DIR}/sub2um_switch_search_results.csv" ./sub2um_switch_search_results.csv

# Phase 3: Archive Results
ARCHIVE_NAME="sub2um_switch_results_${TIMESTAMP}.tar.gz"
tar -czf "${ARCHIVE_NAME}" "${RESULTS_DIR}" sub2um_3d_switch_results.json sub2um_top_candidates.json sub2um_switch_search_results.csv
echo ""
echo "======================================================================"
echo "  [SUCCESS] All Sub-2 um^2 simulations and search runs completed!"
echo "  Archive Package: ${ARCHIVE_NAME}"
echo "======================================================================"
