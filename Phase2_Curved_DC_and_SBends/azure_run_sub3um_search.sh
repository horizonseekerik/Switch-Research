#!/usr/bin/env bash
# ==============================================================================
# PROJECT JANUS: AZURE 4-CPU SUB-3 um^2 OPTICAL SWITCH SEARCH & 3D SIMULATION
# Target: Footprint < 3.00 um^2 | IL <= 0.15 dB (Am) | IL <= 0.20 dB (Cr) | XT <= -20.0 dB
# ==============================================================================
set -e

AZURE_CPUS=4
TIMESTAMP=$(date -u +"%Y%m%d_%H%M%S")
RESULTS_DIR="sub3um_results_${TIMESTAMP}"
mkdir -p "${RESULTS_DIR}"

echo "======================================================================"
echo "  PROJECT JANUS: AZURE 4-CPU SUB-3 um^2 OPTICAL SWITCH SEARCH"
echo "  Target: Footprint < 3.00 um^2 | IL <= 0.15 dB (Am) | IL <= 0.20 dB (Cr)"
echo "  Timestamp: $(date -u +"%Y-%m-%d %H:%M:%SZ")"
echo "======================================================================"

echo "[INFO] Azure VM CPU count configured: ${AZURE_CPUS} cores"
python3 -c "import sys; print(f'Python runtime: {sys.version}')"

# Phase 1: Primary 3D Maxwell & Modal CMT Simulation
echo ""
echo "======================================================================"
echo "  [PRIMARY] EXECUTING 3D VECTOR MAXWELL SWITCH SIMULATION (SUB-3 um^2)"
echo "  Target: 2.368 um^2 Cell (w=240nm, gap=25nm, L_c=2.638um, L_bend=0.55um)"
echo "======================================================================"
python3 sub3um_switch_3d_meep.py --state both --output-dir "${RESULTS_DIR}"
cp "${RESULTS_DIR}/sub3um_3d_switch_results.json" ./sub3um_3d_switch_results.json

# Phase 2: Combinatorial Parameter Grid Search (8,000 designs)
echo ""
echo "======================================================================"
echo "  [SECONDARY] RUNNING 4-CPU SUB-3 um^2 PARAMETER GRID SEARCH"
echo "======================================================================"
python3 sub3um_switch_search.py --cpus "${AZURE_CPUS}" --output-dir "${RESULTS_DIR}"
cp "${RESULTS_DIR}/sub3um_top_candidates.json" ./sub3um_top_candidates.json
cp "${RESULTS_DIR}/sub3um_switch_search_results.csv" ./sub3um_switch_search_results.csv

# Phase 3: Archive Results
ARCHIVE_NAME="sub3um_switch_results_${TIMESTAMP}.tar.gz"
tar -czf "${ARCHIVE_NAME}" "${RESULTS_DIR}" sub3um_3d_switch_results.json sub3um_top_candidates.json sub3um_switch_search_results.csv
echo ""
echo "======================================================================"
echo "  [SUCCESS] All Sub-3 um^2 simulations and search runs completed!"
echo "  Archive Package: ${ARCHIVE_NAME}"
echo "======================================================================"
