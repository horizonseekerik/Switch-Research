#!/usr/bin/env bash
# ==============================================================================
# PROJECT JANUS: AZURE 4-CPU PARALLEL SWITCH SEARCH LAUNCHER
# Target: Sub-5 um^2 Optical Switch Exploration on Azure Cloud VM
# ==============================================================================
set -euo pipefail

echo "======================================================================"
echo "  PROJECT JANUS: AZURE 4-CPU SUB-5 um^2 OPTICAL SWITCH SEARCH"
echo "  Target: Footprint < 5.0 um^2 | IL <= 0.15 dB (Am) | IL <= 0.20 dB (Cr)"
echo "  Timestamp: $(date -u +"%Y-%m-%d %H:%M:%SZ")"
echo "======================================================================"

# Determine script directory
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WORKSPACE_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

# Configure CPU affinity and threads for 4-core Azure VM
export OMP_NUM_THREADS=4
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export PYTHONUNBUFFERED=1

NUM_CPUS=4
echo "[INFO] Azure VM CPU count configured: ${NUM_CPUS} cores"

# Activate Python environment
if [ -n "${VIRTUAL_ENV:-}" ]; then
    echo "[INFO] Using active virtualenv at ${VIRTUAL_ENV}..."
elif [ -d "/opt/janus_env" ]; then
    echo "[INFO] Activating production Azure environment at /opt/janus_env..."
    # shellcheck disable=SC1091
    source /opt/janus_env/bin/activate
elif [ -d "${SCRIPT_DIR}/.venv" ]; then
    echo "[INFO] Activating repo virtualenv at ${SCRIPT_DIR}/.venv..."
    # shellcheck disable=SC1091
    source "${SCRIPT_DIR}/.venv/bin/activate"
elif [ -d "${WORKSPACE_ROOT}/.venv" ]; then
    echo "[INFO] Activating workspace virtualenv at ${WORKSPACE_ROOT}/.venv..."
    # shellcheck disable=SC1091
    source "${WORKSPACE_ROOT}/.venv/bin/activate"
elif [ -d "${HOME}/janus_env" ]; then
    echo "[INFO] Activating user environment at ${HOME}/janus_env..."
    # shellcheck disable=SC1091
    source "${HOME}/janus_env/bin/activate"
elif [ -d "${HOME}/.venv" ]; then
    echo "[INFO] Activating user virtualenv at ${HOME}/.venv..."
    # shellcheck disable=SC1091
    source "${HOME}/.venv/bin/activate"
fi

# Run pre-flight check
echo "[INFO] Running Python environment verification..."
python3 -c "import sys; print(f'Python runtime: {sys.version}')"
python3 -c "try:
    import meep as mp
    print('MEEP FDTD Status: Available (Full Maxwell Mode)')
except ImportError:
    print('MEEP FDTD Status: Analytical CMT Mode (Install conda-forge/pymeep for full wave)')"

# ==============================================================================
# PHASE 1 (PRIMARY & MANDATORY): FULL-WAVE 3D SWITCH SIMULATION
# ==============================================================================
echo ""
echo "======================================================================"
echo "  [PRIMARY] EXECUTING 3D VECTOR MAXWELL SWITCH SIMULATION"
echo "  Target: 3D Silicon Wire (220nm) + Sb2S3 Thin Film (20nm) + Hermite S-Bends"
echo "======================================================================"
python3 "${SCRIPT_DIR}/sub5um_switch_3d_meep.py" --state both --resolution 20 --output-dir "${SCRIPT_DIR}"

# ==============================================================================
# PHASE 2 (OPTIONAL / SWEEP): 4-CPU PARAMETER GRID SEARCH
# ==============================================================================
echo ""
echo "======================================================================"
echo "  [OPTIONAL] RUNNING 4-CPU 2D PARAMETER GRID SEARCH (5,040 DESIGNS)"
echo "======================================================================"
python3 "${SCRIPT_DIR}/azure_4cpu_grid_search.py" --cpus "${NUM_CPUS}" --output-dir "${SCRIPT_DIR}" "$@"

# ==============================================================================
# PHASE 3: COMPREHENSIVE VERIFICATION SUITE (3D & 2D TESTS)
# ==============================================================================
echo ""
echo "======================================================================"
echo "  PHASE 3: RUNNING VERIFICATION TEST SUITE"
echo "======================================================================"
if python3 -c "import pytest" 2>/dev/null; then
    python3 -m pytest "${SCRIPT_DIR}/test_sub5um_switch.py" -v
else
    echo "[INFO] pytest not installed, executing suite via standalone test runner..."
    python3 "${SCRIPT_DIR}/test_sub5um_switch.py"
fi

# ==============================================================================
# ARCHIVE RESULTS
# ==============================================================================
TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
ARCHIVE_NAME="sub5um_switch_results_${TIMESTAMP}.tar.gz"
echo ""
echo "[INFO] Archiving results to ${ARCHIVE_NAME}..."
tar -czf "${SCRIPT_DIR}/${ARCHIVE_NAME}" \
    -C "${SCRIPT_DIR}" \
    sub5um_3d_switch_results.json \
    sub5um_switch_search_results.csv \
    sub5um_top_candidates.json

echo ""
echo "======================================================================"
echo "  SUCCESS: 3D Switch Simulation & Verification Completed Successfully!"
echo "  3D Model Results : ${SCRIPT_DIR}/sub5um_3d_switch_results.json"
echo "  Parameter Search : ${SCRIPT_DIR}/sub5um_switch_search_results.csv"
echo "  Top Candidates   : ${SCRIPT_DIR}/sub5um_top_candidates.json"
echo "  Archive Package  : ${SCRIPT_DIR}/${ARCHIVE_NAME}"
echo "======================================================================"
