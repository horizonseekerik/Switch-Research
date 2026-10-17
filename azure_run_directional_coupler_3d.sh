#!/usr/bin/env bash
# ==============================================================================
# PROJECT JANUS: AZURE 4-CPU 3D DIRECTIONAL COUPLER SWITCH SIMULATION LAUNCHER
# Platform: Standard SOI (220 nm Si, SiO2 Clad) | Wavelength: 1064 nm
# Targets: Footprint < 5.00 um^2 | IL <= 0.08 dB | XT <= -20.0 dB (Both States)
# ==============================================================================
set -e

AZURE_CPUS=4
TIMESTAMP=$(date -u +"%Y%m%d_%H%M%S")
RESULTS_DIR="dc_results_${TIMESTAMP}"
mkdir -p "${RESULTS_DIR}"

echo "======================================================================"
echo "  PROJECT JANUS: AZURE 4-CPU 3D DIRECTIONAL COUPLER SWITCH EXECUTION"
echo "  Targets: Footprint < 5.00 um^2 | IL <= 0.08 dB | XT <= -20.0 dB"
echo "  Timestamp: $(date -u +"%Y-%m-%d %H:%M:%SZ")"
echo "======================================================================"

# Auto-activate Miniforge Conda environment if available
if [ -f "$HOME/miniforge3/etc/profile.d/conda.sh" ]; then
    source "$HOME/miniforge3/etc/profile.d/conda.sh"
    conda activate base || true
fi

# Configure multi-core threading for MEEP 3D FDTD grid acceleration
export OMP_NUM_THREADS="${AZURE_CPUS}"
export MKL_NUM_THREADS="${AZURE_CPUS}"
export NUMEXPR_NUM_THREADS="${AZURE_CPUS}"
export OPENBLAS_NUM_THREADS="${AZURE_CPUS}"

echo "[INFO] Azure VM CPU count configured: ${AZURE_CPUS} cores"
python3 -c "import sys; print(f'Python runtime: {sys.version}')"
python3 -c "import meep as mp; print(f'MEEP 3D Engine: Active (v{mp.__version__}) - 4-Core Accelerated')"

# Clean previous local outputs
rm -f ./directional_coupler_3d_results.json ./directional_coupler_results.csv

echo ""
echo "======================================================================"
echo "  [STARTING] FULL-WAVE 3D VECTOR MAXWELL FDTD (BOTH STATES)"
echo "  Architecture: Balanced Symmetrical Sb2S3 Directional Coupler"
echo "======================================================================"
python3 master_directional_coupler_3d.py --state both --resolution 20 --output-dir "${RESULTS_DIR}"

# Copy results to current working directory for immediate inspection
cp "${RESULTS_DIR}/directional_coupler_3d_results.json" ./directional_coupler_3d_results.json
cp "${RESULTS_DIR}/directional_coupler_results.csv" ./directional_coupler_results.csv

echo ""
echo "======================================================================"
echo "  SIMULATION RESULTS SUMMARY (3D VECTOR MAXWELL FDTD):"
echo "======================================================================"
if command -v column &> /dev/null; then
    column -s, -t < directional_coupler_results.csv
else
    cat directional_coupler_results.csv
fi

# Archive Package
ARCHIVE_NAME="directional_coupler_results_${TIMESTAMP}.tar.gz"
tar -czf "${ARCHIVE_NAME}" "${RESULTS_DIR}" directional_coupler_3d_results.json directional_coupler_results.csv

echo ""
echo "======================================================================"
echo "  [SUCCESS] 3D Directional Coupler Switch Simulation Completed!"
echo "  Primary Results : ./directional_coupler_3d_results.json"
echo "  Tabular CSV     : ./directional_coupler_results.csv"
echo "  Archive Package : ${ARCHIVE_NAME}"
echo "======================================================================"
