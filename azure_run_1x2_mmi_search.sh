#!/usr/bin/env bash
# ==============================================================================
# PROJECT JANUS: AZURE MULTI-CORE SINGLE-CAVITY 1x2 MMI SWITCH SEARCH LAUNCHER
# Target: Single-Cavity 1x2 MMI | IL <= 0.08 dB | XT <= -20 dB | Footprint < 5.0 um^2
# Platform: Free Tier / Cloud Compute (Azure 4-vCPU or HPC Node)
# ==============================================================================
set -e

TIMESTAMP=$(date -u +"%Y%m%d_%H%M%S")
RESULTS_DIR="mmi_1x2_results_${TIMESTAMP}"
mkdir -p "${RESULTS_DIR}"

# Auto-detect CPU cores (default to 4 if undetected)
if command -v nproc > /dev/null 2>&1; then
    CPU_CORES=$(nproc)
else
    CPU_CORES=4
fi

echo "======================================================================"
echo "  PROJECT JANUS: 1x2 SINGLE-CAVITY MMI SWITCH CLOUD SEARCH"
echo "  Target: Single Multimode Body + Asymmetric Sb2S3 Patch"
echo "  Optical Targets: IL <= 0.08 dB | Crosstalk <= -20 dB | Footprint < 5.0 um^2"
echo "  Timestamp: $(date -u +"%Y-%m-%d %H:%M:%SZ")"
echo "  Detected CPU Cores: ${CPU_CORES}"
echo "======================================================================"

# 1. Activate conda environment if available (Miniforge, Miniconda, Anaconda)
for CONDA_PATH in "$HOME/miniforge3" "$HOME/miniconda3" "$HOME/anaconda3" "/opt/conda"; do
    if [ -f "${CONDA_PATH}/etc/profile.d/conda.sh" ]; then
        source "${CONDA_PATH}/etc/profile.d/conda.sh"
        conda activate meep 2>/dev/null || conda activate base 2>/dev/null || true
        break
    fi
done

# 2. Check Python & MEEP runtime status
python3 -c "import sys; print(f'Python runtime: {sys.version.split()[0]}')"
python3 -c "try:
    import meep as mp
    print(f'MEEP FDTD Status: Active (v{mp.__version__}) - Full-Wave 3D Maxwell Mode')
except ImportError:
    print('WARNING: MEEP not found in active environment. Run ./install_meep_azure.sh first.')"

# 3. Launch the Master 1x2 MMI Cloud Search Script
echo ""
echo "======================================================================"
echo "  LAUNCHING MASTER PARALLEL SEARCH ACROSS ${CPU_CORES} CPU CORES..."
echo "======================================================================"

python3 master_1x2_mmi_cloud_search.py \
    --cpus "${CPU_CORES}" \
    --resolution 20 \
    --mode both \
    --output-dir "${RESULTS_DIR}"

# 4. Copy primary output files to repository root for immediate inspection
cp "${RESULTS_DIR}/top_1x2_mmi_candidates.json" ./top_1x2_mmi_candidates.json || true
cp "${RESULTS_DIR}/single_cavity_1x2_mmi_results.csv" ./single_cavity_1x2_mmi_results.csv || true

# 5. Create timestamped archive package
ARCHIVE_NAME="mmi_1x2_search_results_${TIMESTAMP}.tar.gz"
tar -czf "${ARCHIVE_NAME}" "${RESULTS_DIR}" top_1x2_mmi_candidates.json single_cavity_1x2_mmi_results.csv 2>/dev/null || true

echo ""
echo "======================================================================"
echo "  [SUCCESS] 1x2 MMI Switch Search & 3D FDTD Verification Completed!"
echo "  Primary Results : ./top_1x2_mmi_candidates.json"
echo "  Tabular CSV     : ./single_cavity_1x2_mmi_results.csv"
echo "  Archive Package : ${ARCHIVE_NAME}"
echo "======================================================================"
