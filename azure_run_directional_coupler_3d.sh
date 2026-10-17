#!/bin/bash
# ================================================================
# PROJECT JANUS: Azure VM Launcher - 1x2 Directional Coupler Switch
# ================================================================
# Usage:
#   bash azure_run_directional_coupler_3d.sh              # both states, defaults
#   bash azure_run_directional_coupler_3d.sh sweep        # parameter sweep
#   bash azure_run_directional_coupler_3d.sh amorphous    # single state
#   bash azure_run_directional_coupler_3d.sh crystalline  # single state
#
# Requirements: conda activate meep_env (or equivalent)
# Platform:     Azure Standard_D4s_v5 (4 vCPUs)
# ================================================================

set -e

# ── Auto-activate Conda environment if present ───────────────────
for CONDA_PATH in "$HOME/miniforge3" "$HOME/miniconda3" "$HOME/anaconda3" "/opt/conda"; do
    if [ -f "${CONDA_PATH}/etc/profile.d/conda.sh" ]; then
        source "${CONDA_PATH}/etc/profile.d/conda.sh"
        conda activate meep 2>/dev/null || conda activate meep_env 2>/dev/null || conda activate base 2>/dev/null || true
        break
    fi
done

PYTHON_EXEC="python3"
if ! $PYTHON_EXEC -c "import meep" >/dev/null 2>&1; then
    if python -c "import meep" >/dev/null 2>&1; then
        PYTHON_EXEC="python"
    fi
fi

# ── Multi-core CPU settings ───────────────────────────────────────
export OMP_NUM_THREADS=4
export OPENBLAS_NUM_THREADS=4
export MKL_NUM_THREADS=4
export NUMEXPR_NUM_THREADS=4

# ── Parameters ───────────────────────────────────────────────────
STATE="${1:-both}"
LC="${2:-3.55}"
GAP="${3:-0.095}"
RES="${4:-20}"

TIMESTAMP=$(date +%Y%m%d_%H%M%S)
LOG_DIR="./logs"
mkdir -p "${LOG_DIR}"
LOG_FILE="${LOG_DIR}/dc_run_${STATE}_${TIMESTAMP}.log"

echo "================================================================"
echo "  PROJECT JANUS: 1x2 Directional Coupler Switch (3D MEEP FDTD)"
echo "================================================================"
echo "  State      : ${STATE}"
echo "  L_c        : ${LC} um"
echo "  Gap        : $(echo "${GAP} * 1000" | bc -l | xargs printf '%.0f') nm"
echo "  Resolution : ${RES} px/um"
echo "  OMP Threads: ${OMP_NUM_THREADS}"
echo "  Log        : ${LOG_FILE}"
echo "================================================================"
echo ""

START_TIME=$(date +%s)

if [ "${STATE}" = "sweep" ]; then
    echo "Running parameter sweep..."
    "${PYTHON_EXEC}" master_directional_coupler_3d.py \
        --sweep \
        --lc_min 3.30 --lc_max 3.80 --lc_step 0.10 \
        --gaps 0.090 0.095 \
        --resolution "${RES}" 2>&1 | tee "${LOG_FILE}"
else
    "${PYTHON_EXEC}" master_directional_coupler_3d.py \
        --state "${STATE}" \
        --lc "${LC}" \
        --gap "${GAP}" \
        --resolution "${RES}" 2>&1 | tee "${LOG_FILE}"
fi

END_TIME=$(date +%s)
ELAPSED=$((END_TIME - START_TIME))

echo ""
echo "================================================================"
echo "  COMPLETED! Elapsed: ${ELAPSED} s ($(echo "scale=2; ${ELAPSED}/60" | bc) min)"
echo "  Results  : ./top_directional_coupler_candidates.json"
echo "  CSV      : ./directional_coupler_results.csv"
echo "  Log      : ${LOG_FILE}"

# Package results
ARCHIVE="dc_results_${STATE}_${TIMESTAMP}.tar.gz"
tar czf "${ARCHIVE}" \
    top_directional_coupler_candidates.json \
    directional_coupler_results.csv \
    "${LOG_FILE}" \
    2>/dev/null || true
echo "  Archive  : ./${ARCHIVE}"
echo "================================================================"
