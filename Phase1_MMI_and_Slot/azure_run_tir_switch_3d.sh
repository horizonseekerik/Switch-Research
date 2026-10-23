#!/usr/bin/env bash
# ==============================================================================
# PROJECT JANUS: AZURE 4-CPU PARALLEL LAUNCHER — SWITCHED TIR OPTICAL SWITCH
# Architecture: Switched Total Internal Reflection (TIR) Deflector
# Platform: Azure Standard_D4s_v5 (4 vCPUs)
# ==============================================================================
# Usage:
#   bash azure_run_tir_switch_3d.sh              # both states, defaults
#   bash azure_run_tir_switch_3d.sh sweep        # parameter sweep (angles & widths)
#   bash azure_run_tir_switch_3d.sh amorphous    # single state test
#   bash azure_run_tir_switch_3d.sh crystalline  # single state test
# ==============================================================================

set -e

# ── 1. Auto-activate Conda environment if present ─────────────────────────────
for CONDA_PATH in "$HOME/miniforge3" "$HOME/miniconda3" "$HOME/anaconda3" "/opt/conda"; do
    if [ -f "${CONDA_PATH}/etc/profile.d/conda.sh" ]; then
        source "${CONDA_PATH}/etc/profile.d/conda.sh"
        conda activate meep 2>/dev/null || conda activate meep_env 2>/dev/null || conda activate base 2>/dev/null || true
        break
    fi
done

# ── 2. Select working Python executable with PyMEEP ───────────────────────────
PYTHON_EXEC="python3"
if ! $PYTHON_EXEC -c "import meep" >/dev/null 2>&1; then
    if python -c "import meep" >/dev/null 2>&1; then
        PYTHON_EXEC="python"
    fi
fi

# ── 3. Configure Multi-CPU Execution (MPI + OpenMP) ───────────────────────────
export OMP_NUM_THREADS=4
export OPENBLAS_NUM_THREADS=4
export MKL_NUM_THREADS=4
export NUMEXPR_NUM_THREADS=4

NUM_CPUS=4
RUN_CMD="${PYTHON_EXEC}"

# Check if mpirun works with MEEP for multi-core speedup
if command -v mpirun >/dev/null 2>&1; then
    if mpirun -np 2 "${PYTHON_EXEC}" -c "import meep" >/dev/null 2>&1; then
        echo "[INFO] OpenMPI detected! Running across ${NUM_CPUS} CPU cores with mpirun."
        RUN_CMD="mpirun -np ${NUM_CPUS} ${PYTHON_EXEC}"
    else
        echo "[INFO] Running with OpenMP multi-threading (OMP_NUM_THREADS=${NUM_CPUS})."
    fi
else
    echo "[INFO] Running with OpenMP multi-threading (OMP_NUM_THREADS=${NUM_CPUS})."
fi

# ── 4. Parameters ────────────────────────────────────────────────────────────
MODE="${1:-both}"
THETA_SLOT="${2:-30.0}"
W_SLOT="${3:-0.250}"
RES="${4:-20}"

TIMESTAMP=$(date +%Y%m%d_%H%M%S)
LOG_DIR="./logs"
mkdir -p "${LOG_DIR}"
LOG_FILE="${LOG_DIR}/tir_switch_${MODE}_${TIMESTAMP}.log"

echo "================================================================"
echo "  PROJECT JANUS: Switched TIR Optical Switch (3D MEEP FDTD)"
echo "  Principle   : Total Internal Reflection Deflector"
echo "  theta_slot  : ${THETA_SLOT} deg (theta_i = $(echo "90 - ${THETA_SLOT}" | bc -l) deg)"
echo "  w_slot      : $(echo "${W_SLOT} * 1000" | bc -l | xargs printf '%.0f') nm"
echo "  Resolution  : ${RES} px/um"
echo "  Multi-CPU   : ${RUN_CMD}"
echo "  Log File    : ${LOG_FILE}"
echo "================================================================"
echo ""

START_TIME=$(date +%s)

if [ "${MODE}" = "sweep" ]; then
    echo "Running TIR parameter sweep across angles and slot widths..."
    ${RUN_CMD} tir_switch_3d.py \
        --sweep \
        --resolution "${RES}" 2>&1 | tee "${LOG_FILE}"
elif [ "${MODE}" = "amorphous" ] || [ "${MODE}" = "crystalline" ]; then
    ${RUN_CMD} tir_switch_3d.py \
        --state "${MODE}" \
        --theta_slot "${THETA_SLOT}" \
        --w_slot "${W_SLOT}" \
        --resolution "${RES}" 2>&1 | tee "${LOG_FILE}"
else
    # both states
    ${RUN_CMD} tir_switch_3d.py \
        --state both \
        --theta_slot "${THETA_SLOT}" \
        --w_slot "${W_SLOT}" \
        --resolution "${RES}" 2>&1 | tee "${LOG_FILE}"
fi

END_TIME=$(date +%s)
ELAPSED=$((END_TIME - START_TIME))

echo ""
echo "================================================================"
echo "  COMPLETED! Elapsed: ${ELAPSED} s ($(echo "scale=2; ${ELAPSED}/60" | bc) min)"
echo "  JSON Results : ./top_tir_switch_candidates.json"
echo "  CSV Results  : ./tir_switch_results.csv"
echo "  Log File     : ${LOG_FILE}"

ARCHIVE="tir_switch_results_${MODE}_${TIMESTAMP}.tar.gz"
tar czf "${ARCHIVE}" \
    top_tir_switch_candidates.json \
    tir_switch_results.csv \
    "${LOG_FILE}" \
    2>/dev/null || \
tar czf "${ARCHIVE}" \
    tir_switch_sweep_results.csv \
    top_tir_switch_candidates.json \
    "${LOG_FILE}" \
    2>/dev/null || true

echo "  Archive Package: ./${ARCHIVE}"
echo "================================================================"
