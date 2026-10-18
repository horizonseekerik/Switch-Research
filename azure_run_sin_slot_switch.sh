#!/usr/bin/env bash
# ==============================================================================
# PROJECT JANUS: AZURE VM LAUNCHER — Si3N4 EMBEDDED PCM SLOT SWITCH (GEOMETRY B)
# Platform: Azure Standard_D4s_v5 (4 vCPUs)
# ==============================================================================
# Usage:
#   bash azure_run_sin_slot_switch.sh              # Dual-state verification (L_c=2.56 um)
#   bash azure_run_sin_slot_switch.sh sweep        # Fine parameter sweep [2.40 to 2.70 um]
#   bash azure_run_sin_slot_switch.sh amorphous    # Single state run
#   bash azure_run_sin_slot_switch.sh crystalline  # Single state run
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
LC="${2:-2.56}"
WSLOT="${3:-0.065}"
WRAIL="${4:-0.550}"
RES="${5:-25}"

TIMESTAMP=$(date +%Y%m%d_%H%M%S)
LOG_DIR="./logs"
mkdir -p "${LOG_DIR}"
LOG_FILE="${LOG_DIR}/sin_slot_switch_${MODE}_${TIMESTAMP}.log"

echo "================================================================"
echo "  PROJECT JANUS: Si3N4 EMBEDDED PCM SLOT SWITCH (GEOMETRY B)"
echo "  Platform    : Silicon Nitride (Si3N4, n=2.016) on SiO2 (n=1.449)"
echo "  Active PCM  : Sb2S3 directly embedded in-plane between dual rails"
echo "  Dual S-Bends: S-bend tapers capture bent beam in both states"
echo "================================================================"
echo "  Mode        : ${MODE}"
echo "  L_c         : ${LC} um (2:1 resonance)"
echo "  W_slot      : ${WSLOT} um (65 nm central embedded slot)"
echo "  W_rail      : ${WRAIL} um (550 nm Si3N4 rails)"
echo "  Resolution  : ${RES} px/um"
echo "  Log File    : ${LOG_FILE}"
echo "================================================================"
echo ""

START_TIME=$(date +%s)

if [ "${MODE}" = "sweep" ]; then
    echo "Launching fine parameter sweep (L_c=2.40 to 2.70 um, step 0.05 um)..."
    ${RUN_CMD} sin_embedded_slot_switch_3d.py \
        --sweep \
        --lc_min 2.40 --lc_max 2.70 --lc_step 0.05 \
        --w_slot "${WSLOT}" \
        --w_rail "${WRAIL}" \
        --resolution "${RES}" 2>&1 | tee "${LOG_FILE}"
elif [ "${MODE}" = "amorphous" ] || [ "${MODE}" = "crystalline" ]; then
    ${RUN_CMD} sin_embedded_slot_switch_3d.py \
        --state "${MODE}" \
        --lc "${LC}" \
        --w_slot "${WSLOT}" \
        --w_rail "${WRAIL}" \
        --resolution "${RES}" 2>&1 | tee "${LOG_FILE}"
else
    ${RUN_CMD} sin_embedded_slot_switch_3d.py \
        --state both \
        --lc "${LC}" \
        --w_slot "${WSLOT}" \
        --w_rail "${WRAIL}" \
        --resolution "${RES}" 2>&1 | tee "${LOG_FILE}"
fi

END_TIME=$(date +%s)
DURATION=$((END_TIME - START_TIME))

echo ""
echo "================================================================"
echo "  Simulation completed in ${DURATION} seconds."
echo "  Full log saved to: ${LOG_FILE}"
echo "================================================================"
