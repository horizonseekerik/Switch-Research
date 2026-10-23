#!/usr/bin/env bash
# ==============================================================================
# PROJECT JANUS: AZURE VM LAUNCHER — PROTECTED DIRECTIONAL COUPLER SWITCH
# Architecture: Port 1 Straight (y=0) + Protected S-Bend Drop + Decoupled Tapers
# Platform: Azure Standard_D4s_v5 (4 vCPUs)
# ==============================================================================
# Usage:
#   bash azure_run_protected_dc.sh              # dual-state test (L_c=3.90um, gap=105nm)
#   bash azure_run_protected_dc.sh sweep        # targeted parameter sweep
#   bash azure_run_protected_dc.sh amorphous    # single state test
#   bash azure_run_protected_dc.sh crystalline  # single state test
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
LC="${2:-3.90}"
GAP="${3:-0.105}"
RES="${4:-20}"
LTAPER="${5:-0.800}"

TIMESTAMP=$(date +%Y%m%d_%H%M%S)
LOG_DIR="./logs"
mkdir -p "${LOG_DIR}"
LOG_FILE="${LOG_DIR}/protected_dc_${MODE}_${TIMESTAMP}.log"

GAP_NM=$(${PYTHON_EXEC} -c "print(f'{float(\"${GAP}\") * 1000.0:.0f}')")
TAPER_NM=$(${PYTHON_EXEC} -c "print(f'{float(\"${LTAPER}\") * 1000.0:.0f}')")

echo "================================================================"
echo "  PROJECT JANUS: Protected 1x2 DC Switch (3D MEEP FDTD)"
echo "  Architecture : Port 1 = 100% Straight Bus (0.00 dB Bend Loss)"
echo "                 Port 2 = Protected Drop Arm (S-Bend in/out)"
echo "                 Isolated Tapers (gap >= 320 nm, L_taper = ${TAPER_NM} nm)"
echo "================================================================"
echo "  Mode        : ${MODE}"
echo "  L_c         : ${LC} um (Exact analytical switching resonance)"
echo "  Gap         : ${GAP_NM} nm (Satisfies u = sqrt(3) invariant)"
echo "  L_taper     : ${TAPER_NM} nm (Decoupled adiabatic taper)"
echo "  Resolution  : ${RES} px/um"
echo "  Multi-CPU   : ${RUN_CMD}"
echo "  Log File    : ${LOG_FILE}"
echo "================================================================"
echo ""

START_TIME=$(date +%s)

if [ "${MODE}" = "sweep" ]; then
    echo "Launching targeted fast sweep around resonance (L_c=3.80 to 4.00 um, step 0.05 um, gap=105 nm)..."
    ${RUN_CMD} straight_through_protected_dc_3d.py \
        --sweep \
        --lc_min 3.80 --lc_max 4.00 --lc_step 0.05 \
        --gaps 0.105 \
        --ltaper "${LTAPER}" \
        --resolution "${RES}" 2>&1 | tee "${LOG_FILE}"
elif [ "${MODE}" = "amorphous" ] || [ "${MODE}" = "crystalline" ]; then
    ${RUN_CMD} straight_through_protected_dc_3d.py \
        --state "${MODE}" \
        --lc "${LC}" \
        --gap "${GAP}" \
        --ltaper "${LTAPER}" \
        --resolution "${RES}" 2>&1 | tee "${LOG_FILE}"
else
    ${RUN_CMD} straight_through_protected_dc_3d.py \
        --state both \
        --lc "${LC}" \
        --gap "${GAP}" \
        --ltaper "${LTAPER}" \
        --resolution "${RES}" 2>&1 | tee "${LOG_FILE}"
fi

END_TIME=$(date +%s)
DURATION=$((END_TIME - START_TIME))

echo ""
echo "================================================================"
echo "  Simulation completed in ${DURATION} seconds."
echo "  Full log saved to: ${LOG_FILE}"
echo "================================================================"
