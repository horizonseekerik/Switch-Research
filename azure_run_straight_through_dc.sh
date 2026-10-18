#!/usr/bin/env bash
# ==============================================================================
# PROJECT JANUS: AZURE VM LAUNCHER — STRAIGHT-THROUGH 1x2 DC SWITCH
# Architecture: Port 1 (Through, y=0, ZERO BEND) + Port 2 (Drop, S-bent)
# Platform: Azure Standard_D4s_v5 (4 vCPUs)
# ==============================================================================
# Usage:
#   bash azure_run_straight_through_dc.sh              # both states, defaults
#   bash azure_run_straight_through_dc.sh sweep        # targeted parameter sweep
#   bash azure_run_straight_through_dc.sh amorphous    # single state test
#   bash azure_run_straight_through_dc.sh crystalline  # single state test
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
LC="${2:-4.25}"
GAP="${3:-0.100}"
RES="${4:-20}"

TIMESTAMP=$(date +%Y%m%d_%H%M%S)
LOG_DIR="./logs"
mkdir -p "${LOG_DIR}"
LOG_FILE="${LOG_DIR}/straight_dc_${MODE}_${TIMESTAMP}.log"

GAP_NM=$(${PYTHON_EXEC} -c "print(f'{float(\"${GAP}\") * 1000.0:.0f}')")
echo "================================================================"
echo "  PROJECT JANUS: Straight-Through 1x2 DC Switch (3D MEEP FDTD)"
echo "  Architecture : Port 1 = Straight Through (0 dB bend loss)"
echo "                 Port 2 = Evanescent Drop (Low-Loss S-Bend)"
echo "================================================================"
echo "  Mode        : ${MODE}"
echo "  L_c         : ${LC} um (Crystalline Refinement Peak at gap=${GAP}um)"
echo "  Gap         : ${GAP_NM} nm"
echo "  Resolution  : ${RES} px/um"
echo "  Multi-CPU   : ${RUN_CMD}"
echo "  Log File    : ${LOG_FILE}"
echo "================================================================"
echo ""

START_TIME=$(date +%s)

if [ "${MODE}" = "sweep" ]; then
    echo "Launching fine-grid parameter sweep around peak (L_c=4.10 to 4.40 um, gap=100nm)..."
    ${RUN_CMD} straight_through_dc_switch_3d.py \
        --sweep \
        --lc_min 4.10 --lc_max 4.40 --lc_step 0.10 \
        --gaps 0.100 \
        --resolution "${RES}" 2>&1 | tee "${LOG_FILE}"
elif [ "${MODE}" = "amorphous" ] || [ "${MODE}" = "crystalline" ]; then
    ${RUN_CMD} straight_through_dc_switch_3d.py \
        --state "${MODE}" \
        --lc "${LC}" \
        --gap "${GAP}" \
        --resolution "${RES}" 2>&1 | tee "${LOG_FILE}"
else
    # default: both states
    ${RUN_CMD} straight_through_dc_switch_3d.py \
        --state both \
        --lc "${LC}" \
        --gap "${GAP}" \
        --resolution "${RES}" 2>&1 | tee "${LOG_FILE}"
fi

END_TIME=$(date +%s)
ELAPSED=$((END_TIME - START_TIME))

ELAPSED_MIN=$(${PYTHON_EXEC} -c "print(f'{$ELAPSED / 60.0:.2f}')")
echo ""
echo "================================================================"
echo "  COMPLETED! Elapsed: ${ELAPSED} s (~${ELAPSED_MIN} min)"
echo "  JSON Results : ./top_straight_through_dc_candidates.json"
echo "  CSV Results  : ./straight_through_dc_results.csv"
echo "  Log File     : ${LOG_FILE}"

ARCHIVE="straight_dc_results_${MODE}_${TIMESTAMP}.tar.gz"
tar czf "${ARCHIVE}" \
    top_straight_through_dc_candidates.json \
    straight_through_dc_results.csv \
    "${LOG_FILE}" \
    2>/dev/null || \
tar czf "${ARCHIVE}" \
    straight_through_dc_sweep_results.csv \
    top_straight_through_dc_candidates.json \
    "${LOG_FILE}" \
    2>/dev/null || true

echo "  Archive Package: ./${ARCHIVE}"
echo "================================================================"
