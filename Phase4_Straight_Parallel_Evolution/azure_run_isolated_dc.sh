#!/usr/bin/env bash
# ==============================================================================
# PROJECT JANUS: AZURE VM LAUNCHER — ISOLATED-TAPERS 1x2 DC SWITCH
# ==============================================================================
# Architecture:
#   - Isolated adiabatic PCM tapers (0.50 um) placed outside the coupling zone
#   - Pure uniform coupling zone (L_c = 4.40 um, gap = 100 nm, w = 240 nm)
#   - Zero inter-waveguide taper scattering
#   - Zero S-bend fanout leakage
# ==============================================================================

set -e

# ── 1. Conda Activation ───────────────────────────────────────────────────────
for CONDA_PATH in "$HOME/miniforge3" "$HOME/miniconda3" "$HOME/anaconda3" "/opt/conda"; do
    if [ -f "${CONDA_PATH}/etc/profile.d/conda.sh" ]; then
        source "${CONDA_PATH}/etc/profile.d/conda.sh"
        conda activate meep 2>/dev/null || conda activate meep_env 2>/dev/null || conda activate base 2>/dev/null || true
        break
    fi
done

# ── 2. Python Selection ───────────────────────────────────────────────────────
PYTHON_EXEC="python3"
if ! $PYTHON_EXEC -c "import meep" >/dev/null 2>&1; then
    if python -c "import meep" >/dev/null 2>&1; then
        PYTHON_EXEC="python"
    fi
fi

# ── 3. Multi-CPU Configuration (4 vCPUs) ──────────────────────────────────────
export OMP_NUM_THREADS=4
export OPENBLAS_NUM_THREADS=4
export MKL_NUM_THREADS=4
export NUMEXPR_NUM_THREADS=4

NUM_CPUS=4
RUN_CMD="${PYTHON_EXEC}"

if command -v mpirun >/dev/null 2>&1; then
    if mpirun -np 2 "${PYTHON_EXEC}" -c "import meep" >/dev/null 2>&1; then
        echo "[INFO] OpenMPI active! Utilizing ${NUM_CPUS} CPU cores with mpirun."
        RUN_CMD="mpirun -np ${NUM_CPUS} ${PYTHON_EXEC}"
    else
        echo "[INFO] Running with OpenMP multi-threading (OMP_NUM_THREADS=${NUM_CPUS})."
    fi
else
    echo "[INFO] Running with OpenMP multi-threading (OMP_NUM_THREADS=${NUM_CPUS})."
fi

# ── 4. Parameters ────────────────────────────────────────────────────────────
MODE="${1:-both}"
LC="${2:-4.40}"
GAP="${3:-0.100}"
RES="${4:-20}"
TAPER="${5:-0.500}"

TIMESTAMP=$(date +%Y%m%d_%H%M%S)
LOG_DIR="./logs"
mkdir -p "${LOG_DIR}"
LOG_FILE="${LOG_DIR}/isolated_dc_${MODE}_${TIMESTAMP}.log"

GAP_NM=$(${PYTHON_EXEC} -c "print(f'{float(\"${GAP}\") * 1000.0:.0f}')")
TAPER_NM=$(${PYTHON_EXEC} -c "print(f'{float(\"${TAPER}\") * 1000.0:.0f}')")

echo "================================================================"
echo "  PROJECT JANUS: Isolated-Tapers 1x2 DC Switch (3D MEEP FDTD)"
echo "  Architecture : Pure Uniform Coupling Zone + Isolated Tapers"
echo "  Features     : Zero Inter-Taper Scattering | Zero Fanout Leakage"
echo "================================================================"
echo "  Mode        : ${MODE}"
echo "  L_c         : ${LC} um (Pure Uniform Coupling Zone)"
echo "  Gap         : ${GAP_NM} nm"
echo "  Taper       : ${TAPER_NM} nm (Isolated from Other Arm)"
echo "  Resolution  : ${RES} px/um"
echo "  Multi-CPU   : ${RUN_CMD}"
echo "  Log File    : ${LOG_FILE}"
echo "================================================================"
echo ""

START_TIME=$(date +%s)

if [ "${MODE}" = "sweep" ]; then
    echo "Launching fine-grid parameter sweep around peak (L_c=4.20 to 4.60 um)..."
    ${RUN_CMD} straight_through_isolated_dc_3d.py \
        --sweep \
        --lc_min 4.20 --lc_max 4.60 --lc_step 0.10 \
        --gaps "${GAP}" \
        --taper "${TAPER}" \
        --resolution "${RES}" 2>&1 | tee "${LOG_FILE}"
elif [ "${MODE}" = "amorphous" ] || [ "${MODE}" = "crystalline" ]; then
    ${RUN_CMD} straight_through_isolated_dc_3d.py \
        --state "${MODE}" \
        --lc "${LC}" \
        --gap "${GAP}" \
        --taper "${TAPER}" \
        --resolution "${RES}" 2>&1 | tee "${LOG_FILE}"
else
    # default: both states
    ${RUN_CMD} straight_through_isolated_dc_3d.py \
        --state both \
        --lc "${LC}" \
        --gap "${GAP}" \
        --taper "${TAPER}" \
        --resolution "${RES}" 2>&1 | tee "${LOG_FILE}"
fi

END_TIME=$(date +%s)
ELAPSED=$((END_TIME - START_TIME))
ELAPSED_MIN=$(${PYTHON_EXEC} -c "print(f'{$ELAPSED / 60.0:.2f}')")

echo ""
echo "================================================================"
echo "  COMPLETED! Elapsed: ${ELAPSED} s (~${ELAPSED_MIN} min)"
echo "  JSON Results : ./top_isolated_dc_candidates.json"
echo "  CSV Results  : ./isolated_dc_results.csv"
echo "  Log File     : ${LOG_FILE}"
echo "================================================================"
