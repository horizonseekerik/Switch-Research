#!/usr/bin/env bash
# ==============================================================================
# PROJECT JANUS: AZURE VM MULTI-CPU LAUNCHER — REFINED SWITCH ARCHITECTURES
# ==============================================================================
# Supported Architectures:
#   1. option1  : Thin PCM Symmetric (w1=w2=240nm, t_pcm=30nm on both arms)
#   2. option2  : Asymmetric Cores (w1=240nm with 150nm PCM, w2=340nm bare Si)
#   3. combined : Thin PCM + Asymmetric Phase-Matched Core (Top Candidate)
#   4. all      : Runs all 3 architectures comparatively
#   5. mpb      : Runs 3D vector MPB phase-matching mode solver (fast verification)
#
# Usage:
#   bash azure_run_refinements.sh combined             # Run combined top candidate (both states)
#   bash azure_run_refinements.sh all                  # Run all 3 architectures & print comparison table
#   bash azure_run_refinements.sh option1              # Run Option 1
#   bash azure_run_refinements.sh option2              # Run Option 2
#   bash azure_run_refinements.sh mpb                  # Run MPB vector mode solver
#   bash azure_run_refinements.sh sweep combined       # Run parameter sweep for combined architecture
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
        echo "[INFO] OpenMPI active! Utilizing ${NUM_CPUS} CPU cores with mpirun."
        RUN_CMD="mpirun -np ${NUM_CPUS} ${PYTHON_EXEC}"
    else
        echo "[INFO] Running with OpenMP multi-threading (OMP_NUM_THREADS=${NUM_CPUS})."
    fi
else
    echo "[INFO] Running with OpenMP multi-threading (OMP_NUM_THREADS=${NUM_CPUS})."
fi

# ── 4. Parse Arguments ────────────────────────────────────────────────────────
CMD_ARG="${1:-combined}"
TIMESTAMP=$(date +%Y%m%d_%H%M%S)
LOG_DIR="./logs"
mkdir -p "${LOG_DIR}"

if [ "${CMD_ARG}" = "mpb" ]; then
    echo "================================================================"
    echo "  PROJECT JANUS: MPB 3D Vector Phase-Matching & Supermode Solver"
    echo "================================================================"
    ${PYTHON_EXEC} solve_exact_phase_matching.py 2>&1 | tee "${LOG_DIR}/mpb_solve_${TIMESTAMP}.log"
    exit 0
fi

if [ "${CMD_ARG}" = "sweep" ]; then
    ARCH="${2:-combined}"
    RES="${3:-20}"
    LOG_FILE="${LOG_DIR}/refined_sweep_${ARCH}_${TIMESTAMP}.log"
    echo "================================================================"
    echo "  PROJECT JANUS: Parameter Sweep for [${ARCH^^}]"
    echo "  Resolution  : ${RES} px/um"
    echo "  Log File    : ${LOG_FILE}"
    echo "================================================================"
    START_TIME=$(date +%s)
    ${RUN_CMD} straight_through_refined_switch_3d.py \
        --arch "${ARCH}" \
        --sweep \
        --resolution "${RES}" 2>&1 | tee "${LOG_FILE}"
    END_TIME=$(date +%s)
    ELAPSED=$((END_TIME - START_TIME))
    echo "Completed in ${ELAPSED} s."
    exit 0
fi

ARCH="${CMD_ARG}"
RES="${2:-20}"
LOG_FILE="${LOG_DIR}/refined_${ARCH}_${TIMESTAMP}.log"

echo "================================================================"
echo "  PROJECT JANUS: 3D FDTD Refined Directional Coupler Switch"
echo "  Architecture: ${ARCH^^}"
echo "  Resolution  : ${RES} px/um"
echo "  Multi-CPU   : ${RUN_CMD}"
echo "  Log File    : ${LOG_FILE}"
echo "================================================================"
echo ""

START_TIME=$(date +%s)

if [ "${ARCH}" = "all" ]; then
    echo "Running all 3 refinement architectures comparatively..."
    ${RUN_CMD} straight_through_refined_switch_3d.py \
        --arch all \
        --resolution "${RES}" 2>&1 | tee "${LOG_FILE}"
else
    ${RUN_CMD} straight_through_refined_switch_3d.py \
        --arch "${ARCH}" \
        --state both \
        --resolution "${RES}" 2>&1 | tee "${LOG_FILE}"
fi

END_TIME=$(date +%s)
ELAPSED=$((END_TIME - START_TIME))
ELAPSED_MIN=$(${PYTHON_EXEC} -c "print(f'{$ELAPSED / 60.0:.1f}')")

echo ""
echo "================================================================"
echo "  SIMULATION COMPLETED!"
echo "  Total Runtime : ${ELAPSED} s (~${ELAPSED_MIN} min)"
echo "  Log File      : ${LOG_FILE}"

ARCHIVE="refined_${ARCH}_results_${TIMESTAMP}.tar.gz"
tar czf "${ARCHIVE}" \
    top_refined_switch_candidates.json \
    "${LOG_FILE}" \
    2>/dev/null || \
tar czf "${ARCHIVE}" \
    "${LOG_FILE}" \
    2>/dev/null || true

echo "  Results Archive: ./${ARCHIVE}"
echo "================================================================"
