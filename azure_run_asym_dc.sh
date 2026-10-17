#!/bin/bash
# ================================================================
# PROJECT JANUS: Azure VM Launcher — Asymmetric DC Switch
# ================================================================
# Usage:
#   bash azure_run_asym_dc.sh                     # both states, default params
#   bash azure_run_asym_dc.sh sweep               # full parameter sweep
#   bash azure_run_asym_dc.sh amorphous 4.20 0.085 0.270
#   bash azure_run_asym_dc.sh crystalline 4.20 0.085 0.270
#
# Design:
#   WG1 = 240 nm, carries Sb2S3 PCM (switchable)
#   WG2 = wider bare Si (w2 tuned so n_eff_WG2 ≈ n_eff_WG1_am → Δβ=0)
#   → Amorphous : resonant cross coupling → Port 2 ≥ 98%
#   → Crystalline: large Δβ quenches coupling → Port 1 ≥ 98%
#
# Platform: Azure Standard_D4s_v5 (4 vCPUs)
# ================================================================

set -e

export OMP_NUM_THREADS=4
export OPENBLAS_NUM_THREADS=4
export MKL_NUM_THREADS=4
export NUMEXPR_NUM_THREADS=4

MODE="${1:-both}"
LC="${2:-4.20}"
GAP="${3:-0.085}"
W2="${4:-0.270}"
RES="${5:-20}"

TIMESTAMP=$(date +%Y%m%d_%H%M%S)
LOG_DIR="./logs"
mkdir -p "${LOG_DIR}"
LOG_FILE="${LOG_DIR}/asym_dc_${MODE}_${TIMESTAMP}.log"

echo "================================================================"
echo "  PROJECT JANUS: Asymmetric 1x2 DC Switch (3D MEEP FDTD)"
echo "================================================================"
echo "  Mode       : ${MODE}"
echo "  L_c        : ${LC} um"
echo "  Gap        : $(echo "${GAP} * 1000" | bc -l | xargs printf '%.0f') nm"
echo "  W2 (bare)  : $(echo "${W2} * 1000" | bc -l | xargs printf '%.0f') nm"
echo "  Resolution : ${RES} px/um"
echo "  OMP Threads: ${OMP_NUM_THREADS}"
echo "  Log        : ${LOG_FILE}"
echo "================================================================"
echo ""

START_TIME=$(date +%s)

if [ "${MODE}" = "sweep" ]; then
    echo "Running full parameter sweep (w2 × gap × L_c × 2 states)..."
    python3 asymmetric_dc_switch_3d.py \
        --sweep \
        --resolution "${RES}" 2>&1 | tee "${LOG_FILE}"
elif [ "${MODE}" = "amorphous" ] || [ "${MODE}" = "crystalline" ]; then
    python3 asymmetric_dc_switch_3d.py \
        --state "${MODE}" \
        --lc "${LC}" \
        --gap "${GAP}" \
        --w2 "${W2}" \
        --resolution "${RES}" 2>&1 | tee "${LOG_FILE}"
else
    # both states
    python3 asymmetric_dc_switch_3d.py \
        --state both \
        --lc "${LC}" \
        --gap "${GAP}" \
        --w2 "${W2}" \
        --resolution "${RES}" 2>&1 | tee "${LOG_FILE}"
fi

END_TIME=$(date +%s)
ELAPSED=$((END_TIME - START_TIME))

echo ""
echo "================================================================"
echo "  COMPLETED! Elapsed: ${ELAPSED} s ($(echo "scale=2; ${ELAPSED}/60" | bc) min)"
echo "  Results  : ./top_asym_dc_candidates.json"
echo "  CSV      : ./asym_dc_results.csv (or asym_dc_sweep_results.csv)"
echo "  Log      : ${LOG_FILE}"

ARCHIVE="asym_dc_results_${MODE}_${TIMESTAMP}.tar.gz"
tar czf "${ARCHIVE}" \
    top_asym_dc_candidates.json \
    asym_dc_results.csv \
    "${LOG_FILE}" \
    2>/dev/null || \
tar czf "${ARCHIVE}" \
    asym_dc_sweep_results.csv \
    top_asym_dc_candidates.json \
    "${LOG_FILE}" \
    2>/dev/null || true
echo "  Archive  : ./${ARCHIVE}"
echo "================================================================"
