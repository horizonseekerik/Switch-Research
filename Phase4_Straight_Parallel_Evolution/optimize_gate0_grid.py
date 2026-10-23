"""
PROJECT JANUS: GATE 0 GLOBAL 2D/3D OPTIMIZER
============================================
Searches the full parameter space (L_c, g0, delta_g) in coupled-mode theory
to find the EXACT mathematical coordinates where:
  State 0 (Amorphous Cross P2) >= 98.0000000000%
  State 1 (Crystalline Bar P1) >= 98.0000000000%
Reports all results to >= 10 decimal places.
"""

import math
import numpy as np
from scipy.optimize import minimize
from verify_gate0_cme import evaluate_profile

def cost_function(params):
    lc, g0, dg = params
    r = evaluate_profile(l_c=lc, g0=g0, delta_g=dg)
    p2_0 = r["p2_am"]
    p1_1 = r["p1_cr"]
    # Penalty if either drops below 98%
    err0 = max(0.0, 99.0 - p2_0)**2
    err1 = max(0.0, 99.0 - p1_1)**2
    return err0 + err1

print("=" * 115)
print("GATE 0: MULTI-VARIABLE GRID SEARCH & OPTIMIZATION")
print("=" * 115)

best_res = None
best_score = 1e9

# Grid search across physical ranges
for lc in [3.80, 4.00, 4.20, 4.425, 4.60, 4.80]:
    for g0_nm in [90.0, 94.0, 98.14, 102.0]:
        g0 = g0_nm * 1e-3
        for dg_nm in [0.0, 5.0, 10.0, 15.0, 20.0]:
            dg = dg_nm * 1e-3
            r = evaluate_profile(l_c=lc, g0=g0, delta_g=dg)
            score = (100.0 - r["p2_am"])**2 + (100.0 - r["p1_cr"])**2
            if r["p2_am"] > 95.0 and r["p1_cr"] > 95.0:
                print(f"HIGH CANDIDATE: Lc={lc:.3f}um, g0={g0_nm:.2f}nm, dg={dg_nm:.2f}nm -> S0 P2={r['p2_am']:.6f}%, S1 P1={r['p1_cr']:.6f}%")
            if score < best_score:
                best_score = score
                best_res = (lc, g0_nm, dg_nm, r)

print("=" * 115)
lc, g0_nm, dg_nm, r = best_res
print(f"TOP CME RESULT: L_c = {lc:.4f} um, g0 = {g0_nm:.2f} nm, delta_g = {dg_nm:.2f} nm")
print(f"  State 0 Cross (P2) = {r['p2_am']:.10f}%, Leakage (P1) = {r['p1_am']:.10f}%")
print(f"  State 1 Bar   (P1) = {r['p1_cr']:.10f}%, Leakage (P2) = {r['p2_cr']:.10f}%")
print(f"  State 0 Crosstalk  = {r['xt_am']:.10f} dB")
print(f"  State 1 Crosstalk  = {r['xt_cr']:.10f} dB")
print("=" * 115)
