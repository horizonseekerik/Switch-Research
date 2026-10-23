import sys, os, math
import numpy as np
import meep as mp
from meep import mpb
from mpb_candidate3_gap_solver import solve_candidate3_supermodes

print("=" * 115)
print("  FINE-GRAINED MPB ZERO-CROSSING SCAN (98.0 - 101.0 nm) AT RESOLUTION 32")
print("=" * 115)
print(f"{'Gap (nm)':<10} {'dn_am':<12} {'L_c(0) (um)':<14} {'dn_cr':<12} {'L_c(1) (um)':<14} {'S/kappa':<12} {'Offset(um)':<12}")
print("-" * 115)

for g_nm in np.arange(98.0, 101.2, 0.4):
    g = g_nm * 1e-3
    res = solve_candidate3_supermodes(g, resolution=32)
    print(f"{g_nm:<10.2f} {res['dn_am']:<12.6f} {res['lc_0']:<14.6f} {res['dn_cr']:<12.6f} {res['lc_1']:<14.6f} {res['ratio_S_kappa']:<12.6f} {res['offset']:<12.6f}")
    sys.stdout.flush()

print("=" * 115)
