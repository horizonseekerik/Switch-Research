import sys, os
from optimize_cmt_exact import simulate_coupler_cmode

for gap_nm in [87, 88, 89, 90, 91]:
    gap = gap_nm / 1000.0
    print(f"\n=== Gap = {gap_nm} nm ===")
    for lc in [3.56, 3.58, 3.60, 3.62, 3.64, 3.66, 3.68]:
        p2_am, p1_cr, phi_ends, k, s = simulate_coupler_cmode(lc, gap)
        both_pass = (p2_am >= 0.98) and (p1_cr >= 0.98)
        print(f"L_c = {lc:.2f} um: AM P2 = {p2_am*100:.2f}%, CR P1 = {p1_cr*100:.2f}% | Both >98%: {both_pass}")
