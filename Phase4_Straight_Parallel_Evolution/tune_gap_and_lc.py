"""
PROJECT JANUS: GAP & LENGTH 2D OPTIMIZER FOR PROTECTED SWITCH
=============================================================
Calibrates gap in [105, 110, 115] nm and L_c in [3.4, 4.2] um
to find the exact point where u = sqrt(3) is satisfied in 3D FDTD,
achieving simultaneous > 95% transmission in BOTH states with < 5.2 um^2 area.
"""

import sys, os, math
import meep as mp
from test_protected_sub6um_dc import simulate_protected_dc

def run_grid():
    grid = [
        # gap, lc
        (0.105, 3.40),
        (0.105, 3.70),
        (0.105, 4.00),
        (0.110, 3.60),
        (0.110, 3.90),
        (0.110, 4.20),
    ]

    print("=" * 95)
    print("  PROJECT JANUS: GAP & LENGTH 2D OPTIMIZER")
    print("=" * 95)
    print(f"{'Gap(nm)':<10} {'L_c(um)':<10} {'AM P2(Cross)%':<15} {'CR P1(Bar)%':<15} {'AM XT(dB)':<12} {'CR XT(dB)':<12} {'Area(um2)':<10}")
    print("-" * 95)

    res_ref = simulate_protected_dc("amorphous", resolution=15, l_c=3.80, gap=0.105, verbose=False)
    p_ref = res_ref["p_ref"]

    for gap, lc in grid:
        res_am = simulate_protected_dc("amorphous", resolution=15, l_c=lc, gap=gap, p_ref=p_ref, verbose=False)
        res_cr = simulate_protected_dc("crystalline", resolution=15, l_c=lc, gap=gap, p_ref=p_ref, verbose=False)

        am_p2 = res_am["p2"] * 100.0
        cr_p1 = res_cr["p1"] * 100.0
        am_xt = res_am["xt_db"]
        cr_xt = res_cr["xt_db"]
        area = res_am["cell_area"]

        print(f"{gap*1000:<10.0f} {lc:<10.2f} {am_p2:<15.2f} {cr_p1:<15.2f} {am_xt:<12.2f} {cr_xt:<12.2f} {area:<10.3f}")
        sys.stdout.flush()

if __name__ == '__main__':
    run_grid()
