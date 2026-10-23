import sys, os, math, time
import numpy as np
import meep as mp

from test_steep_isolated_dc import simulate_steep_dc

def run_sweep():
    gap = 0.112
    resolution = 15

    print("=" * 95)
    print("  PROJECT JANUS: HALF-BEAT KOGELNIK SWEEP (Gap = 112 nm)")
    print("=" * 95)
    print(f"{'L_c(um)':<8} {'AM P2(Cross)%':<15} {'AM P1(Bar)%':<13} {'CR P1(Bar)%':<13} {'CR P2(Cross)%':<15} {'AM Sig/G':<11} {'CR Sig/G':<11} {'Area(um2)':<10}")
    print("-" * 95)

    res_ref = simulate_steep_dc("amorphous", resolution=resolution, l_c=2.25, gap=gap, verbose=False)
    p_ref = res_ref["p_ref"]

    for lc in [2.00, 2.15, 2.25, 2.35, 2.50]:
        res_am = simulate_steep_dc("amorphous", resolution=resolution, l_c=lc, gap=gap, p_ref=p_ref, verbose=False)
        res_cr = simulate_steep_dc("crystalline", resolution=resolution, l_c=lc, gap=gap, p_ref=p_ref, verbose=False)

        am_p2 = res_am["p2"] * 100.0
        am_p1 = res_am["p1"] * 100.0
        cr_p1 = res_cr["p1"] * 100.0
        cr_p2 = res_cr["p2"] * 100.0

        am_tot = am_p2 + am_p1
        cr_tot = cr_p1 + cr_p2

        am_ratio = (am_p2 / am_tot) * 100.0 if am_tot > 0 else 0.0
        cr_ratio = (cr_p1 / cr_tot) * 100.0 if cr_tot > 0 else 0.0
        area = res_am["cell_area"]

        print(f"{lc:<8.2f} {am_p2:<15.2f} {am_p1:<13.2f} {cr_p1:<13.2f} {cr_p2:<15.2f} {am_ratio:<11.2f}% {cr_ratio:<11.2f}% {area:<10.3f}")
        sys.stdout.flush()

    print("=" * 95)

if __name__ == '__main__':
    run_sweep()
