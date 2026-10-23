import sys, os, math, time
import numpy as np
import meep as mp

from test_steep_isolated_dc import simulate_steep_dc

def run_sweep():
    gap = 0.098
    resolution = 15

    print("=" * 95)
    print("  PROJECT JANUS: MAPPING THE EXACT CRYSTALLINE NULL")
    print("  Parameters: gap = 98 nm, t_pcm = 160 nm, res = 15")
    print("=" * 95)
    print(f"{'L_c(um)':<8} {'CR Bar (P1)%':<15} {'CR Leak (P2)%':<15} {'CR Guided%':<12} {'CR XT (dB)':<12} {'CR IL (dB)':<12} {'Area(um2)':<10}")
    print("-" * 95)

    res_ref = simulate_steep_dc("amorphous", resolution=resolution, l_c=4.90, gap=gap, verbose=False)
    p_ref = res_ref["p_ref"]

    for lc in [4.50, 4.70, 4.90, 5.10, 5.30]:
        res_cr = simulate_steep_dc("crystalline", resolution=resolution, l_c=lc, gap=gap, p_ref=p_ref, verbose=False)

        cr_p1 = res_cr["p1"] * 100.0
        cr_p2 = res_cr["p2"] * 100.0
        cr_tot = cr_p1 + cr_p2
        cr_xt = res_cr["xt_db"]
        cr_il = res_cr["il_db"]
        area = res_cr["cell_area"]

        print(f"{lc:<8.2f} {cr_p1:<15.2f} {cr_p2:<15.2f} {cr_tot:<12.2f} {cr_xt:<12.2f} {cr_il:<12.2f} {area:<10.3f}")
        sys.stdout.flush()

    print("=" * 95)

if __name__ == '__main__':
    run_sweep()
