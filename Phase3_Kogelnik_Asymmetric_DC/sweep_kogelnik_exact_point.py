import sys, os, math, time
import numpy as np
import meep as mp

from test_protected_sub6um_dc import simulate_protected_dc

def run_fine_sweep():
    gap = 0.090
    t_pcm = 0.190
    resolution = 15

    print("=" * 105)
    print("  PROJECT JANUS: EXACT KOGELNIK OPERATING POINT SWEEP")
    print(f"  Parameters: Gap = {gap*1000:.0f} nm, t_pcm = {t_pcm*1000:.0f} nm, Res = {resolution}")
    print("=" * 105)
    print(f"{'L_c(um)':<8} {'AM P2(Cross)%':<15} {'AM P1(Bar)%':<13} {'CR P1(Bar)%':<13} {'CR P2(Cross)%':<15} {'AM Sig/G%':<12} {'CR Sig/G%':<12} {'Area(um2)':<10}")
    print("-" * 105)

    res_ref = simulate_protected_dc("amorphous", resolution=resolution, l_c=3.70, gap=gap, t_pcm=t_pcm, verbose=False)
    p_ref = res_ref["p_ref"]

    for lc in [3.30, 3.50, 3.70, 3.90, 4.10]:
        res_am = simulate_protected_dc("amorphous", resolution=resolution, l_c=lc, gap=gap, t_pcm=t_pcm, p_ref=p_ref, verbose=False)
        res_cr = simulate_protected_dc("crystalline", resolution=resolution, l_c=lc, gap=gap, t_pcm=t_pcm, p_ref=p_ref, verbose=False)

        am_p2 = res_am["p2"] * 100.0
        am_p1 = res_am["p1"] * 100.0
        cr_p1 = res_cr["p1"] * 100.0
        cr_p2 = res_cr["p2"] * 100.0

        am_tot = am_p2 + am_p1
        cr_tot = cr_p1 + cr_p2

        am_ratio = (am_p2 / am_tot) * 100.0 if am_tot > 0 else 0.0
        cr_ratio = (cr_p1 / cr_tot) * 100.0 if cr_tot > 0 else 0.0
        area = res_am["cell_area"]

        print(f"{lc:<8.2f} {am_p2:<15.2f} {am_p1:<13.2f} {cr_p1:<13.2f} {cr_p2:<15.2f} {am_ratio:<12.2f} {cr_ratio:<12.2f} {area:<10.3f}")
        sys.stdout.flush()

    print("=" * 105)

if __name__ == '__main__':
    run_fine_sweep()
