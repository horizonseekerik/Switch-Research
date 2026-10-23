import sys, os, math
import numpy as np
import meep as mp

from test_kogelnik_dc_switch import simulate_kogelnik_switch

def sweep_lc():
    # Cache reference flux once
    res_0 = simulate_kogelnik_switch("amorphous", resolution=15, l_c=4.50, gap=0.098, t_pcm=0.160)
    p_ref = res_0["p_ref"]

    print("\n" + "=" * 80)
    print("  SWEEPING L_c FOR MAXIMUM CROSS-COUPLING IN AMORPHOUS STATE (res=15)")
    print("=" * 80)
    print(f"{'L_c (um)':<10} {'Port 1 (Bar)%':<15} {'Port 2 (Cross)%':<16} {'Rad Loss%':<12} {'IL (dB)':<10} {'XT (dB)':<10}")
    print("-" * 80)

    best_p2 = 0.0
    best_lc = 4.50

    for lc in [4.00, 4.20, 4.40, 4.60, 4.80, 5.00, 5.20]:
        res = simulate_kogelnik_switch("amorphous", resolution=15, l_c=lc, gap=0.098, t_pcm=0.160, p_ref=p_ref)
        p1 = res["p1"] * 100
        p2 = res["p2"] * 100
        rad = (1.0 - res["p_tot"]) * 100
        il = res["il_db"]
        xt = res["xt_db"]
        print(f"{lc:<10.2f} {p1:<15.2f} {p2:<16.2f} {rad:<12.2f} {il:<10.3f} {xt:<10.2f}")
        if p2 > best_p2:
            best_p2 = p2
            best_lc = lc

    print("=" * 80)
    print(f"Optimal L_c for Amorphous state: L_c = {best_lc:.2f} um (Port 2 = {best_p2:.2f}%)")
    print("=" * 80)

if __name__ == '__main__':
    sweep_lc()
