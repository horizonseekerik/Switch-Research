import sys, os, math
import numpy as np
import meep as mp

from test_kogelnik_dc_switch import simulate_kogelnik_switch

def sweep_peak():
    # Cache reference flux once at lc=3.70
    res_0 = simulate_kogelnik_switch("amorphous", resolution=15, l_c=3.70, gap=0.098, t_pcm=0.160)
    p_ref = res_0["p_ref"]

    print("\n" + "=" * 80)
    print("  PEAK HUNT: SWEEPING L_c IN [3.40, 4.00] um (Amorphous State, res=15)")
    print("=" * 80)
    print(f"{'L_c (um)':<10} {'Port 1 (Bar)%':<15} {'Port 2 (Cross)%':<16} {'Rad Loss%':<12} {'IL (dB)':<10} {'XT (dB)':<10}")
    print("-" * 80)

    best_p2 = 0.0
    best_lc = 3.70

    for lc in [3.40, 3.50, 3.60, 3.70, 3.80, 3.90, 4.00]:
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

    # Now verify crystalline state at the best_lc
    print(f"\nVerifying Crystalline state at optimal L_c = {best_lc:.2f} um...")
    res_cr = simulate_kogelnik_switch("crystalline", resolution=15, l_c=best_lc, gap=0.098, t_pcm=0.160, p_ref=p_ref)
    p1_cr = res_cr["p1"] * 100
    p2_cr = res_cr["p2"] * 100
    il_cr = res_cr["il_db"]
    xt_cr = res_cr["xt_db"]

    print("\n" + "=" * 80)
    print(f"DUAL-STATE PERFORMANCE AT L_c = {best_lc:.2f} um:")
    print(f"  State 0 (Amorphous   -> Port 2 Cross) : Signal = {best_p2:.2f}%, Leak = {100-best_p2-rad:.2f}%")
    print(f"  State 1 (Crystalline -> Port 1 Bar)   : Signal = {p1_cr:.2f}%, Leak = {p2_cr:.2f}%, XT = {xt_cr:.2f} dB")
    print("=" * 80)

if __name__ == '__main__':
    sweep_peak()
