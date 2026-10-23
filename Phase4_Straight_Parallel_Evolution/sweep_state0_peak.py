"""
FAST SWEEP FOR STATE 0 RESONANCE PEAK AT GAP = 98.14 nm
======================================================
Finds the exact physical L_c where State 0 reaches > 98% modal power.
"""

import sys, os, math, time
import numpy as np

from test_exact_kogelnik_zero_crossing_3d import simulate_kogelnik_exact

if __name__ == "__main__":
    t0 = time.time()
    gap = 0.09814
    print("=" * 105)
    print(f"SWEEPING STATE 0 (AMORPHOUS) ACROSS L_c AT GAP = {gap*1e3:.2f} nm")
    print("=" * 105)
    print(f"{'L_c (um)':<10} {'P2 Modal %':<16} {'P1 Modal %':<16} {'Total Modal %':<16} {'P2 Flux %':<16} {'IL (dB)':<10}")
    print("-" * 105)

    res_ref = simulate_kogelnik_exact(state="amorphous", resolution=15, l_c=3.80, gap=gap, verbose=False)
    p_ref_flux = res_ref["p_ref_flux"]
    p_ref_mode = res_ref["p_ref_mode"]

    for lc in [3.60, 3.80, 4.00, 4.20, 4.425]:
        res = simulate_kogelnik_exact(state="amorphous", resolution=15, l_c=lc, gap=gap,
                                      p_ref_flux=p_ref_flux, p_ref_mode=p_ref_mode, verbose=False)
        p2_m = res["P2_mode"] * 100.0
        p1_m = res["P1_mode"] * 100.0
        tot_m = res["tot_mode"] * 100.0
        p2_f = res["P2_flux"] * 100.0
        il = res["il_mode"]
        print(f"{lc:<10.3f} {p2_m:<16.10f} {p1_m:<16.10f} {tot_m:<16.10f} {p2_f:<16.10f} {il:<10.4f}")
        sys.stdout.flush()

    print("=" * 105)
    print(f"Elapsed time = {time.time() - t0:.2f} s")
