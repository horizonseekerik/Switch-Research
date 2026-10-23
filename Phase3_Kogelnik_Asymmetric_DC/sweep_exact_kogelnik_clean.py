import sys, os, math, time
import numpy as np
from test_exact_kogelnik_zero_crossing_3d import simulate_kogelnik_exact

if __name__ == "__main__":
    t0 = time.time()
    lc = 4.425
    l_tip = 0.500
    print("=" * 115)
    print(f"CLEAN KOGELNIK GAP SWEEP (105 - 110 nm) AT L_c = {lc:.4f} um, l_tip = {l_tip*1e3:.0f} nm")
    print("=" * 115)
    print(f"{'Gap (nm)':<10} {'AM P2 Modal %':<16} {'AM P1 Modal %':<16} {'CR P1 Modal %':<16} {'CR P2 Modal %':<16} {'AM Tot %':<12} {'CR Tot %':<12}")
    print("-" * 115)

    res_ref = simulate_kogelnik_exact(state="amorphous", resolution=15, l_c=lc, gap=0.108, l_tip=l_tip, verbose=False)
    p_ref_flux = res_ref["p_ref_flux"]
    p_ref_mode = res_ref["p_ref_mode"]

    for g_nm in [106.0, 108.0, 110.0]:
        gap = g_nm * 1e-3
        res_am = simulate_kogelnik_exact(state="amorphous", resolution=15, l_c=lc, gap=gap, l_tip=l_tip,
                                         p_ref_flux=p_ref_flux, p_ref_mode=p_ref_mode, verbose=False)
        res_cr = simulate_kogelnik_exact(state="crystalline", resolution=15, l_c=lc, gap=gap, l_tip=l_tip,
                                         p_ref_flux=p_ref_flux, p_ref_mode=p_ref_mode, verbose=False)

        am_p2_m = res_am["P2_mode"] * 100.0
        am_p1_m = res_am["P1_mode"] * 100.0
        cr_p1_m = res_cr["P1_mode"] * 100.0
        cr_p2_m = res_cr["P2_mode"] * 100.0

        am_tot_m = res_am["tot_mode"] * 100.0
        cr_tot_m = res_cr["tot_mode"] * 100.0

        print(f"{g_nm:<10.2f} {am_p2_m:<16.10f} {am_p1_m:<16.10f} {cr_p1_m:<16.10f} {cr_p2_m:<16.10f} {am_tot_m:<12.6f} {cr_tot_m:<12.6f}")
        sys.stdout.flush()

    print("=" * 115)
    print(f"Elapsed time = {time.time() - t0:.2f} s")
