"""
PROJECT JANUS: HIGH-PRECISION GAP SCAN (99.0 to 104.0 nm) AT L_c = 4.40 um
==========================================================================
Scans Gap in [99.5, 101.0, 102.5, 104.0] nm to find the exact zero-crossing
where State 0 (amorphous) peak locks onto State 1 (crystalline) null at L_c = 4.40 um.
"""

import sys, os, math, time
import numpy as np

from test_symmetric_pcm_aligned_3d import simulate_symmetric_dc

if __name__ == "__main__":
    t0 = time.time()
    lc = 4.40
    print("=" * 115)
    print(f"HIGH-PRECISION GAP SCAN AT L_c = {lc:.2f} um")
    print("=" * 115)
    print(f"{'Gap (nm)':<10} {'AM P2 Modal %':<16} {'AM P1 Modal %':<16} {'CR P1 Modal %':<16} {'CR P2 Modal %':<16} {'AM P2 Flux %':<16} {'CR P1 Flux %':<16}")
    print("-" * 115)

    res_ref = simulate_symmetric_dc(state="amorphous", resolution=15, l_c=lc, gap=0.101, verbose=False)
    p_ref_flux = res_ref["p_ref_flux"]
    p_ref_mode = res_ref["p_ref_mode"]

    for g_nm in [99.5, 101.0, 102.5, 104.0]:
        gap = g_nm * 1e-3
        res_am = simulate_symmetric_dc(state="amorphous", resolution=15, l_c=lc, gap=gap,
                                       p_ref_flux=p_ref_flux, p_ref_mode=p_ref_mode, verbose=False)
        res_cr = simulate_symmetric_dc(state="crystalline", resolution=15, l_c=lc, gap=gap,
                                       p_ref_flux=p_ref_flux, p_ref_mode=p_ref_mode, verbose=False)

        am_p2_m = res_am["P2_mode"] * 100.0
        am_p1_m = res_am["P1_mode"] * 100.0
        cr_p1_m = res_cr["P1_mode"] * 100.0
        cr_p2_m = res_cr["P2_mode"] * 100.0

        am_p2_f = res_am["P2_flux"] * 100.0
        cr_p1_f = res_cr["P1_flux"] * 100.0

        print(f"{g_nm:<10.2f} {am_p2_m:<16.10f} {am_p1_m:<16.10f} {cr_p1_m:<16.10f} {cr_p2_m:<16.10f} {am_p2_f:<16.10f} {cr_p1_f:<16.10f}")
        sys.stdout.flush()

    print("=" * 115)
    print(f"Elapsed time = {time.time() - t0:.2f} s")
