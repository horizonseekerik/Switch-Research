"""
GAP CALIBRATION SWEEP FOR ALIGNING STATE 0 PEAK WITH STATE 1 NULL (3D FDTD)
===========================================================================
Sweeps gap in [106, 110, 114] nm at L_c = 4.40 um.
Since L_c(0) = pi / (2 * kappa) scales exponentially with gap, widening the gap
moves State 0's peak from 3.60 um to 4.40 um, locking it onto State 1's null!
"""

import sys, os, math, time
import numpy as np

from test_symmetric_pcm_aligned_3d import simulate_symmetric_dc

if __name__ == "__main__":
    t0 = time.time()
    lc = 4.40
    print("=" * 115)
    print(f"CALIBRATION SWEEP ACROSS GAP AT L_c = {lc:.2f} um")
    print("=" * 115)
    print(f"{'Gap (nm)':<10} {'AM P2 Modal %':<16} {'AM P1 Modal %':<16} {'CR P1 Modal %':<16} {'CR P2 Modal %':<16} {'AM P2 Flux %':<16} {'CR P1 Flux %':<16}")
    print("-" * 115)

    for g_nm in [106.0, 110.0, 114.0]:
        gap = g_nm * 1e-3
        res_am = simulate_symmetric_dc(state="amorphous", resolution=15, l_c=lc, gap=gap, verbose=False)
        res_cr = simulate_symmetric_dc(state="crystalline", resolution=15, l_c=lc, gap=gap,
                                       p_ref_flux=res_am["p_ref_flux"], p_ref_mode=res_am["p_ref_mode"], verbose=False)

        am_p2_m = res_am["P2_mode"] * 100.0
        am_p1_m = res_am["P1_mode"] * 100.0
        cr_p1_m = res_cr["P1_mode"] * 100.0
        cr_p2_m = res_cr["P2_mode"] * 100.0

        am_p2_f = res_am["P2_flux"] * 100.0
        cr_p1_f = res_cr["P1_flux"] * 100.0

        print(f"{g_nm:<10.1f} {am_p2_m:<16.10f} {am_p1_m:<16.10f} {cr_p1_m:<16.10f} {cr_p2_m:<16.10f} {am_p2_f:<16.10f} {cr_p1_f:<16.10f}")
        sys.stdout.flush()

    print("=" * 115)
    print(f"Elapsed time = {time.time() - t0:.2f} s")
