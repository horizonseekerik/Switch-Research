"""
PROJECT JANUS: 3D FDTD MEEP SIMULATION - CANDIDATE 3 AT RESONANCE GAP (98.2 nm)
=================================================================================
Exact MPB Supermode Solution:
  Gap = 98.2000000000 nm -> S / kappa = 2.0012481023 (Theoretical Ideal = 2.0000000000)
  L_c(0) = 4.4250000000 um
  L_c(1) = 4.4260000000 um
  Beat-length offset = 0.0010000000 um (1.0 nm!)

Sweep range: L_c = [4.30, 4.38, 4.44, 4.50] um in full 3D Maxwell FDTD (PyMEEP).
"""

import sys, os, math, time, argparse
import numpy as np
import meep as mp
from test_aligned_detuned_switch import simulate_aligned_switch

def run_aligned_resonance_sweep():
    gap = 0.0982000000
    t_pcm = 0.2000000000
    res = 15
    y_out = 0.4200000000
    l_bend = 0.9000000000
    w_bend = 0.3200000000

    print("=" * 115)
    print("  PROJECT JANUS: 3D FDTD MEEP VERIFICATION OF CANDIDATE 3 AT RESONANCE GAP (98.2 nm)")
    print(f"  Configuration: Gap = {gap*1e3:.4f} nm, t_pcm = {t_pcm*1e3:.1f} nm, Resolution = {res}")
    print("=" * 115)
    print(f"{'L_c(um)':<8} {'AM P2(Mod)%':<14} {'AM P1(Mod)%':<13} {'CR P1(Mod)%':<13} {'CR P2(Mod)%':<14} {'AM Sig/G%':<12} {'CR Sig/G%':<12} {'Area(um2)':<10}")
    print("-" * 115)

    # Reference straight waveguide at center of sweep
    res_ref = simulate_aligned_switch("amorphous", resolution=res, l_c=4.40, gap=gap, t_pcm=t_pcm,
                                      y_out=y_out, l_bend=l_bend, w_bend=w_bend, verbose=False)
    p_ref_flux = res_ref["p_ref_flux"]
    p_ref_mode = res_ref["p_ref_mode"]

    for lc in [4.32, 4.38, 4.44, 4.50]:
        res_am = simulate_aligned_switch("amorphous", resolution=res, l_c=lc, gap=gap, t_pcm=t_pcm,
                                         y_out=y_out, l_bend=l_bend, w_bend=w_bend,
                                         p_ref_flux=p_ref_flux, p_ref_mode=p_ref_mode, verbose=False)
        res_cr = simulate_aligned_switch("crystalline", resolution=res, l_c=lc, gap=gap, t_pcm=t_pcm,
                                         y_out=y_out, l_bend=l_bend, w_bend=w_bend,
                                         p_ref_flux=p_ref_flux, p_ref_mode=p_ref_mode, verbose=False)

        am_p2 = res_am["p2_mode"] * 100.0
        am_p1 = res_am["p1_mode"] * 100.0
        cr_p1 = res_cr["p1_mode"] * 100.0
        cr_p2 = res_cr["p2_mode"] * 100.0
        am_ratio = res_am["ratio_guided"]
        cr_ratio = res_cr["ratio_guided"]
        area = res_am["cell_area"]

        print(f"{lc:<8.2f} {am_p2:<14.6f} {am_p1:<13.6f} {cr_p1:<13.6f} {cr_p2:<14.6f} {am_ratio:<12.6f} {cr_ratio:<12.6f} {area:<10.4f}")
        sys.stdout.flush()

    print("=" * 115)

if __name__ == '__main__':
    run_aligned_resonance_sweep()
