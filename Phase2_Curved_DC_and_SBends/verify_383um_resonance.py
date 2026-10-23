"""
PROJECT JANUS: 3D FDTD MEEP SIMULATION AT EXACT KOGELNIK POINT (L_c ~ 3.83 um, Gap = 94.5 nm)
=============================================================================================
Analytical Derivation:
  Delta_beta = 1.4190771803 rad/um
  kappa_exact = Delta_beta / (2 * sqrt(3)) = 0.4096522960 rad/um
  L_c_exact = pi / (2 * kappa_exact) = 3.8344640000 um
  Gap for kappa = 0.40965 rad/um: 94.5000000000 nm.
  At this exact point:
    State 0 (Amorphous): kappa * L_c = pi/2 -> P2 = 100.0000000000%, P1 = 0.0000000000%
    State 1 (Crystalline): S * L_c = pi -> P1 = 100.0000000000%, P2 = 0.0000000000%
"""

import sys, os, math, time
import numpy as np
import meep as mp
from test_decoupled_transitions_dc import simulate_decoupled_dc

def run_kogelnik_exact_sweep():
    gap = 0.0945000000
    t_pcm = 0.2000000000
    res = 15

    print("=" * 115)
    print("  PROJECT JANUS: 3D FDTD VERIFICATION AT EXACT KOGELNIK RESONANCE POINT")
    print(f"  Parameters: Gap = {gap*1e3:.2f} nm, t_pcm = {t_pcm*1e3:.1f} nm, Res = {res}")
    print("=" * 115)
    print(f"{'L_c(um)':<8} {'AM P2(Cross)%':<16} {'AM P1(Bar)%':<14} {'CR P1(Bar)%':<14} {'CR P2(Cross)%':<16} {'AM Sig/G%':<12} {'CR Sig/G%':<12} {'Area(um2)':<10}")
    print("-" * 115)

    res_ref = simulate_decoupled_dc("amorphous", resolution=res, l_c=3.83, gap=gap, t_pcm=t_pcm, verbose=False)
    p_ref = res_ref["p_ref"]

    for lc in [3.70, 3.78, 3.84, 3.90]:
        res_am = simulate_decoupled_dc("amorphous", resolution=res, l_c=lc, gap=gap, t_pcm=t_pcm, p_ref=p_ref, verbose=False)
        res_cr = simulate_decoupled_dc("crystalline", resolution=res, l_c=lc, gap=gap, t_pcm=t_pcm, p_ref=p_ref, verbose=False)

        am_p2 = res_am["p2"] * 100.0
        am_p1 = res_am["p1"] * 100.0
        cr_p1 = res_cr["p1"] * 100.0
        cr_p2 = res_cr["p2"] * 100.0

        am_tot = am_p2 + am_p1
        cr_tot = cr_p1 + cr_p2

        am_ratio = (am_p2 / am_tot) * 100.0 if am_tot > 0 else 0.0
        cr_ratio = (cr_p1 / cr_tot) * 100.0 if cr_tot > 0 else 0.0
        area = res_am["cell_area"]

        print(f"{lc:<8.2f} {am_p2:<16.6f} {am_p1:<14.6f} {cr_p1:<14.6f} {cr_p2:<16.6f} {am_ratio:<12.6f} {cr_ratio:<12.6f} {area:<10.4f}")
        sys.stdout.flush()

    print("=" * 115)

if __name__ == '__main__':
    run_kogelnik_exact_sweep()
