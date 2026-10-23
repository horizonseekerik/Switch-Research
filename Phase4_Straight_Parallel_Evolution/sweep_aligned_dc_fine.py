import sys, os, math, time
import numpy as np
import meep as mp
from test_aligned_detuned_switch import simulate_aligned_switch

def run_fine_aligned_sweep():
    gap = 0.089
    t_pcm = 0.200
    res = 15

    print("=" * 115)
    print("  PROJECT JANUS: HIGH-PRECISION DUAL-STATE 3D FDTD SWEEP")
    print(f"  Configuration: Gap = {gap*1000:.0f} nm, t_pcm = {t_pcm*1000:.0f} nm, Res = {res}")
    print("  Design: Straight continuous bus, widened S-bend (w_bend = 360 nm), synchronous 0.25 um tapers")
    print("=" * 115)
    print(f"{'L_c(um)':<8} {'AM P2(Mod)%':<14} {'AM P1(Mod)%':<13} {'CR P1(Mod)%':<13} {'CR P2(Mod)%':<14} {'AM Sig/G%':<12} {'CR Sig/G%':<12} {'AM XT(dB)':<11} {'CR XT(dB)':<11} {'Area(um2)':<9}")
    print("-" * 115)

    res_ref = simulate_aligned_switch("amorphous", resolution=res, l_c=3.62, gap=gap, t_pcm=t_pcm, verbose=False)
    p_ref_flux = res_ref["p_ref_flux"]
    p_ref_mode = res_ref["p_ref_mode"]

    for lc in [3.56, 3.62, 3.68]:
        t0 = time.time()
        res_am = simulate_aligned_switch("amorphous", resolution=res, l_c=lc, gap=gap, t_pcm=t_pcm,
                                         p_ref_flux=p_ref_flux, p_ref_mode=p_ref_mode, verbose=False)
        res_cr = simulate_aligned_switch("crystalline", resolution=res, l_c=lc, gap=gap, t_pcm=t_pcm,
                                         p_ref_flux=p_ref_flux, p_ref_mode=p_ref_mode, verbose=False)

        am_p2 = res_am["p2_mode"] * 100.0
        am_p1 = res_am["p1_mode"] * 100.0
        cr_p1 = res_cr["p1_mode"] * 100.0
        cr_p2 = res_cr["p2_mode"] * 100.0

        am_ratio = res_am["ratio_guided"]
        cr_ratio = res_cr["ratio_guided"]
        am_xt = res_am["xt_db"]
        cr_xt = res_cr["xt_db"]
        area = res_am["cell_area"]
        elapsed = time.time() - t0

        print(f"{lc:<8.2f} {am_p2:<14.2f} {am_p1:<13.2f} {cr_p1:<13.2f} {cr_p2:<14.2f} {am_ratio:<12.2f} {cr_ratio:<12.2f} {am_xt:<11.2f} {cr_xt:<11.2f} {area:<9.3f}")
        sys.stdout.flush()

    print("=" * 115)

if __name__ == '__main__':
    run_fine_aligned_sweep()
