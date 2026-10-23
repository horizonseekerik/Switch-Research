"""
FIND EXACT GAP FOR EACH PCM THICKNESS TO SATISFY Lc0 == Lc1 (DELTA L = 0)
========================================================================
For each t_pcm, finds the coupling gap where State 0 half-beat (Lc0) equals
State 1 full-beat (Lc1).
"""

import math
import numpy as np
import meep as mp
from meep import mpb

fcen = 1.0 / 1.064
clad = mp.Medium(index=1.449)
core = mp.Medium(index=3.480)
h_wg = 0.220
w = 0.240
cell = mp.Vector3(0, 4.0, 2.5)

def solve_coupler(gap, t_pcm):
    d = gap + w
    # State 0 (both amorphous)
    geom_0 = [
        mp.Block(center=mp.Vector3(0, 0, -0.6), size=mp.Vector3(mp.inf, mp.inf, 1.2), material=clad),
        mp.Block(center=mp.Vector3(0, +d/2.0, h_wg/2.0), size=mp.Vector3(mp.inf, w, h_wg), material=core),
        mp.Block(center=mp.Vector3(0, +d/2.0, h_wg + t_pcm/2.0), size=mp.Vector3(mp.inf, w, t_pcm), material=mp.Medium(index=2.70)),
        mp.Block(center=mp.Vector3(0, -d/2.0, h_wg/2.0), size=mp.Vector3(mp.inf, w, h_wg), material=core),
        mp.Block(center=mp.Vector3(0, -d/2.0, h_wg + t_pcm/2.0), size=mp.Vector3(mp.inf, w, t_pcm), material=mp.Medium(index=2.70)),
    ]
    ms = mpb.ModeSolver(geometry=geom_0, geometry_lattice=mp.Lattice(size=cell),
                        k_points=[mp.Vector3(2.2*fcen, 0, 0)], resolution=32,
                        default_material=clad, num_bands=2, verbose=False)
    k_e0 = ms.find_k(mp.NO_PARITY, fcen, 1, 1, mp.Vector3(1, 0, 0), 1e-4, 2.0*fcen, 0.5*fcen, 3.48*fcen)[0]
    k_o0 = ms.find_k(mp.NO_PARITY, fcen, 2, 2, mp.Vector3(1, 0, 0), 1e-4, 2.0*fcen, 0.5*fcen, 3.48*fcen)[0]
    lc_0 = 1.0 / (2.0 * (k_e0 - k_o0))

    # State 1 (WG1 crystalline, WG2 amorphous)
    geom_1 = [
        mp.Block(center=mp.Vector3(0, 0, -0.6), size=mp.Vector3(mp.inf, mp.inf, 1.2), material=clad),
        mp.Block(center=mp.Vector3(0, +d/2.0, h_wg/2.0), size=mp.Vector3(mp.inf, w, h_wg), material=core),
        mp.Block(center=mp.Vector3(0, +d/2.0, h_wg + t_pcm/2.0), size=mp.Vector3(mp.inf, w, t_pcm), material=mp.Medium(index=3.30)),
        mp.Block(center=mp.Vector3(0, -d/2.0, h_wg/2.0), size=mp.Vector3(mp.inf, w, h_wg), material=core),
        mp.Block(center=mp.Vector3(0, -d/2.0, h_wg + t_pcm/2.0), size=mp.Vector3(mp.inf, w, t_pcm), material=mp.Medium(index=2.70)),
    ]
    ms.geometry = geom_1
    k_e1 = ms.find_k(mp.NO_PARITY, fcen, 1, 1, mp.Vector3(1, 0, 0), 1e-4, 2.0*fcen, 0.5*fcen, 3.48*fcen)[0]
    k_o1 = ms.find_k(mp.NO_PARITY, fcen, 2, 2, mp.Vector3(1, 0, 0), 1e-4, 2.0*fcen, 0.5*fcen, 3.48*fcen)[0]
    lc_1 = 1.0 / (k_e1 - k_o1)

    return lc_0, lc_1

def find_zero_diff_gap(t_pcm_nm):
    t_pcm = t_pcm_nm * 1e-3
    # Bisection search over gap in [0.060, 0.160] um
    g_low = 0.060
    g_high = 0.160

    lc0_l, lc1_l = solve_coupler(g_low, t_pcm)
    diff_l = lc0_l - lc1_l

    lc0_h, lc1_h = solve_coupler(g_high, t_pcm)
    diff_h = lc0_h - lc1_h

    print(f"t_pcm={t_pcm_nm}nm: at {g_low*1e3:.1f}nm diff={diff_l*1e3:+.1f}nm, at {g_high*1e3:.1f}nm diff={diff_h*1e3:+.1f}nm")

    # If diff doesn't cross zero, report
    if diff_l * diff_h > 0:
        return None, None, None

    for _ in range(8):
        g_mid = (g_low + g_high) / 2.0
        lc0_m, lc1_m = solve_coupler(g_mid, t_pcm)
        diff_m = lc0_m - lc1_m
        if diff_l * diff_m <= 0:
            g_high = g_mid
            diff_h = diff_m
        else:
            g_low = g_mid
            diff_l = diff_m

    g_opt = (g_low + g_high) / 2.0
    lc0, lc1 = solve_coupler(g_opt, t_pcm)
    return g_opt, lc0, lc1

if __name__ == "__main__":
    print("=" * 85)
    print(f"{'t_pcm (nm)':<12} {'Opt Gap (nm)':<14} {'Lc0 (um)':<12} {'Lc1 (um)':<12} {'Diff (nm)':<12} {'Area (um^2)':<14}")
    print("-" * 85)
    for t_nm in [80.0, 100.0, 120.0, 140.0, 160.0, 190.0, 237.5]:
        g_opt, lc0, lc1 = find_zero_diff_gap(t_nm)
        if g_opt is not None:
            diff = (lc0 - lc1) * 1e3
            area = (0.35 + lc0 + 0.85 + 0.35) * (0.45 + 0.24 + 2*0.08)
            print(f"{t_nm:<12.1f} {g_opt*1e3:<14.2f} {lc0:<12.4f} {lc1:<12.4f} {diff:<12.2f} {area:<14.4f}")
        else:
            print(f"{t_nm:<12.1f} {'No zero-crossing in [60, 160]nm':<50}")
    print("=" * 85)
