"""
PROJECT JANUS: Si3N4 PHASE-MISMATCHED DIRECTIONAL COUPLER OPERATING POINT FINDER
================================================================================
Calculates exact MPB vector supermodes at lambda_0 = 1064 nm on Si3N4 (n = 2.016)
clad in SiO2 (n = 1.449) with Sb2S3 PCM (n_am = 2.70, n_cr = 3.30).
"""

import sys, os, math
import numpy as np
import meep as mp
from meep import mpb

def solve_sin_mode(w=0.450, h=0.300, t_pcm=0.030, n_pcm=2.70, resolution=40):
    lambda_0 = 1.064
    fcen = 1.0 / lambda_0
    clad = mp.Medium(index=1.449)
    core = mp.Medium(index=2.016)
    pcm = mp.Medium(index=n_pcm)

    sy, sz = 3.0, 3.0
    cell = mp.Vector3(0, sy, sz)
    geom = [
        mp.Block(center=mp.Vector3(0, 0, -sz/4.0), size=mp.Vector3(mp.inf, mp.inf, sz/2.0), material=clad),
        mp.Block(center=mp.Vector3(0, 0, h/2.0), size=mp.Vector3(mp.inf, w, h), material=core),
    ]
    if t_pcm > 1e-4:
        geom.append(mp.Block(center=mp.Vector3(0, 0, h + t_pcm/2.0), size=mp.Vector3(mp.inf, w, t_pcm), material=pcm))

    k_guess = 1.70 * fcen
    ms = mpb.ModeSolver(geometry=geom, geometry_lattice=mp.Lattice(size=cell),
                        k_points=[mp.Vector3(k_guess, 0, 0)], resolution=resolution,
                        default_material=clad, num_bands=2)
    old_stdout = sys.stdout
    sys.stdout = open(os.devnull, 'w')
    try:
        k_res = ms.find_k(mp.NO_PARITY, fcen, 1, 1, mp.Vector3(1, 0, 0), 1e-4, k_guess, 0.5*fcen, 2.016*fcen)
    finally:
        sys.stdout.close()
        sys.stdout = old_stdout
    return k_res[0] / fcen

def solve_sin_supermodes(w=0.450, gap=0.100, h=0.300, t_pcm=0.030, n_pcm=2.70, resolution=40):
    lambda_0 = 1.064
    fcen = 1.0 / lambda_0
    clad = mp.Medium(index=1.449)
    core = mp.Medium(index=2.016)
    pcm = mp.Medium(index=n_pcm)

    y_sep = w + gap
    y1, y2 = y_sep / 2.0, -y_sep / 2.0
    sy, sz = 4.0, 3.0
    cell = mp.Vector3(0, sy, sz)
    geom = [
        mp.Block(center=mp.Vector3(0, 0, -sz/4.0), size=mp.Vector3(mp.inf, mp.inf, sz/2.0), material=clad),
        mp.Block(center=mp.Vector3(0, y1, h/2.0), size=mp.Vector3(mp.inf, w, h), material=core),
        mp.Block(center=mp.Vector3(0, y2, h/2.0), size=mp.Vector3(mp.inf, w, h), material=core),
    ]
    if t_pcm > 1e-4:
        geom.append(mp.Block(center=mp.Vector3(0, y1, h + t_pcm/2.0), size=mp.Vector3(mp.inf, w, t_pcm), material=pcm))
        geom.append(mp.Block(center=mp.Vector3(0, y2, h + t_pcm/2.0), size=mp.Vector3(mp.inf, w, t_pcm), material=pcm))

    k_guess = 1.70 * fcen
    ms = mpb.ModeSolver(geometry=geom, geometry_lattice=mp.Lattice(size=cell),
                        k_points=[mp.Vector3(k_guess, 0, 0)], resolution=resolution,
                        default_material=clad, num_bands=4)
    old_stdout = sys.stdout
    sys.stdout = open(os.devnull, 'w')
    try:
        k_res = ms.find_k(mp.NO_PARITY, fcen, 1, 2, mp.Vector3(1, 0, 0), 1e-4, k_guess, 0.5*fcen, 2.016*fcen)
    finally:
        sys.stdout.close()
        sys.stdout = old_stdout
    ne = max(k_res[0], k_res[1]) / fcen
    no = min(k_res[0], k_res[1]) / fcen
    return ne, no

if __name__ == '__main__':
    w = 0.450
    h = 0.300
    lambda_0 = 1.064
    print("=" * 80)
    print("Si3N4 DIRECTIONAL COUPLER EXACT OPERATING POINT FINDER")
    print("=" * 80)

    for t_pcm in [0.020, 0.030, 0.040]:
        print(f"\n--- Testing Sb2S3 Patch Thickness t_pcm = {t_pcm*1000:.0f} nm ---")
        n_am = solve_sin_mode(w, h, t_pcm=t_pcm, n_pcm=2.70)
        n_cr = solve_sin_mode(w, h, t_pcm=t_pcm, n_pcm=3.30)
        dn_eff = n_cr - n_am
        target_split = dn_eff / math.sqrt(3.0)
        optimal_lc = math.sqrt(3.0) * lambda_0 / (2.0 * dn_eff)

        print(f"  n_eff(Amorphous)   = {n_am:.5f}")
        print(f"  n_eff(Crystalline) = {n_cr:.5f}")
        print(f"  delta_n_eff        = {dn_eff:.5f}")
        print(f"  Target Splitting   = {target_split:.5f}")
        print(f"  Target Coupling Lc = {optimal_lc:.3f} um")

        print(f"\n{'Gap (nm)':<10} {'n_even':<10} {'n_odd':<10} {'dn_split':<10} {'L_c (um)':<10} {'u':<8} {'P_bar(Cr)%':<12} {'XT(Cr) dB':<12}")
        print("-" * 82)
        for g in [0.060, 0.080, 0.100, 0.120, 0.140, 0.160]:
            ne, no = solve_sin_supermodes(w, g, h, t_pcm=t_pcm, n_pcm=2.70)
            split = ne - no
            lc = lambda_0 / (2.0 * max(split, 1e-6))
            u = dn_eff / max(split, 1e-6)
            arg = 0.5 * math.pi * math.sqrt(1.0 + u**2)
            p_cross_cr = (1.0 / (1.0 + u**2)) * (math.sin(arg)**2)
            p_bar_cr = 1.0 - p_cross_cr
            xt_cr = 10.0 * math.log10(max(p_cross_cr / max(p_bar_cr, 1e-9), 1e-9))
            print(f"{g*1000:<10.1f} {ne:<10.5f} {no:<10.5f} {split:<10.5f} {lc:<10.3f} {u:<8.3f} {p_bar_cr*100:<11.2f}% {xt_cr:<12.2f}")
