"""
PROJECT JANUS: EXACT 1x2 SWITCH PHYSICS SOLVER (CLEAN OUTPUT)
"""
import sys
import os
import math
import numpy as np
import meep as mp
from meep import mpb

# Redirect MPB stdout to devnull during solve
def solve_single_wg_neff(w_wg, h_wg=0.220, t_pcm=0.0, n_pcm=2.70, resolution=40):
    lambda_0 = 1.064
    fcen = 1.0 / lambda_0
    n_si = 3.480
    n_sio2 = 1.449

    clad_mat = mp.Medium(index=n_sio2)
    core_mat = mp.Medium(index=n_si)
    pcm_mat = mp.Medium(index=n_pcm)

    sy = 2.0
    sz = 2.0
    cell = mp.Vector3(0, sy, sz)

    geometry = [
        mp.Block(center=mp.Vector3(0, 0, -sz/4.0), size=mp.Vector3(mp.inf, mp.inf, sz/2.0), material=clad_mat),
        mp.Block(center=mp.Vector3(0, 0, h_wg/2.0), size=mp.Vector3(mp.inf, w_wg, h_wg), material=core_mat),
    ]
    if t_pcm > 0:
        geometry.append(mp.Block(center=mp.Vector3(0, 0, h_wg + t_pcm/2.0), size=mp.Vector3(mp.inf, w_wg, t_pcm), material=pcm_mat))

    k_guess = 2.40 * fcen
    k_points = [mp.Vector3(k_guess, 0, 0)]

    ms = mpb.ModeSolver(
        geometry=geometry,
        geometry_lattice=mp.Lattice(size=cell),
        k_points=k_points,
        resolution=resolution,
        default_material=clad_mat,
        num_bands=2
    )
    # Silence MPB output
    old_stdout = sys.stdout
    sys.stdout = open(os.devnull, 'w')
    try:
        k_res = ms.find_k(mp.NO_PARITY, fcen, 1, 1, mp.Vector3(1, 0, 0), 1e-4, k_guess, 0.5*fcen, 3.48*fcen)
    finally:
        sys.stdout.close()
        sys.stdout = old_stdout
    return k_res[0] / fcen

def solve_coupler_supermodes(w1, w2, gap, h_wg=0.220, t_pcm=0.0, n_pcm=2.70, resolution=40):
    lambda_0 = 1.064
    fcen = 1.0 / lambda_0
    n_si = 3.480
    n_sio2 = 1.449

    clad_mat = mp.Medium(index=n_sio2)
    core_mat = mp.Medium(index=n_si)
    pcm_mat = mp.Medium(index=n_pcm)

    y_sep = (w1 + w2)/2.0 + gap
    y1 = y_sep / 2.0
    y2 = -y_sep / 2.0

    sy = 3.0
    sz = 2.0
    cell = mp.Vector3(0, sy, sz)

    geometry = [
        mp.Block(center=mp.Vector3(0, 0, -sz/4.0), size=mp.Vector3(mp.inf, mp.inf, sz/2.0), material=clad_mat),
        mp.Block(center=mp.Vector3(0, y1, h_wg/2.0), size=mp.Vector3(mp.inf, w1, h_wg), material=core_mat),
        mp.Block(center=mp.Vector3(0, y2, h_wg/2.0), size=mp.Vector3(mp.inf, w2, h_wg), material=core_mat),
    ]
    if t_pcm > 0:
        geometry.append(mp.Block(center=mp.Vector3(0, y1, h_wg + t_pcm/2.0), size=mp.Vector3(mp.inf, w1, t_pcm), material=pcm_mat))

    k_guess = 2.45 * fcen
    k_points = [mp.Vector3(k_guess, 0, 0)]

    ms = mpb.ModeSolver(
        geometry=geometry,
        geometry_lattice=mp.Lattice(size=cell),
        k_points=k_points,
        resolution=resolution,
        default_material=clad_mat,
        num_bands=4
    )
    old_stdout = sys.stdout
    sys.stdout = open(os.devnull, 'w')
    try:
        k_res = ms.find_k(mp.NO_PARITY, fcen, 1, 2, mp.Vector3(1, 0, 0), 1e-4, k_guess, 0.5*fcen, 3.48*fcen)
    finally:
        sys.stdout.close()
        sys.stdout = old_stdout
    
    n_even = max(k_res[0], k_res[1]) / fcen
    n_odd = min(k_res[0], k_res[1]) / fcen
    return n_even, n_odd

if __name__ == "__main__":
    print("=" * 75)
    print("PROJECT JANUS: HIGH-ACCURACY MODE & PHASE-MATCHING SOLVER")
    print("=" * 75)

    print("\n1. BARE SILICON WAVEGUIDES (h = 220 nm):")
    bare_widths = np.arange(0.180, 0.360, 0.020)
    bare_neff = {}
    for w in bare_widths:
        neff = solve_single_wg_neff(w, t_pcm=0.0)
        bare_neff[round(w, 3)] = neff
        print(f"  w = {w*1000:3.0f} nm -> n_eff = {neff:.5f}")

    print("\n2. ACTIVE WAVEGUIDE 1 (w1 = 200, 220, 240 nm with Sb2S3 PCM):")
    for w1 in [0.200, 0.220, 0.240]:
        print(f"\n  --- Waveguide 1 width = {w1*1000:.0f} nm ---")
        for t_pcm in [0.030, 0.050, 0.080, 0.100, 0.150]:
            n_am = solve_single_wg_neff(w1, t_pcm=t_pcm, n_pcm=2.70)
            n_cr = solve_single_wg_neff(w1, t_pcm=t_pcm, n_pcm=3.30)
            dn = n_cr - n_am
            # Find which bare waveguide width w2 matches n_am
            print(f"    t_pcm={t_pcm*1000:3.0f} nm | n_am={n_am:.5f} | n_cr={n_cr:.5f} | delta_n={dn:.5f}")

    print("=" * 75)
