"""
Check supermode coupling for identical waveguides with Sb2S3 PCM (amorphous)
"""
import sys, os, math
import numpy as np
import meep as mp
from meep import mpb

def solve_symmetric_coupler(w, gap, t_pcm, n_pcm=2.70, resolution=40):
    lambda_0 = 1.064
    fcen = 1.0 / lambda_0
    h_wg = 0.220
    clad_mat = mp.Medium(index=1.449)
    core_mat = mp.Medium(index=3.480)
    pcm_mat = mp.Medium(index=n_pcm)

    y_sep = w + gap
    y1 = y_sep / 2.0
    y2 = -y_sep / 2.0

    sy, sz = 3.0, 2.5
    cell = mp.Vector3(0, sy, sz)
    geom = [
        mp.Block(center=mp.Vector3(0, 0, -sz/4.0), size=mp.Vector3(mp.inf, mp.inf, sz/2.0), material=clad_mat),
        mp.Block(center=mp.Vector3(0, y1, h_wg/2.0), size=mp.Vector3(mp.inf, w, h_wg), material=core_mat),
        mp.Block(center=mp.Vector3(0, y2, h_wg/2.0), size=mp.Vector3(mp.inf, w, h_wg), material=core_mat),
        mp.Block(center=mp.Vector3(0, y1, h_wg + t_pcm/2.0), size=mp.Vector3(mp.inf, w, t_pcm), material=pcm_mat),
        mp.Block(center=mp.Vector3(0, y2, h_wg + t_pcm/2.0), size=mp.Vector3(mp.inf, w, t_pcm), material=pcm_mat),
    ]

    k_guess = 2.50 * fcen
    ms = mpb.ModeSolver(geometry=geom, geometry_lattice=mp.Lattice(size=cell),
                        k_points=[mp.Vector3(k_guess, 0, 0)], resolution=resolution,
                        default_material=clad_mat, num_bands=4)
    old_stdout = sys.stdout
    sys.stdout = open(os.devnull, 'w')
    try:
        k_res = ms.find_k(mp.NO_PARITY, fcen, 1, 2, mp.Vector3(1, 0, 0), 1e-4, k_guess, 0.5*fcen, 3.48*fcen)
    finally:
        sys.stdout.close()
        sys.stdout = old_stdout

    n_even = max(k_res[0], k_res[1]) / fcen
    n_odd = min(k_res[0], k_res[1]) / fcen
    dn = n_even - n_odd
    l_c = lambda_0 / (2.0 * max(dn, 1e-6))
    return n_even, n_odd, dn, l_c

if __name__ == "__main__":
    w = 0.240
    t_pcm = 0.150
    print(f"Sweeping gap for symmetric coupler (w={w*1000:.0f}nm, t_pcm={t_pcm*1000:.0f}nm, Amorphous):")
    for gap in [0.030, 0.035, 0.040, 0.045, 0.050, 0.055, 0.060]:
        ne, no, dn, lc = solve_symmetric_coupler(w, gap, t_pcm)
        print(f"  gap = {gap*1000:2.0f} nm -> n_even = {ne:.5f}, n_odd = {no:.5f}, dn = {dn:.5f} -> L_c = {lc:.3f} um")
