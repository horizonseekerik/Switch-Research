"""
Find ultra-compact phase-matched designs for 1x2 switch
"""
import sys, os, math
import numpy as np
import meep as mp
from meep import mpb

def solve_wg_neff(w_wg, h_wg=0.220, t_pcm=0.0, n_pcm=2.70, resolution=40):
    lambda_0 = 1.064
    fcen = 1.0 / lambda_0
    clad_mat = mp.Medium(index=1.449)
    core_mat = mp.Medium(index=3.480)
    pcm_mat = mp.Medium(index=n_pcm)
    sy, sz = 2.5, 2.5
    cell = mp.Vector3(0, sy, sz)
    geom = [
        mp.Block(center=mp.Vector3(0, 0, -sz/4.0), size=mp.Vector3(mp.inf, mp.inf, sz/2.0), material=clad_mat),
        mp.Block(center=mp.Vector3(0, 0, h_wg/2.0), size=mp.Vector3(mp.inf, w_wg, h_wg), material=core_mat),
    ]
    if t_pcm > 0:
        geom.append(mp.Block(center=mp.Vector3(0, 0, h_wg + t_pcm/2.0), size=mp.Vector3(mp.inf, w_wg, t_pcm), material=pcm_mat))
    k_guess = 2.30 * fcen
    ms = mpb.ModeSolver(geometry=geom, geometry_lattice=mp.Lattice(size=cell),
                        k_points=[mp.Vector3(k_guess, 0, 0)], resolution=resolution,
                        default_material=clad_mat, num_bands=2)
    old_stdout = sys.stdout
    sys.stdout = open(os.devnull, 'w')
    try:
        k_res = ms.find_k(mp.NO_PARITY, fcen, 1, 1, mp.Vector3(1, 0, 0), 1e-4, k_guess, 0.5*fcen, 3.48*fcen)
    finally:
        sys.stdout.close()
        sys.stdout = old_stdout
    return k_res[0] / fcen

if __name__ == "__main__":
    print("Testing narrower w1 and thicker t_pcm to maximize delta_n and minimize L_c:")
    for w1 in [0.160, 0.180, 0.200]:
        for t_pcm in [0.150, 0.180, 0.220]:
            n_am = solve_wg_neff(w1, t_pcm=t_pcm, n_pcm=2.70)
            n_cr = solve_wg_neff(w1, t_pcm=t_pcm, n_pcm=3.30)
            dn = n_cr - n_am
            lc_opt = 0.9214 / dn
            print(f"  w1={w1*1000:3.0f}nm, t_pcm={t_pcm*1000:3.0f}nm -> n_am={n_am:.5f}, n_cr={n_cr:.5f}, dn={dn:.5f} -> L_c(ideal)={lc_opt:.2f} um")
