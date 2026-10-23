import sys, os, math
import numpy as np
import meep as mp
from meep import mpb

lambda_0 = 1.064
fcen = 1.0 / lambda_0
h_wg = 0.220
clad_mat = mp.Medium(index=1.449)
core_mat = mp.Medium(index=3.480)

def solve_wg(w_wg, t_pcm, n_pcm, resolution=32):
    pcm_mat = mp.Medium(index=n_pcm)
    sy, sz = 2.5, 2.5
    cell = mp.Vector3(0, sy, sz)
    geom = [
        mp.Block(center=mp.Vector3(0, 0, -sz/4.0), size=mp.Vector3(mp.inf, mp.inf, sz/2.0), material=clad_mat),
        mp.Block(center=mp.Vector3(0, 0, h_wg/2.0), size=mp.Vector3(mp.inf, w_wg, h_wg), material=core_mat),
        mp.Block(center=mp.Vector3(0, 0, h_wg + t_pcm/2.0), size=mp.Vector3(mp.inf, w_wg, t_pcm), material=pcm_mat)
    ]
    k_guess = 2.40 * fcen
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

print("="*105)
print(f"{'w_core(nm)':<12} {'t_pcm(nm)':<12} {'n_am':<14} {'n_cr':<14} {'dn_eff':<14} {'ideal_Lc(um)':<16} {'Delta_beta(rad/um)':<20}")
print("="*105)

for w_nm in [180, 200, 220, 240]:
    w = w_nm * 1e-3
    for t_nm in [200, 240, 280, 320]:
        t = t_nm * 1e-3
        nam = solve_wg(w, t, 2.70)
        ncr = solve_wg(w, t, 3.30)
        dn = ncr - nam
        dbeta = (2.0 * math.pi / lambda_0) * dn
        lc_kogelnik = (math.sqrt(3.0) * lambda_0) / (2.0 * dn)
        print(f"{w_nm:<12} {t_nm:<12} {nam:<14.10f} {ncr:<14.10f} {dn:<14.10f} {lc_kogelnik:<16.10f} {dbeta:<20.10f}")
