import sys, os, math
import numpy as np
import meep as mp
from meep import mpb

lambda_0 = 1.064
fcen = 1.0 / lambda_0
h_wg = 0.220
clad_mat = mp.Medium(index=1.449)
core_mat = mp.Medium(index=3.480)

def solve_isolated_wg(w_wg, t_pcm, n_pcm, resolution=40):
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

def solve_supermodes(w_wg, gap, t_pcm, n_pcm=2.70, resolution=40):
    pcm_mat = mp.Medium(index=n_pcm)
    y_sep = w_wg + gap
    y1 = y_sep / 2.0
    y2 = -y_sep / 2.0
    sy, sz = 3.0, 2.5
    cell = mp.Vector3(0, sy, sz)
    geom = [
        mp.Block(center=mp.Vector3(0, 0, -sz/4.0), size=mp.Vector3(mp.inf, mp.inf, sz/2.0), material=clad_mat),
        mp.Block(center=mp.Vector3(0, y1, h_wg/2.0), size=mp.Vector3(mp.inf, w_wg, h_wg), material=core_mat),
        mp.Block(center=mp.Vector3(0, y2, h_wg/2.0), size=mp.Vector3(mp.inf, w_wg, h_wg), material=core_mat),
        mp.Block(center=mp.Vector3(0, y1, h_wg + t_pcm/2.0), size=mp.Vector3(mp.inf, w_wg, t_pcm), material=pcm_mat),
        mp.Block(center=mp.Vector3(0, y2, h_wg + t_pcm/2.0), size=mp.Vector3(mp.inf, w_wg, t_pcm), material=pcm_mat),
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
    dn_super = n_even - n_odd
    l_pi = lambda_0 / (2.0 * max(dn_super, 1e-6))
    kappa = math.pi / (2.0 * l_pi)
    return n_even, n_odd, dn_super, l_pi, kappa

if __name__ == "__main__":
    for w_wg in [0.240, 0.220, 0.200]:
        t_pcm = 0.150
        n_am = solve_isolated_wg(w_wg, t_pcm, 2.70)
        n_cr = solve_isolated_wg(w_wg, t_pcm, 3.30)
        dn_eff = n_cr - n_am
        target_dn_super = dn_eff / math.sqrt(3.0)
        ideal_lc = (math.sqrt(3.0) * lambda_0) / (2.0 * dn_eff)
        print("="*75)
        print(f"ISOLATED WG: w={w_wg*1000:.0f} nm, t_pcm={t_pcm*1000:.0f} nm")
        print(f"  n_am = {n_am:.5f}, n_cr = {n_cr:.5f} -> dn_eff = {dn_eff:.5f}")
        print(f"  TARGET SUPERMODE SPLIT: dn_super = {target_dn_super:.5f}")
        print(f"  IDEAL COUPLING LENGTH : L_c = {ideal_lc:.3f} um")
        print("="*75)
        print(f"{'Gap (nm)':<10}{'n_even':<10}{'n_odd':<10}{'dn_super':<12}{'L_pi (um)':<12}{'Ratio (dn_s/tgt)':<15}")
        for g_nm in range(50, 150, 10):
            gap = g_nm * 1e-3
            ne, no, dns, lpi, kap = solve_supermodes(w_wg, gap, t_pcm)
            ratio = dns / target_dn_super
            print(f"{g_nm:<10}{ne:<10.5f}{no:<10.5f}{dns:<12.5f}{lpi:<12.3f}{ratio:<15.3f}")
