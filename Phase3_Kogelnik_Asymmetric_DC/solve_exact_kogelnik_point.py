import sys, os, math
import numpy as np
import meep as mp
from meep import mpb

def solve_exact_kogelnik():
    lambda_0 = 1.064
    fcen = 1.0 / lambda_0
    clad = mp.Medium(index=1.449)
    core = mp.Medium(index=3.480)
    w_wg = 0.240
    h_wg = 0.220
    pcm_am = mp.Medium(index=2.70)
    pcm_cr = mp.Medium(index=3.30)

    sy, sz = 3.0, 2.5
    cell = mp.Vector3(0, sy, sz)

    print("=" * 80)
    print("  EXACT KOGELNIK OPERATING POINT SOLVER FOR 1x2 SWITCH")
    print("=" * 80)

    for t_pcm in [0.100, 0.120, 0.140, 0.150, 0.160, 0.180]:
        geom_am = [
            mp.Block(center=mp.Vector3(0, 0, -sz/4.0), size=mp.Vector3(mp.inf, mp.inf, sz/2.0), material=clad),
            mp.Block(center=mp.Vector3(0, 0, h_wg/2.0), size=mp.Vector3(mp.inf, w_wg, h_wg), material=core),
            mp.Block(center=mp.Vector3(0, 0, h_wg + t_pcm/2.0), size=mp.Vector3(mp.inf, w_wg, t_pcm), material=pcm_am),
        ]
        ms = mpb.ModeSolver(geometry=geom_am, geometry_lattice=mp.Lattice(size=cell),
                            k_points=[mp.Vector3(2.45*fcen, 0, 0)], resolution=40,
                            default_material=clad, num_bands=2)
        old_stdout = sys.stdout
        sys.stdout = open(os.devnull, 'w')
        k_am = ms.find_k(mp.NO_PARITY, fcen, 1, 1, mp.Vector3(1, 0, 0), 1e-4, 2.45*fcen, 0.5*fcen, 3.48*fcen)
        
        geom_cr = [
            mp.Block(center=mp.Vector3(0, 0, -sz/4.0), size=mp.Vector3(mp.inf, mp.inf, sz/2.0), material=clad),
            mp.Block(center=mp.Vector3(0, 0, h_wg/2.0), size=mp.Vector3(mp.inf, w_wg, h_wg), material=core),
            mp.Block(center=mp.Vector3(0, 0, h_wg + t_pcm/2.0), size=mp.Vector3(mp.inf, w_wg, t_pcm), material=pcm_cr),
        ]
        ms.geometry = geom_cr
        k_cr = ms.find_k(mp.NO_PARITY, fcen, 1, 1, mp.Vector3(1, 0, 0), 1e-4, 2.45*fcen, 0.5*fcen, 3.48*fcen)
        sys.stdout.close()
        sys.stdout = old_stdout

        nam = k_am[0] / fcen
        ncr = k_cr[0] / fcen
        dn_eff = ncr - nam
        delta_beta = (2.0 * math.pi / lambda_0) * dn_eff
        target_kappa = delta_beta / (2.0 * math.sqrt(3.0))
        target_lc = math.pi / (2.0 * target_kappa)

        print(f"\nt_pcm = {t_pcm*1000:.0f} nm: n_am = {nam:.4f}, n_cr = {ncr:.4f}, Delta_n_eff = {dn_eff:.4f}")
        print(f"  Delta_beta = {delta_beta:.4f} rad/um")
        print(f"  Target kappa for exact Kogelnik null: {target_kappa:.4f} rad/um")
        print(f"  Exact required L_c for both states:   {target_lc:.3f} um")
        cell_len = 0.35 + target_lc + 1.15 + 0.35
        cell_w = 0.52 + 0.30/2.0 + 0.24/2.0 + 2 * 0.12
        area = cell_len * cell_w
        print(f"  Cell Footprint: {cell_len:.2f} um x {cell_w:.2f} um = {area:.3f} um^2 (Limit <= 6.0 um^2: {'PASS' if area <= 6.0 else 'FAIL'})")

if __name__ == '__main__':
    solve_exact_kogelnik()
