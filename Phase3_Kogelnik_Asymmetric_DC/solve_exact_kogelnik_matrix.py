import sys, os, math
import numpy as np
import meep as mp
from meep import mpb

def analyze():
    lambda_0 = 1.064
    fcen = 1.0 / lambda_0
    clad = mp.Medium(index=1.449)
    core = mp.Medium(index=3.480)
    w_wg = 0.240
    h_wg = 0.220
    pcm_am = mp.Medium(index=2.70)
    pcm_cr = mp.Medium(index=3.30)

    sy, sz = 4.0, 3.0
    cell = mp.Vector3(0, sy, sz)

    print("=" * 85)
    print("  PROJECT JANUS: MPB VECTOR SUPERMODE DISPERSION MAPPING")
    print("=" * 85)

    for t_pcm in [0.160, 0.180, 0.200, 0.220]:
        geom_am = [
            mp.Block(center=mp.Vector3(0, 0, -sz/4.0), size=mp.Vector3(mp.inf, mp.inf, sz/2.0), material=clad),
            mp.Block(center=mp.Vector3(0, 0, h_wg/2.0), size=mp.Vector3(mp.inf, w_wg, h_wg), material=core),
            mp.Block(center=mp.Vector3(0, 0, h_wg + t_pcm/2.0), size=mp.Vector3(mp.inf, w_wg, t_pcm), material=pcm_am),
        ]
        ms = mpb.ModeSolver(geometry=geom_am, geometry_lattice=mp.Lattice(size=cell),
                            k_points=[mp.Vector3(2.5*fcen, 0, 0)], resolution=40,
                            default_material=clad, num_bands=2)
        old_stdout = sys.stdout
        sys.stdout = open(os.devnull, 'w')
        k_am = ms.find_k(mp.NO_PARITY, fcen, 1, 1, mp.Vector3(1, 0, 0), 1e-4, 2.5*fcen, 0.5*fcen, 3.48*fcen)
        
        geom_cr = [
            mp.Block(center=mp.Vector3(0, 0, -sz/4.0), size=mp.Vector3(mp.inf, mp.inf, sz/2.0), material=clad),
            mp.Block(center=mp.Vector3(0, 0, h_wg/2.0), size=mp.Vector3(mp.inf, w_wg, h_wg), material=core),
            mp.Block(center=mp.Vector3(0, 0, h_wg + t_pcm/2.0), size=mp.Vector3(mp.inf, w_wg, t_pcm), material=pcm_cr),
        ]
        ms.geometry = geom_cr
        k_cr = ms.find_k(mp.NO_PARITY, fcen, 1, 1, mp.Vector3(1, 0, 0), 1e-4, 2.5*fcen, 0.5*fcen, 3.48*fcen)
        sys.stdout.close()
        sys.stdout = old_stdout

        nam = k_am[0] / fcen
        ncr = k_cr[0] / fcen
        dn = ncr - nam
        delta_beta = (2.0 * math.pi / lambda_0) * dn

        print(f"\n--- t_pcm = {t_pcm*1000:.0f} nm: n_am = {nam:.4f}, n_cr = {ncr:.4f}, dn = {dn:.4f}, delta_beta = {delta_beta:.4f} rad/um ---")
        print(f"{'Gap (nm)':<10} {'dn_split':<10} {'kappa':<10} {'L_pi/2 (um)':<12} {'u = dbeta/(2k)':<16} {'Match to sqrt(3)?':<20}")
        print("-" * 80)

        for gap in [0.090, 0.100, 0.110, 0.120, 0.130, 0.140]:
            y_sep = w_wg + gap
            y1, y2 = y_sep / 2.0, -y_sep / 2.0
            geom_sym = [
                mp.Block(center=mp.Vector3(0, 0, -sz/4.0), size=mp.Vector3(mp.inf, mp.inf, sz/2.0), material=clad),
                mp.Block(center=mp.Vector3(0, y1, h_wg/2.0), size=mp.Vector3(mp.inf, w_wg, h_wg), material=core),
                mp.Block(center=mp.Vector3(0, y2, h_wg/2.0), size=mp.Vector3(mp.inf, w_wg, h_wg), material=core),
                mp.Block(center=mp.Vector3(0, y1, h_wg + t_pcm/2.0), size=mp.Vector3(mp.inf, w_wg, t_pcm), material=pcm_am),
                mp.Block(center=mp.Vector3(0, y2, h_wg + t_pcm/2.0), size=mp.Vector3(mp.inf, w_wg, t_pcm), material=pcm_am),
            ]
            ms = mpb.ModeSolver(geometry=geom_sym, geometry_lattice=mp.Lattice(size=cell),
                                k_points=[mp.Vector3(2.5*fcen, 0, 0)], resolution=40,
                                default_material=clad, num_bands=4)
            old_stdout = sys.stdout
            sys.stdout = open(os.devnull, 'w')
            k_sym = ms.find_k(mp.NO_PARITY, fcen, 1, 2, mp.Vector3(1, 0, 0), 1e-4, 2.5*fcen, 0.5*fcen, 3.48*fcen)
            sys.stdout.close()
            sys.stdout = old_stdout

            ne = max(k_sym[0], k_sym[1]) / fcen
            no = min(k_sym[0], k_sym[1]) / fcen
            dn_split = ne - no
            kappa = (math.pi / lambda_0) * dn_split
            lc = math.pi / (2.0 * kappa)
            u = delta_beta / (2.0 * kappa)
            match_str = f"EXACT MATCH!" if abs(u - math.sqrt(3)) < 0.05 else f"diff = {u - math.sqrt(3):+.3f}"

            print(f"{gap*1000:<10.0f} {dn_split:<10.5f} {kappa:<10.4f} {lc:<12.3f} {u:<16.3f} {match_str:<20}")

if __name__ == '__main__':
    analyze()
