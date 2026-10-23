"""
PROJECT JANUS: SLOT COUPLER WITH PCM IN THE GAP
Find gap width where kappa_cr / kappa_am = 2.000 exactly!
"""
import sys, os, math
import numpy as np
import meep as mp
from meep import mpb

def solve_slot_coupler(w_si=0.240, h_si=0.220, gap=0.060, resolution=40):
    lambda_0 = 1.064
    fcen = 1.0 / lambda_0
    clad_mat = mp.Medium(index=1.449)
    si_mat = mp.Medium(index=3.480)

    y_sep = w_si + gap
    y1 = y_sep / 2.0
    y2 = -y_sep / 2.0

    sy, sz = 3.0, 2.5
    cell = mp.Vector3(0, sy, sz)

    results = {}
    for state, n_pcm in [("amorphous", 2.70), ("crystalline", 3.30)]:
        pcm_mat = mp.Medium(index=n_pcm)
        geom = [
            mp.Block(center=mp.Vector3(0, 0, -sz/4.0), size=mp.Vector3(mp.inf, mp.inf, sz/2.0), material=clad_mat),
            # WG1 and WG2
            mp.Block(center=mp.Vector3(0, y1, h_si/2.0), size=mp.Vector3(mp.inf, w_si, h_si), material=si_mat),
            mp.Block(center=mp.Vector3(0, y2, h_si/2.0), size=mp.Vector3(mp.inf, w_si, h_si), material=si_mat),
            # PCM inside the gap (slot) between WG1 and WG2
            mp.Block(center=mp.Vector3(0, 0, h_si/2.0), size=mp.Vector3(mp.inf, gap, h_si), material=pcm_mat),
        ]
        k_guess = 2.60 * fcen
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

        ne = max(k_res[0], k_res[1]) / fcen
        no = min(k_res[0], k_res[1]) / fcen
        dn = ne - no
        kappa = math.pi * dn / lambda_0
        lc = lambda_0 / (2.0 * max(dn, 1e-6))
        results[state] = (ne, no, dn, kappa, lc)
    
    dn_am, kappa_am, lc_am = results["amorphous"][2], results["amorphous"][3], results["amorphous"][4]
    dn_cr, kappa_cr, lc_cr = results["crystalline"][2], results["crystalline"][3], results["crystalline"][4]
    ratio = kappa_cr / max(kappa_am, 1e-6)
    return lc_am, lc_cr, ratio

if __name__ == "__main__":
    print("=" * 70)
    print("SWEEPING GAP FOR SLOT COUPLER (Sb2S3 IN GAP, w_si=240nm, h=220nm):")
    print(f"{'Gap (nm)':<10} {'L_c (Am) [um]':<16} {'L_c (Cr) [um]':<16} {'kappa_cr/kappa_am':<18}")
    print("-" * 70)
    for g in [0.040, 0.050, 0.060, 0.070, 0.080, 0.090, 0.100, 0.120]:
        lc_am, lc_cr, ratio = solve_slot_coupler(w_si=0.240, gap=g)
        print(f"{g*1000:<10.0f} {lc_am:<16.3f} {lc_cr:<16.3f} {ratio:<18.3f}")
    print("=" * 70)
