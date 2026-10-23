"""
PROJECT JANUS: SYMMETRIC COUPLING MODULATION SOLVER
Calculate kappa_cr / kappa_am when PCM is on TOP of both waveguides
"""
import sys, os, math
import numpy as np
import meep as mp
from meep import mpb

def solve_cm_modes(w=0.240, h=0.220, gap=0.040, t_pcm=0.150, resolution=40):
    lambda_0 = 1.064
    fcen = 1.0 / lambda_0
    clad_mat = mp.Medium(index=1.449)
    core_mat = mp.Medium(index=3.480)

    y_sep = w + gap
    y1 = y_sep / 2.0
    y2 = -y_sep / 2.0

    sy, sz = 3.0, 2.5
    cell = mp.Vector3(0, sy, sz)

    results = {}
    for state, n_pcm in [("amorphous", 2.70), ("crystalline", 3.30)]:
        pcm_mat = mp.Medium(index=n_pcm)
        geom = [
            mp.Block(center=mp.Vector3(0, 0, -sz/4.0), size=mp.Vector3(mp.inf, mp.inf, sz/2.0), material=clad_mat),
            mp.Block(center=mp.Vector3(0, y1, h/2.0), size=mp.Vector3(mp.inf, w, h), material=core_mat),
            mp.Block(center=mp.Vector3(0, y2, h/2.0), size=mp.Vector3(mp.inf, w, h), material=core_mat),
            mp.Block(center=mp.Vector3(0, y1, h + t_pcm/2.0), size=mp.Vector3(mp.inf, w, t_pcm), material=pcm_mat),
            mp.Block(center=mp.Vector3(0, y2, h + t_pcm/2.0), size=mp.Vector3(mp.inf, w, t_pcm), material=pcm_mat),
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

    lc_am = results["amorphous"][4]
    lc_cr = results["crystalline"][4]
    ratio = lc_am / lc_cr  # kappa_cr / kappa_am
    return lc_am, lc_cr, ratio

if __name__ == "__main__":
    print("=" * 70)
    print("SWEEPING COUPLING MODULATION RATIO (PCM ON TOP OF BOTH WAVEGUIDES):")
    print(f"{'w (nm)':<8} {'Gap (nm)':<10} {'t_pcm (nm)':<12} {'L_c(Am)':<12} {'L_c(Cr)':<12} {'Ratio (k_cr/k_am)':<18}")
    print("-" * 70)
    for w in [0.200, 0.220, 0.240]:
        for gap in [0.030, 0.040, 0.050]:
            for t_pcm in [0.080, 0.120, 0.150]:
                lc_am, lc_cr, ratio = solve_cm_modes(w=w, gap=gap, t_pcm=t_pcm)
                print(f"{w*1000:<8.0f} {gap*1000:<10.0f} {t_pcm*1000:<12.0f} {lc_am:<12.3f} {lc_cr:<12.3f} {ratio:<18.3f}")
    print("=" * 70)
