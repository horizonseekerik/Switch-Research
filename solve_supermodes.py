import math
import numpy as np
import meep as mp
from meep import mpb
import contextlib, os, sys

@contextlib.contextmanager
def nostdout():
    with open(os.devnull, "w") as devnull:
        old = sys.stdout; sys.stdout = devnull
        try: yield
        finally: sys.stdout = old

def solve_coupler_supermodes(w1, w2, gap, t_pcm, pcm_state="crystalline", resolution=32):
    lambda_0 = 1.064
    fcen = 1.0 / lambda_0
    h = 0.220
    sy = 3.0
    sz = 2.0
    cell = mp.Vector3(0, sy, sz)
    clad = mp.Medium(index=1.449)
    core = mp.Medium(index=3.480)

    y_sep = (w1 + w2) / 2.0 + gap
    y1 = y_sep / 2.0
    y2 = -y_sep / 2.0

    geom = [
        mp.Block(center=mp.Vector3(0, 0, -sz/4.0), size=mp.Vector3(mp.inf, mp.inf, sz/2.0), material=clad),
        # WG1 (with PCM)
        mp.Block(center=mp.Vector3(0, y1, h/2.0), size=mp.Vector3(mp.inf, w1, h), material=core),
        # WG2 (bare Si)
        mp.Block(center=mp.Vector3(0, y2, h/2.0), size=mp.Vector3(mp.inf, w2, h), material=core),
    ]
    if t_pcm > 0:
        n_pcm = 2.70 if pcm_state == "amorphous" else 3.30
        geom.append(mp.Block(
            center=mp.Vector3(0, y1, h + t_pcm/2.0),
            size=mp.Vector3(mp.inf, w1, t_pcm),
            material=mp.Medium(index=n_pcm)
        ))

    k_guess = 2.35 * fcen
    ms = mpb.ModeSolver(
        geometry=geom,
        geometry_lattice=mp.Lattice(size=cell),
        k_points=[mp.Vector3(k_guess, 0, 0)],
        resolution=resolution,
        default_material=clad,
        num_bands=4
    )
    # Find mode 1 (even supermode) and mode 2 (odd supermode)
    k1 = ms.find_k(mp.NO_PARITY, fcen, 1, 1, mp.Vector3(1, 0, 0), 1e-4, k_guess, 0.5*fcen, 3.48*fcen)[0] / fcen
    k2 = ms.find_k(mp.NO_PARITY, fcen, 2, 2, mp.Vector3(1, 0, 0), 1e-4, k_guess, 0.5*fcen, 3.48*fcen)[0] / fcen
    return k1, k2

if __name__ == "__main__":
    print("=" * 75)
    print("COUPLED SUPERMODE ANALYSIS FOR RESONANT PHASE-MATCHED SWITCH")
    print("=" * 75)
    # For tp = 100 nm: w1 = 148 nm, w2 = 240 nm
    w1, w2, tp = 0.148, 0.240, 0.100
    print(f"\nConfiguration: w1={w1*1000}nm, w2={w2*1000}nm, t_pcm={tp*1000}nm")
    print(f"{'Gap (nm)':<10} {'n_even':<10} {'n_odd':<10} {'Delta_n':<10} {'kappa (um^-1)':<15} {'L_c (um)':<10}")
    print("-" * 75)
    for gap in [0.030, 0.040, 0.050, 0.060, 0.070, 0.080, 0.100]:
        with nostdout():
            ne, no = solve_coupler_supermodes(w1, w2, gap, tp, pcm_state="crystalline")
        dn = ne - no
        kappa = (math.pi / 1.064) * dn
        lc = 1.064 / (2.0 * dn)
        print(f"{gap*1000:<10.0f} {ne:<10.4f} {no:<10.4f} {dn:<10.5f} {kappa:<15.4f} {lc:<10.2f}")
