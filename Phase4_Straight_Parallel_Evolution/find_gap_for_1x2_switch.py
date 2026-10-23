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

def solve_asym_supermodes(w1, w2, gap, t_pcm=0.150):
    lambda_0 = 1.064
    fcen = 1.0 / lambda_0
    sy, sz = 3.0, 2.0
    cell = mp.Vector3(0, sy, sz)
    clad = mp.Medium(index=1.449)
    core = mp.Medium(index=3.480)

    y_sep = (w1 + w2) / 2.0 + gap
    y1 = y_sep / 2.0
    y2 = -y_sep / 2.0

    geom = [
        mp.Block(center=mp.Vector3(0, 0, -sz/4.0), size=mp.Vector3(mp.inf, mp.inf, sz/2.0), material=clad),
        # WG1 (with PCM)
        mp.Block(center=mp.Vector3(0, y1, 0.220/2.0), size=mp.Vector3(mp.inf, w1, 0.220), material=core),
        # WG2 (bare Si)
        mp.Block(center=mp.Vector3(0, y2, 0.220/2.0), size=mp.Vector3(mp.inf, w2, 0.220), material=core),
        # PCM on top of WG1 (Amorphous state)
        mp.Block(center=mp.Vector3(0, y1, 0.220 + t_pcm/2.0), size=mp.Vector3(mp.inf, w1, t_pcm), material=mp.Medium(index=2.70))
    ]

    k_guess = 2.35 * fcen
    ms = mpb.ModeSolver(
        geometry=geom,
        geometry_lattice=mp.Lattice(size=cell),
        k_points=[mp.Vector3(k_guess, 0, 0)],
        resolution=32,
        default_material=clad,
        num_bands=4
    )
    k1 = ms.find_k(mp.NO_PARITY, fcen, 1, 1, mp.Vector3(1, 0, 0), 1e-4, k_guess, 0.5*fcen, 3.48*fcen)[0] / fcen
    k2 = ms.find_k(mp.NO_PARITY, fcen, 2, 2, mp.Vector3(1, 0, 0), 1e-4, k_guess, 0.5*fcen, 3.48*fcen)[0] / fcen
    return k1, k2

print("=" * 65)
print("FINDING EXACT GAP FOR L_c = 4.60 um (w1 = 185 nm, w2 = 250 nm)")
print("=" * 65)
print(f"{'gap (nm)':<10} {'n_even':<10} {'n_odd':<10} {'Delta_n':<10} {'kappa (um^-1)':<15} {'L_c (um)':<10}")
print("-" * 65)
for g in [0.020, 0.025, 0.030, 0.035, 0.040, 0.045, 0.050]:
    with nostdout():
        k1, k2 = solve_asym_supermodes(w1=0.185, w2=0.250, gap=g)
    dn = k1 - k2
    kappa = (math.pi / 1.064) * dn
    lc = 1.064 / (2.0 * dn)
    print(f"{g*1000:<10.0f} {k1:<10.4f} {k2:<10.4f} {dn:<10.5f} {kappa:<15.4f} {lc:<10.2f}")
