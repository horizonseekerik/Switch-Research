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

def solve_3db_coupler(w=0.250, gap=0.045, h=0.220):
    lambda_0 = 1.064
    fcen = 1.0 / lambda_0
    sy, sz = 3.0, 2.0
    cell = mp.Vector3(0, sy, sz)
    clad = mp.Medium(index=1.449)
    core = mp.Medium(index=3.480)

    y_sep = w + gap
    y1 = y_sep / 2.0
    y2 = -y_sep / 2.0

    geom = [
        mp.Block(center=mp.Vector3(0, 0, -sz/4.0), size=mp.Vector3(mp.inf, mp.inf, sz/2.0), material=clad),
        mp.Block(center=mp.Vector3(0, y1, h/2.0), size=mp.Vector3(mp.inf, w, h), material=core),
        mp.Block(center=mp.Vector3(0, y2, h/2.0), size=mp.Vector3(mp.inf, w, h), material=core),
    ]

    k_guess = 2.40 * fcen
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

if __name__ == "__main__":
    print("=" * 60)
    print("CALCULATING EXACT 3dB COUPLER DIMENSIONS (MPB)")
    print("=" * 60)
    for g in [0.035, 0.040, 0.045, 0.050, 0.055]:
        with nostdout():
            k1, k2 = solve_3db_coupler(w=0.250, gap=g)
        dn = k1 - k2
        lc = 1.064 / (2.0 * dn)
        l_3db = lc / 2.0
        print(f"gap = {g*1000:.0f} nm: dn = {dn:.5f}, L_c = {lc:.3f} um, L_3dB = {l_3db:.3f} um")
