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

def solve_slot_coupler(w=0.240, gap=0.030, h=0.220, state="amorphous"):
    lambda_0 = 1.064
    fcen = 1.0 / lambda_0
    sy, sz = 3.0, 2.0
    cell = mp.Vector3(0, sy, sz)
    clad = mp.Medium(index=1.449)
    core = mp.Medium(index=3.480)
    n_pcm = 2.70 if state == "amorphous" else 3.30

    y_sep = w + gap
    y1 = y_sep / 2.0
    y2 = -y_sep / 2.0

    geom = [
        mp.Block(center=mp.Vector3(0, 0, -sz/4.0), size=mp.Vector3(mp.inf, mp.inf, sz/2.0), material=clad),
        # Silicon cores
        mp.Block(center=mp.Vector3(0, y1, h/2.0), size=mp.Vector3(mp.inf, w, h), material=core),
        mp.Block(center=mp.Vector3(0, y2, h/2.0), size=mp.Vector3(mp.inf, w, h), material=core),
        # Slot PCM in the gap between the two waveguides
        mp.Block(center=mp.Vector3(0, 0, h/2.0), size=mp.Vector3(mp.inf, gap, h), material=mp.Medium(index=n_pcm)),
    ]

    k_guess = 2.45 * fcen
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
    print("=" * 75)
    print("SLOT-COUPLED DIRECTIONAL SWITCH (PCM IN THE GAP)")
    print("=" * 75)
    print(f"{'w (nm)':<10} {'gap (nm)':<10} {'L_c(am) um':<12} {'L_c(cr) um':<12} {'Ratio':<10} {'Area um^2':<10}")
    print("-" * 75)
    for w in [0.200, 0.220, 0.240]:
        for g in [0.020, 0.025, 0.030, 0.035, 0.040]:
            with nostdout():
                k1_am, k2_am = solve_slot_coupler(w=w, gap=g, state="amorphous")
                k1_cr, k2_cr = solve_slot_coupler(w=w, gap=g, state="crystalline")
            dn_am = k1_am - k2_am
            dn_cr = k1_cr - k2_cr
            lc_am = 1.064 / (2.0 * dn_am)
            lc_cr = 1.064 / (2.0 * dn_cr)
            ratio = lc_am / lc_cr
            area = (lc_am + 0.95) * 0.61
            print(f"{w*1000:<10.0f} {g*1000:<10.0f} {lc_am:<12.2f} {lc_cr:<12.2f} {ratio:<10.3f} {area:<10.2f}")
