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

def solve_symmetric_pcm_coupler(w=0.240, gap=0.035, h=0.220, t_pcm=0.100, state="amorphous"):
    lambda_0 = 1.064
    fcen = 1.0 / lambda_0
    sy = 3.0
    sz = 2.0
    cell = mp.Vector3(0, sy, sz)
    clad = mp.Medium(index=1.449)
    core = mp.Medium(index=3.480)
    n_pcm = 2.70 if state == "amorphous" else 3.30

    y_sep = w + gap
    y1 = y_sep / 2.0
    y2 = -y_sep / 2.0

    geom = [
        mp.Block(center=mp.Vector3(0, 0, -sz/4.0), size=mp.Vector3(mp.inf, mp.inf, sz/2.0), material=clad),
        # Both WG1 and WG2
        mp.Block(center=mp.Vector3(0, y1, h/2.0), size=mp.Vector3(mp.inf, w, h), material=core),
        mp.Block(center=mp.Vector3(0, y2, h/2.0), size=mp.Vector3(mp.inf, w, h), material=core),
        # PCM on top of BOTH waveguides
        mp.Block(center=mp.Vector3(0, y1, h + t_pcm/2.0), size=mp.Vector3(mp.inf, w, t_pcm), material=mp.Medium(index=n_pcm)),
        mp.Block(center=mp.Vector3(0, y2, h + t_pcm/2.0), size=mp.Vector3(mp.inf, w, t_pcm), material=mp.Medium(index=n_pcm)),
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
    print("=" * 70)
    print("SYMMETRIC DUAL-COATED COUPLER (COUPLING-MODULATED)")
    print("=" * 70)
    for tp in [0.050, 0.080, 0.100, 0.150]:
        with nostdout():
            k1_am, k2_am = solve_symmetric_pcm_coupler(t_pcm=tp, state="amorphous")
            k1_cr, k2_cr = solve_symmetric_pcm_coupler(t_pcm=tp, state="crystalline")
        dn_am = k1_am - k2_am
        dn_cr = k1_cr - k2_cr
        lc_am = 1.064 / (2.0 * dn_am)
        lc_cr = 1.064 / (2.0 * dn_cr)
        print(f"t_pcm = {tp*1000:3.0f} nm:")
        print(f"  Amorphous  : dn = {dn_am:.5f}, L_c = {lc_am:.2f} um")
        print(f"  Crystalline: dn = {dn_cr:.5f}, L_c = {lc_cr:.2f} um")
        print(f"  Ratio L_c(am)/L_c(cr) = {lc_am/lc_cr:.3f}")
