"""
MPB SOLVER: SCAN CORE WIDTH (w = 200 - 240 nm) TO HIT u = sqrt(3) ~ 1.73205
===========================================================================
Finds the exact core width where:
  u = Delta_beta / (2 * kappa) == sqrt(3) ~ 1.7320508
  L_c(0) == L_null(1) in 3D FDTD!
"""

import math
import numpy as np
import meep as mp
from meep import mpb

fcen = 1.0 / 1.064
clad = mp.Medium(index=1.449)
core = mp.Medium(index=3.480)
h_wg = 0.220
t_pcm = 0.190
cell = mp.Vector3(0, 3.5, 2.5)

print("=" * 115)
print(f"{'Core w (nm)':<12} {'n_am':<12} {'n_cr':<12} {'Delta_beta (rad/um)':<22} {'kappa (rad/um)':<16} {'u = Dbeta/2k':<14} {'Target sqrt(3)':<14}")
print("-" * 115)

for w_nm in [200.0, 210.0, 220.0, 230.0, 240.0]:
    w = w_nm * 1e-3
    gap = 0.09814

    # 1. Single waveguide amorphous
    geom_am = [
        mp.Block(center=mp.Vector3(0, 0, -0.6), size=mp.Vector3(mp.inf, mp.inf, 1.2), material=clad),
        mp.Block(center=mp.Vector3(0, 0, h_wg/2.0), size=mp.Vector3(mp.inf, w, h_wg), material=core),
        mp.Block(center=mp.Vector3(0, 0, h_wg + t_pcm/2.0), size=mp.Vector3(mp.inf, w, t_pcm), material=mp.Medium(index=2.70)),
    ]
    ms_am = mpb.ModeSolver(geometry=geom_am, geometry_lattice=mp.Lattice(size=cell),
                           k_points=[mp.Vector3(2.2*fcen, 0, 0)], resolution=30,
                           default_material=clad, num_bands=1)
    k_am = ms_am.find_k(mp.NO_PARITY, fcen, 1, 1, mp.Vector3(1, 0, 0), 1e-4, 2.2*fcen, 0.5*fcen, 3.48*fcen)[0]
    n_am = k_am / fcen

    # 2. Single waveguide crystalline
    geom_cr = [
        mp.Block(center=mp.Vector3(0, 0, -0.6), size=mp.Vector3(mp.inf, mp.inf, 1.2), material=clad),
        mp.Block(center=mp.Vector3(0, 0, h_wg/2.0), size=mp.Vector3(mp.inf, w, h_wg), material=core),
        mp.Block(center=mp.Vector3(0, 0, h_wg + t_pcm/2.0), size=mp.Vector3(mp.inf, w, t_pcm), material=mp.Medium(index=3.30)),
    ]
    ms_am.geometry = geom_cr
    k_cr = ms_am.find_k(mp.NO_PARITY, fcen, 1, 1, mp.Vector3(1, 0, 0), 1e-4, 2.2*fcen, 0.5*fcen, 3.48*fcen)[0]
    n_cr = k_cr / fcen

    dbeta = 2.0 * math.pi * fcen * (n_cr - n_am)

    # 3. Directional coupler supermodes in amorphous state
    d = gap + w
    geom_dc_am = [
        mp.Block(center=mp.Vector3(0, 0, -0.6), size=mp.Vector3(mp.inf, mp.inf, 1.2), material=clad),
        mp.Block(center=mp.Vector3(0, +d/2.0, h_wg/2.0), size=mp.Vector3(mp.inf, w, h_wg), material=core),
        mp.Block(center=mp.Vector3(0, +d/2.0, h_wg + t_pcm/2.0), size=mp.Vector3(mp.inf, w, t_pcm), material=mp.Medium(index=2.70)),
        mp.Block(center=mp.Vector3(0, -d/2.0, h_wg/2.0), size=mp.Vector3(mp.inf, w, h_wg), material=core),
        mp.Block(center=mp.Vector3(0, -d/2.0, h_wg + t_pcm/2.0), size=mp.Vector3(mp.inf, w, t_pcm), material=mp.Medium(index=2.70)),
    ]
    ms_dc = mpb.ModeSolver(geometry=geom_dc_am, geometry_lattice=mp.Lattice(size=cell),
                           k_points=[mp.Vector3(2.2*fcen, 0, 0)], resolution=30,
                           default_material=clad, num_bands=2)
    k_even = ms_dc.find_k(mp.NO_PARITY, fcen, 1, 1, mp.Vector3(1, 0, 0), 1e-4, 2.2*fcen, 0.5*fcen, 3.48*fcen)[0]
    k_odd  = ms_dc.find_k(mp.NO_PARITY, fcen, 2, 2, mp.Vector3(1, 0, 0), 1e-4, 2.2*fcen, 0.5*fcen, 3.48*fcen)[0]

    kappa = (k_even - k_odd) * math.pi * fcen
    u = dbeta / (2.0 * kappa) if kappa > 0 else 0.0

    print(f"{w_nm:<12.1f} {n_am:<12.5f} {n_cr:<12.5f} {dbeta:<22.5f} {kappa:<16.5f} {u:<14.5f} {math.sqrt(3):<14.5f}")

print("=" * 115)
