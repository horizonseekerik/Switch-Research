"""
EXACT SCAN AT w = 223.8 nm
"""
import math
import meep as mp
from meep import mpb

fcen = 1.0 / 1.064
clad = mp.Medium(index=1.449)
core = mp.Medium(index=3.480)
h_wg = 0.220
t_pcm = 0.190
cell = mp.Vector3(0, 3.5, 2.5)

gap = 0.09814

for w_nm in [223.0, 223.5, 223.8, 224.0, 225.0]:
    w = w_nm * 1e-3
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

    geom_cr = [
        mp.Block(center=mp.Vector3(0, 0, -0.6), size=mp.Vector3(mp.inf, mp.inf, 1.2), material=clad),
        mp.Block(center=mp.Vector3(0, 0, h_wg/2.0), size=mp.Vector3(mp.inf, w, h_wg), material=core),
        mp.Block(center=mp.Vector3(0, 0, h_wg + t_pcm/2.0), size=mp.Vector3(mp.inf, w, t_pcm), material=mp.Medium(index=3.30)),
    ]
    ms_am.geometry = geom_cr
    k_cr = ms_am.find_k(mp.NO_PARITY, fcen, 1, 1, mp.Vector3(1, 0, 0), 1e-4, 2.2*fcen, 0.5*fcen, 3.48*fcen)[0]
    n_cr = k_cr / fcen
    dbeta = 2.0 * math.pi * fcen * (n_cr - n_am)

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
    u = dbeta / (2.0 * kappa)
    lc_0 = math.pi / (2.0 * kappa)
    S = math.sqrt(kappa**2 + (dbeta/2.0)**2)
    lc_1 = math.pi / S

    print(f"w = {w_nm:.1f} nm: dbeta = {dbeta:.5f}, kappa = {kappa:.5f}, u = {u:.5f} (target {math.sqrt(3):.5f}), Lc(0) = {lc_0:.4f} um, Lc(1) = {lc_1:.4f} um, diff = {abs(lc_0 - lc_1)*1e3:.1f} nm")
