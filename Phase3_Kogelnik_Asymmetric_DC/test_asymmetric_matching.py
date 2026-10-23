import math
import meep as mp
from meep import mpb

fcen = 1.0 / 1.064
clad = mp.Medium(index=1.449)
core = mp.Medium(index=3.480)
h_wg = 0.220
t_pcm = 0.190
w1 = 0.240
w2 = 0.28687
gap = 0.09814
cell = mp.Vector3(0, 4.0, 2.5)

d = (w1 + w2)/2.0 + gap
y1 = +d/2.0
y2 = -d/2.0

# 1. State 0: WG1 with a-Sb2S3, WG2 bare Si
geom_s0 = [
    mp.Block(center=mp.Vector3(0, 0, -0.6), size=mp.Vector3(mp.inf, mp.inf, 1.2), material=clad),
    mp.Block(center=mp.Vector3(0, y1, h_wg/2.0), size=mp.Vector3(mp.inf, w1, h_wg), material=core),
    mp.Block(center=mp.Vector3(0, y1, h_wg + t_pcm/2.0), size=mp.Vector3(mp.inf, w1, t_pcm), material=mp.Medium(index=2.70)),
    mp.Block(center=mp.Vector3(0, y2, h_wg/2.0), size=mp.Vector3(mp.inf, w2, h_wg), material=core),
]
ms0 = mpb.ModeSolver(geometry=geom_s0, geometry_lattice=mp.Lattice(size=cell),
                     k_points=[mp.Vector3(2.2*fcen, 0, 0)], resolution=32,
                     default_material=clad, num_bands=2)
k_even_0 = ms0.find_k(mp.NO_PARITY, fcen, 1, 1, mp.Vector3(1, 0, 0), 1e-4, 2.0*fcen, 0.5*fcen, 3.48*fcen)[0]
k_odd_0  = ms0.find_k(mp.NO_PARITY, fcen, 2, 2, mp.Vector3(1, 0, 0), 1e-4, 2.0*fcen, 0.5*fcen, 3.48*fcen)[0]
lc_0 = 1.0 / (2.0 * (k_even_0 - k_odd_0))

# 2. State 1: WG1 with c-Sb2S3, WG2 bare Si
geom_s1 = [
    mp.Block(center=mp.Vector3(0, 0, -0.6), size=mp.Vector3(mp.inf, mp.inf, 1.2), material=clad),
    mp.Block(center=mp.Vector3(0, y1, h_wg/2.0), size=mp.Vector3(mp.inf, w1, h_wg), material=core),
    mp.Block(center=mp.Vector3(0, y1, h_wg + t_pcm/2.0), size=mp.Vector3(mp.inf, w1, t_pcm), material=mp.Medium(index=3.30)),
    mp.Block(center=mp.Vector3(0, y2, h_wg/2.0), size=mp.Vector3(mp.inf, w2, h_wg), material=core),
]
ms1 = mpb.ModeSolver(geometry=geom_s1, geometry_lattice=mp.Lattice(size=cell),
                     k_points=[mp.Vector3(2.2*fcen, 0, 0)], resolution=32,
                     default_material=clad, num_bands=2)
k_even_1 = ms1.find_k(mp.NO_PARITY, fcen, 1, 1, mp.Vector3(1, 0, 0), 1e-4, 2.0*fcen, 0.5*fcen, 3.48*fcen)[0]
k_odd_1  = ms1.find_k(mp.NO_PARITY, fcen, 2, 2, mp.Vector3(1, 0, 0), 1e-4, 2.0*fcen, 0.5*fcen, 3.48*fcen)[0]
lc_1 = 1.0 / (k_even_1 - k_odd_1)

delta_n_0 = (k_even_0 - k_odd_0) / fcen
delta_n_1 = (k_even_1 - k_odd_1) / fcen

print("=" * 80)
print(f"ASYMMETRIC CORE MATCHING RESULTS (w1 = {w1*1e3:.1f} nm, w2 = {w2*1e3:.2f} nm, Gap = {gap*1e3:.2f} nm)")
print("=" * 80)
print(f"  State 0 Beat Length L_c(0) = {lc_0:.10f} um  (Delta_n = {delta_n_0:.10f})")
print(f"  State 1 Null Period L_c(1) = {lc_1:.10f} um  (Delta_n = {delta_n_1:.10f})")
print(f"  Offset L_c(0) - L_c(1)     = {(lc_0 - lc_1)*1e3:.10f} nm")
print("=" * 80)
