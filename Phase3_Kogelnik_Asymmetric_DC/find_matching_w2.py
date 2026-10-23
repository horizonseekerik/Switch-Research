import math
import meep as mp
from meep import mpb

fcen = 1.0 / 1.064
clad = mp.Medium(index=1.449)
core = mp.Medium(index=3.480)
h_wg = 0.220
t_pcm = 0.190
cell = mp.Vector3(0, 3.0, 2.5)

# Target: Mode index of WG1 (w=240 nm, h=220 nm, t_pcm=190 nm a-Sb2S3 n=2.70)
geom_wg1_am = [
    mp.Block(center=mp.Vector3(0, 0, -0.6), size=mp.Vector3(mp.inf, mp.inf, 1.2), material=clad),
    mp.Block(center=mp.Vector3(0, 0, h_wg/2.0), size=mp.Vector3(mp.inf, 0.240, h_wg), material=core),
    mp.Block(center=mp.Vector3(0, 0, h_wg + t_pcm/2.0), size=mp.Vector3(mp.inf, 0.240, t_pcm), material=mp.Medium(index=2.70)),
]
ms1 = mpb.ModeSolver(geometry=geom_wg1_am, geometry_lattice=mp.Lattice(size=cell),
                     k_points=[mp.Vector3(2.2*fcen, 0, 0)], resolution=32,
                     default_material=clad, num_bands=1)
k_am = ms1.find_k(mp.NO_PARITY, fcen, 1, 1, mp.Vector3(1, 0, 0), 1e-4, 2.0*fcen, 0.5*fcen, 3.48*fcen)[0]
n_eff_target = k_am / fcen

print("=" * 80)
print(f"TARGET MODE INDEX (WG1 with a-Sb2S3): n_eff = {n_eff_target:.10f}")
print("=" * 80)
print(f"{'Bare Si Width w2 (nm)':<25} {'n_eff (bare Si)':<25} {'Delta_n':<20}")
print("-" * 80)

for w2_nm in [250.0, 260.0, 270.0, 280.0, 290.0, 300.0, 310.0, 320.0]:
    w2 = w2_nm * 1e-3
    geom_bare = [
        mp.Block(center=mp.Vector3(0, 0, -0.6), size=mp.Vector3(mp.inf, mp.inf, 1.2), material=clad),
        mp.Block(center=mp.Vector3(0, 0, h_wg/2.0), size=mp.Vector3(mp.inf, w2, h_wg), material=core),
    ]
    ms1.geometry = geom_bare
    k_bare = ms1.find_k(mp.NO_PARITY, fcen, 1, 1, mp.Vector3(1, 0, 0), 1e-4, 2.0*fcen, 0.5*fcen, 3.48*fcen)[0]
    n_bare = k_bare / fcen
    dn = n_bare - n_eff_target
    print(f"{w2_nm:<25.2f} {n_bare:<25.10f} {dn:<20.10f}")
print("=" * 80)
