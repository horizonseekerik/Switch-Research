import math, meep as mp
from meep import mpb

fcen = 1.0 / 1.064
clad = mp.Medium(index=1.449)
core = mp.Medium(index=3.480)
h_wg = 0.220
cell = mp.Vector3(0, 4.0, 3.0)

# Target: n_eff of WG1 with Amorphous PCM (w=240nm, h=220nm, tpcm=160nm)
geom_am = [
    mp.Block(center=mp.Vector3(0, 0, -0.75), size=mp.Vector3(mp.inf, mp.inf, 1.5), material=clad),
    mp.Block(center=mp.Vector3(0, 0, h_wg/2.0), size=mp.Vector3(mp.inf, 0.240, h_wg), material=core),
    mp.Block(center=mp.Vector3(0, 0, h_wg + 0.160/2.0), size=mp.Vector3(mp.inf, 0.240, 0.160), material=mp.Medium(index=2.70)),
]
ms = mpb.ModeSolver(geometry=geom_am, geometry_lattice=mp.Lattice(size=cell),
                    k_points=[mp.Vector3(2.5*fcen, 0, 0)], resolution=30,
                    default_material=clad, num_bands=1)
k_am = ms.find_k(mp.NO_PARITY, fcen, 1, 1, mp.Vector3(1, 0, 0), 1e-4, 2.5*fcen, 0.5*fcen, 3.48*fcen)
target_nam = k_am[0] / fcen
print(f"Target WG1 (Amorphous PCM) n_eff = {target_nam:.4f}")
print("-" * 55)
print(f"{'w_bare (nm)':<15} {'n_eff bare Si':<15} {'Diff (n_bare - target)':<20}")
print("-" * 55)

for w_bare in [0.280, 0.290, 0.300, 0.310, 0.320, 0.330, 0.340]:
    geom_bare = [
        mp.Block(center=mp.Vector3(0, 0, -0.75), size=mp.Vector3(mp.inf, mp.inf, 1.5), material=clad),
        mp.Block(center=mp.Vector3(0, 0, h_wg/2.0), size=mp.Vector3(mp.inf, w_bare, h_wg), material=core),
    ]
    ms.geometry = geom_bare
    k_bare = ms.find_k(mp.NO_PARITY, fcen, 1, 1, mp.Vector3(1, 0, 0), 1e-4, 2.5*fcen, 0.5*fcen, 3.48*fcen)
    n_bare = k_bare[0] / fcen
    diff = n_bare - target_nam
    print(f"{w_bare*1000:<15.0f} {n_bare:<15.4f} {diff:<20.4f}")
