import math
import meep as mp
from meep import mpb

fcen = 1.0 / 1.064
clad = mp.Medium(index=1.449)
core = mp.Medium(index=3.480)
h_wg = 0.220
w = 0.240
cell = mp.Vector3(0, 3.0, 2.5)

print(f"{'t_pcm (nm)':<12} {'n_eff_am':<14} {'n_eff_cr':<14} {'Delta_n_eff':<14} {'Delta_beta (rad/um)':<20}")
print("-" * 74)

for t_nm in [40, 60, 80, 100, 120, 140, 160, 180, 200, 237.5]:
    t_pcm = t_nm * 1e-3
    # Amorphous
    geom_am = [
        mp.Block(center=mp.Vector3(0, 0, -0.6), size=mp.Vector3(mp.inf, mp.inf, 1.2), material=clad),
        mp.Block(center=mp.Vector3(0, 0, h_wg/2.0), size=mp.Vector3(mp.inf, w, h_wg), material=core),
        mp.Block(center=mp.Vector3(0, 0, h_wg + t_pcm/2.0), size=mp.Vector3(mp.inf, w, t_pcm), material=mp.Medium(index=2.70))
    ]
    ms = mpb.ModeSolver(geometry=geom_am, geometry_lattice=mp.Lattice(size=cell),
                        k_points=[mp.Vector3(2.2*fcen, 0, 0)], resolution=32,
                        default_material=clad, num_bands=1, verbose=False)
    k_am = ms.find_k(mp.NO_PARITY, fcen, 1, 1, mp.Vector3(1, 0, 0), 1e-4, 2.0*fcen, 0.5*fcen, 3.48*fcen)[0]
    n_am = k_am / fcen

    # Crystalline
    geom_cr = [
        mp.Block(center=mp.Vector3(0, 0, -0.6), size=mp.Vector3(mp.inf, mp.inf, 1.2), material=clad),
        mp.Block(center=mp.Vector3(0, 0, h_wg/2.0), size=mp.Vector3(mp.inf, w, h_wg), material=core),
        mp.Block(center=mp.Vector3(0, 0, h_wg + t_pcm/2.0), size=mp.Vector3(mp.inf, w, t_pcm), material=mp.Medium(index=3.30))
    ]
    ms.geometry = geom_cr
    k_cr = ms.find_k(mp.NO_PARITY, fcen, 1, 1, mp.Vector3(1, 0, 0), 1e-4, 2.0*fcen, 0.5*fcen, 3.48*fcen)[0]
    n_cr = k_cr / fcen

    dn = n_cr - n_am
    dbeta = 2 * math.pi * dn / 1.064
    print(f"{t_nm:<12.1f} {n_am:<14.6f} {n_cr:<14.6f} {dn:<14.6f} {dbeta:<20.6f}")
