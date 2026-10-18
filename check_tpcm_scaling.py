import math, meep as mp
from meep import mpb

fcen = 1.0 / 1.064
clad = mp.Medium(index=1.449)
core = mp.Medium(index=3.480)
w_wg, h_wg = 0.240, 0.220
cell = mp.Vector3(0, 4.0, 3.0)

print(f"{'t_pcm (nm)':<10} | {'n_am':<7} | {'n_cr':<7} | {'dn':<7} | {'dbeta (rad/um)':<15} | {'target_kappa':<12}")
print('-' * 75)

for tp in [0.140, 0.160, 0.180, 0.200, 0.220]:
    geom_am = [
        mp.Block(center=mp.Vector3(0, 0, -0.75), size=mp.Vector3(mp.inf, mp.inf, 1.5), material=clad),
        mp.Block(center=mp.Vector3(0, 0, h_wg/2.0), size=mp.Vector3(mp.inf, w_wg, h_wg), material=core),
        mp.Block(center=mp.Vector3(0, 0, h_wg + tp/2.0), size=mp.Vector3(mp.inf, w_wg, tp), material=mp.Medium(index=2.70)),
    ]
    ms = mpb.ModeSolver(geometry=geom_am, geometry_lattice=mp.Lattice(size=cell),
                        k_points=[mp.Vector3(2.5*fcen, 0, 0)], resolution=30,
                        default_material=clad, num_bands=1)
    k_am = ms.find_k(mp.NO_PARITY, fcen, 1, 1, mp.Vector3(1, 0, 0), 1e-3, 2.5*fcen, 0.5*fcen, 3.48*fcen)
    
    geom_cr = [
        mp.Block(center=mp.Vector3(0, 0, -0.75), size=mp.Vector3(mp.inf, mp.inf, 1.5), material=clad),
        mp.Block(center=mp.Vector3(0, 0, h_wg/2.0), size=mp.Vector3(mp.inf, w_wg, h_wg), material=core),
        mp.Block(center=mp.Vector3(0, 0, h_wg + tp/2.0), size=mp.Vector3(mp.inf, w_wg, tp), material=mp.Medium(index=3.30)),
    ]
    ms.geometry = geom_cr
    k_cr = ms.find_k(mp.NO_PARITY, fcen, 1, 1, mp.Vector3(1, 0, 0), 1e-3, 2.5*fcen, 0.5*fcen, 3.48*fcen)
    
    nam = k_am[0] / fcen
    ncr = k_cr[0] / fcen
    dn = ncr - nam
    dbeta = (2 * math.pi / 1.064) * dn
    target_k = dbeta / (2 * math.sqrt(3.0))
    print(f"{tp*1000:<10.0f} | {nam:<7.4f} | {ncr:<7.4f} | {dn:<7.4f} | {dbeta:<15.4f} | {target_k:<12.4f}")
