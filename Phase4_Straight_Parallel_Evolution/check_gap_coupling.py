import sys, os, math
import numpy as np
import meep as mp
from meep import mpb

fcen = 1.0 / 1.064
w_wg, h_wg = 0.240, 0.220
clad = mp.Medium(index=1.449)
core = mp.Medium(index=3.480)
cell = mp.Vector3(0, 5.0, 3.0)

print(f"{'Gap (nm)':<10} {'ne':<10} {'no':<10} {'dn':<10} {'kappa (rad/um)':<15} {'L_pi/2 (um)':<15}")
print("-" * 75)

for gap in [0.090, 0.150, 0.210, 0.300, 0.400, 0.500, 0.600, 0.800]:
    y_sep = w_wg + gap
    y1, y2 = y_sep / 2.0, -y_sep / 2.0
    geom = [
        mp.Block(center=mp.Vector3(0, 0, -1.5 / 2.0), size=mp.Vector3(mp.inf, mp.inf, 1.5), material=clad),
        mp.Block(center=mp.Vector3(0, y1, h_wg / 2.0), size=mp.Vector3(mp.inf, w_wg, h_wg), material=core),
        mp.Block(center=mp.Vector3(0, y2, h_wg / 2.0), size=mp.Vector3(mp.inf, w_wg, h_wg), material=core),
    ]
    ms = mpb.ModeSolver(geometry=geom, geometry_lattice=mp.Lattice(size=cell),
                        k_points=[mp.Vector3(2.5 * fcen, 0, 0)], resolution=30,
                        default_material=clad, num_bands=4)
    old_stdout = sys.stdout
    sys.stdout = open(os.devnull, 'w')
    k_res = ms.find_k(mp.NO_PARITY, fcen, 1, 2, mp.Vector3(1, 0, 0), 1e-4, 2.5 * fcen, 0.5 * fcen, 3.48 * fcen)
    sys.stdout.close()
    sys.stdout = old_stdout

    ne = max(k_res[0], k_res[1]) / fcen
    no = min(k_res[0], k_res[1]) / fcen
    dn = ne - no
    kappa = (math.pi / 1.064) * dn
    lc = math.pi / (2.0 * kappa) if kappa > 0 else 1e9

    print(f"{gap*1000:<10.0f} {ne:<10.4f} {no:<10.4f} {dn:<10.5f} {kappa:<15.4f} {lc:<15.2f}")
