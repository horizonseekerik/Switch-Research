import sys, os, math
import meep as mp
from meep import mpb

fcen = 1.0 / 1.064
h_wg = 0.220
clad = mp.Medium(index=1.449)
core = mp.Medium(index=3.480)
cell = mp.Vector3(0, 4.0, 3.0)

print(f"{'Width (nm)':<12} {'n_eff':<10} {'dn vs 240nm':<15} {'delta_beta (rad/um)':<22} {'u = dbeta/(2*0.2)':<20} {'F_max = 1/(1+u^2)':<18}")
print("-" * 100)

geom_ref = [
    mp.Block(center=mp.Vector3(0, 0, -1.5 / 2.0), size=mp.Vector3(mp.inf, mp.inf, 1.5), material=clad),
    mp.Block(center=mp.Vector3(0, 0, h_wg / 2.0), size=mp.Vector3(mp.inf, 0.240, h_wg), material=core),
]
ms_ref = mpb.ModeSolver(geometry=geom_ref, geometry_lattice=mp.Lattice(size=cell),
                        k_points=[mp.Vector3(2.5 * fcen, 0, 0)], resolution=30,
                        default_material=clad, num_bands=2)
old_stdout = sys.stdout
sys.stdout = open(os.devnull, 'w')
k_ref = ms_ref.find_k(mp.NO_PARITY, fcen, 1, 1, mp.Vector3(1, 0, 0), 1e-4, 2.5 * fcen, 0.5 * fcen, 3.48 * fcen)
sys.stdout.close()
sys.stdout = old_stdout
n_ref = k_ref[0] / fcen

for w in [0.240, 0.300, 0.350, 0.400, 0.450, 0.500]:
    geom = [
        mp.Block(center=mp.Vector3(0, 0, -1.5 / 2.0), size=mp.Vector3(mp.inf, mp.inf, 1.5), material=clad),
        mp.Block(center=mp.Vector3(0, 0, h_wg / 2.0), size=mp.Vector3(mp.inf, w, h_wg), material=core),
    ]
    ms = mpb.ModeSolver(geometry=geom, geometry_lattice=mp.Lattice(size=cell),
                        k_points=[mp.Vector3(2.5 * fcen, 0, 0)], resolution=30,
                        default_material=clad, num_bands=2)
    old_stdout = sys.stdout
    sys.stdout = open(os.devnull, 'w')
    k = ms.find_k(mp.NO_PARITY, fcen, 1, 1, mp.Vector3(1, 0, 0), 1e-4, 2.5 * fcen, 0.5 * fcen, 3.48 * fcen)
    sys.stdout.close()
    sys.stdout = old_stdout
    n = k[0] / fcen
    dn = n - n_ref
    dbeta = (2.0 * math.pi / 1.064) * dn
    u = dbeta / (2.0 * 0.2043)
    f_max = 1.0 / (1.0 + u**2)
    print(f"{w*1000:<12.0f} {n:<10.4f} {dn:<15.4f} {dbeta:<22.4f} {u:<20.2f} {f_max*100:<18.2f}%")
