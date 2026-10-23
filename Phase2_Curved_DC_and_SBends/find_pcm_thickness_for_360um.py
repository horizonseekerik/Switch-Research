import math
import meep as mp
from meep import mpb

fcen = 1.0 / 1.064
clad = mp.Medium(index=1.449)
core = mp.Medium(index=3.480)
h_wg = 0.220
w = 0.240
gap = 0.09814
d = gap + w
cell = mp.Vector3(0, 4.0, 2.5)

print("=" * 85)
print(f"PCM THICKNESS SCAN AT GAP = {gap*1e3:.2f} nm, w = {w*1e3:.1f} nm")
print("=" * 85)
print(f"{'t_pcm (nm)':<12} {'Delta_n':<14} {'Lc_0 (um)':<14} {'Lc_1 (um)':<14} {'Diff (nm)':<14}")
print("-" * 85)

for t_nm in [190.0, 210.0, 230.0, 250.0, 270.0, 300.0]:
    t_pcm = t_nm * 1e-3

    # State 0 supermodes (both amorphous)
    geom_0 = [
        mp.Block(center=mp.Vector3(0, 0, -0.6), size=mp.Vector3(mp.inf, mp.inf, 1.2), material=clad),
        mp.Block(center=mp.Vector3(0, +d/2.0, h_wg/2.0), size=mp.Vector3(mp.inf, w, h_wg), material=core),
        mp.Block(center=mp.Vector3(0, +d/2.0, h_wg + t_pcm/2.0), size=mp.Vector3(mp.inf, w, t_pcm), material=mp.Medium(index=2.70)),
        mp.Block(center=mp.Vector3(0, -d/2.0, h_wg/2.0), size=mp.Vector3(mp.inf, w, h_wg), material=core),
        mp.Block(center=mp.Vector3(0, -d/2.0, h_wg + t_pcm/2.0), size=mp.Vector3(mp.inf, w, t_pcm), material=mp.Medium(index=2.70)),
    ]
    ms = mpb.ModeSolver(geometry=geom_0, geometry_lattice=mp.Lattice(size=cell),
                        k_points=[mp.Vector3(2.2*fcen, 0, 0)], resolution=32,
                        default_material=clad, num_bands=2, verbose=False)
    k_e0 = ms.find_k(mp.NO_PARITY, fcen, 1, 1, mp.Vector3(1, 0, 0), 1e-4, 2.0*fcen, 0.5*fcen, 3.48*fcen)[0]
    k_o0 = ms.find_k(mp.NO_PARITY, fcen, 2, 2, mp.Vector3(1, 0, 0), 1e-4, 2.0*fcen, 0.5*fcen, 3.48*fcen)[0]
    lc_0 = 1.0 / (2.0 * (k_e0 - k_o0))

    # State 1 supermodes (WG1 crystalline, WG2 amorphous)
    geom_1 = [
        mp.Block(center=mp.Vector3(0, 0, -0.6), size=mp.Vector3(mp.inf, mp.inf, 1.2), material=clad),
        mp.Block(center=mp.Vector3(0, +d/2.0, h_wg/2.0), size=mp.Vector3(mp.inf, w, h_wg), material=core),
        mp.Block(center=mp.Vector3(0, +d/2.0, h_wg + t_pcm/2.0), size=mp.Vector3(mp.inf, w, t_pcm), material=mp.Medium(index=3.30)),
        mp.Block(center=mp.Vector3(0, -d/2.0, h_wg/2.0), size=mp.Vector3(mp.inf, w, h_wg), material=core),
        mp.Block(center=mp.Vector3(0, -d/2.0, h_wg + t_pcm/2.0), size=mp.Vector3(mp.inf, w, t_pcm), material=mp.Medium(index=2.70)),
    ]
    ms.geometry = geom_1
    k_e1 = ms.find_k(mp.NO_PARITY, fcen, 1, 1, mp.Vector3(1, 0, 0), 1e-4, 2.0*fcen, 0.5*fcen, 3.48*fcen)[0]
    k_o1 = ms.find_k(mp.NO_PARITY, fcen, 2, 2, mp.Vector3(1, 0, 0), 1e-4, 2.0*fcen, 0.5*fcen, 3.48*fcen)[0]
    lc_1 = 1.0 / (k_e1 - k_o1)

    diff = (lc_0 - lc_1) * 1e3
    dn = (k_e1 - k_o1) / fcen
    print(f"{t_nm:<12.1f} {dn:<14.6f} {lc_0:<14.5f} {lc_1:<14.5f} {diff:<14.2f}")
