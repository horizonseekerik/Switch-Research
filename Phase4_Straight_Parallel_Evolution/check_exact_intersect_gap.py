import math
import meep as mp
from meep import mpb

fcen = 1.0 / 1.064
clad = mp.Medium(index=1.449)
core = mp.Medium(index=3.480)
h_wg = 0.220
t_pcm = 0.190
w = 0.240
cell = mp.Vector3(0, 4.0, 2.5)

print("=" * 105)
print("INTERSECTION GAP FINDER IN MPB (Incorporating Fanout Ends L_ends0 = 0.70um, L_ends1 = 0.285um)")
print("=" * 105)
print(f"{'Gap (nm)':<10} {'Lc_0 (um)':<12} {'Lc_1 (um)':<12} {'Lc_phys0 (um)':<15} {'Lc_phys1 (um)':<15} {'Mismatch (nm)':<15}")
print("-" * 105)

for gap_nm in [103.0, 105.0, 107.0, 108.5, 110.0]:
    gap = gap_nm * 1e-3
    d = w + gap
    
    # State 0 supermodes (both waveguides have amorphous Sb2S3)
    geom_0 = [
        mp.Block(center=mp.Vector3(0, 0, -0.6), size=mp.Vector3(mp.inf, mp.inf, 1.2), material=clad),
        mp.Block(center=mp.Vector3(0, +d/2.0, h_wg/2.0), size=mp.Vector3(mp.inf, w, h_wg), material=core),
        mp.Block(center=mp.Vector3(0, +d/2.0, h_wg + t_pcm/2.0), size=mp.Vector3(mp.inf, w, t_pcm), material=mp.Medium(index=2.70)),
        mp.Block(center=mp.Vector3(0, -d/2.0, h_wg/2.0), size=mp.Vector3(mp.inf, w, h_wg), material=core),
        mp.Block(center=mp.Vector3(0, -d/2.0, h_wg + t_pcm/2.0), size=mp.Vector3(mp.inf, w, t_pcm), material=mp.Medium(index=2.70)),
    ]
    ms0 = mpb.ModeSolver(geometry=geom_0, geometry_lattice=mp.Lattice(size=cell),
                         k_points=[mp.Vector3(2.2*fcen, 0, 0)], resolution=32,
                         default_material=clad, num_bands=2)
    k_even_0 = ms0.find_k(mp.NO_PARITY, fcen, 1, 1, mp.Vector3(1, 0, 0), 1e-4, 2.0*fcen, 0.5*fcen, 3.48*fcen)[0]
    k_odd_0  = ms0.find_k(mp.NO_PARITY, fcen, 2, 2, mp.Vector3(1, 0, 0), 1e-4, 2.0*fcen, 0.5*fcen, 3.48*fcen)[0]
    lc_0 = 1.0 / (2.0 * (k_even_0 - k_odd_0))
    
    # State 1 supermodes (WG1 crystalline Sb2S3, WG2 amorphous Sb2S3)
    geom_1 = [
        mp.Block(center=mp.Vector3(0, 0, -0.6), size=mp.Vector3(mp.inf, mp.inf, 1.2), material=clad),
        mp.Block(center=mp.Vector3(0, +d/2.0, h_wg/2.0), size=mp.Vector3(mp.inf, w, h_wg), material=core),
        mp.Block(center=mp.Vector3(0, +d/2.0, h_wg + t_pcm/2.0), size=mp.Vector3(mp.inf, w, t_pcm), material=mp.Medium(index=3.30)),
        mp.Block(center=mp.Vector3(0, -d/2.0, h_wg/2.0), size=mp.Vector3(mp.inf, w, h_wg), material=core),
        mp.Block(center=mp.Vector3(0, -d/2.0, h_wg + t_pcm/2.0), size=mp.Vector3(mp.inf, w, t_pcm), material=mp.Medium(index=2.70)),
    ]
    ms0.geometry = geom_1
    k_even_1 = ms0.find_k(mp.NO_PARITY, fcen, 1, 1, mp.Vector3(1, 0, 0), 1e-4, 2.0*fcen, 0.5*fcen, 3.48*fcen)[0]
    k_odd_1  = ms0.find_k(mp.NO_PARITY, fcen, 2, 2, mp.Vector3(1, 0, 0), 1e-4, 2.0*fcen, 0.5*fcen, 3.48*fcen)[0]
    lc_1 = 1.0 / (k_even_1 - k_odd_1)
    
    lc_phys0 = lc_0 - 0.700
    lc_phys1 = lc_1 - 0.285
    mismatch_nm = (lc_phys0 - lc_phys1) * 1e3
    
    print(f"{gap_nm:<10.1f} {lc_0:<12.5f} {lc_1:<12.5f} {lc_phys0:<15.5f} {lc_phys1:<15.5f} {mismatch_nm:<15.2f}")
print("=" * 105)
