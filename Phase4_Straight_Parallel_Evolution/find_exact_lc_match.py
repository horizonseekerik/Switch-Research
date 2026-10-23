import meep as mp
from meep import mpb

fcen = 1.0 / 1.064
clad = mp.Medium(index=1.449)
core = mp.Medium(index=3.480)
h_wg = 0.220
t_pcm = 0.190
w = 0.240
cell = mp.Vector3(0, 4.0, 2.5)

print("=" * 80)
print("FINE GAP SCAN FOR w = 240 nm")
print("=" * 80)

for gap_nm in [103.0, 103.2, 103.4, 103.5, 103.6, 103.8, 104.0]:
    gap = gap_nm * 1e-3
    d = gap + w
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
    print(f"Gap = {gap_nm:.2f} nm: Lc0 = {lc_0:.5f} um, Lc1 = {lc_1:.5f} um, Diff = {diff:+.2f} nm")
