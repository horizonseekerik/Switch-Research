import math
import meep as mp
from meep import mpb

fcen = 1.0 / 1.064
clad = mp.Medium(index=1.449)
core = mp.Medium(index=3.480)
h_wg = 0.220
cell = mp.Vector3(0, 4.0, 3.0)

# First find exact w2 that matches w1=240nm with Amorphous PCM (tpcm=160nm)
geom_am1 = [
    mp.Block(center=mp.Vector3(0, 0, -0.75), size=mp.Vector3(mp.inf, mp.inf, 1.5), material=clad),
    mp.Block(center=mp.Vector3(0, 0, h_wg/2.0), size=mp.Vector3(mp.inf, 0.240, h_wg), material=core),
    mp.Block(center=mp.Vector3(0, 0, h_wg + 0.160/2.0), size=mp.Vector3(mp.inf, 0.240, 0.160), material=mp.Medium(index=2.70)),
]
ms1 = mpb.ModeSolver(geometry=geom_am1, geometry_lattice=mp.Lattice(size=cell),
                     k_points=[mp.Vector3(2.5*fcen, 0, 0)], resolution=30,
                     default_material=clad, num_bands=1)
k_am1 = ms1.find_k(mp.NO_PARITY, fcen, 1, 1, mp.Vector3(1, 0, 0), 1e-4, 2.5*fcen, 0.5*fcen, 3.48*fcen)
target_nam = k_am1[0] / fcen
print(f"Target WG1 (Amorphous) n_eff = {target_nam:.5f}")

# Crystalline single WG1 n_eff
geom_cr1 = [
    mp.Block(center=mp.Vector3(0, 0, -0.75), size=mp.Vector3(mp.inf, mp.inf, 1.5), material=clad),
    mp.Block(center=mp.Vector3(0, 0, h_wg/2.0), size=mp.Vector3(mp.inf, 0.240, h_wg), material=core),
    mp.Block(center=mp.Vector3(0, 0, h_wg + 0.160/2.0), size=mp.Vector3(mp.inf, 0.240, 0.160), material=mp.Medium(index=3.30)),
]
ms1.geometry = geom_cr1
k_cr1 = ms1.find_k(mp.NO_PARITY, fcen, 1, 1, mp.Vector3(1, 0, 0), 1e-4, 2.5*fcen, 0.5*fcen, 3.48*fcen)
target_ncr = k_cr1[0] / fcen
print(f"Target WG1 (Crystalline) n_eff = {target_ncr:.5f}")
delta_beta = 2 * math.pi * fcen * (target_ncr - target_nam)
print(f"Delta beta (cr - am) = {delta_beta:.4f} rad/um")

# Fine scan for exact w2
best_w2 = 0.289
min_diff = 1e9
for w_test in [0.288, 0.289, 0.290]:
    geom_test = [
        mp.Block(center=mp.Vector3(0, 0, -0.75), size=mp.Vector3(mp.inf, mp.inf, 1.5), material=clad),
        mp.Block(center=mp.Vector3(0, 0, h_wg/2.0), size=mp.Vector3(mp.inf, w_test, h_wg), material=core),
    ]
    ms1.geometry = geom_test
    k_t = ms1.find_k(mp.NO_PARITY, fcen, 1, 1, mp.Vector3(1, 0, 0), 1e-4, 2.5*fcen, 0.5*fcen, 3.48*fcen)
    n_t = k_t[0] / fcen
    diff = abs(n_t - target_nam)
    print(f"w_bare = {w_test*1000:.1f} nm, n_eff = {n_t:.5f}, diff = {diff:.5f}")
    if diff < min_diff:
        min_diff = diff
        best_w2 = w_test

print(f"\nBest matching bare Si width w2 = {best_w2*1000:.1f} nm\n")

# Now sweep gap between WG1 (with PCM) and WG2 (bare Si)
w1 = 0.240
w2 = best_w2
t_pcm = 0.160

print(f"{'gap (nm)':<10} {'dk_am (rad/um)':<16} {'dk_cr (rad/um)':<16} {'Ratio dk_cr/dk_am':<20} {'Lc (um)':<10}")
print("-" * 72)

for gap_nm in [80, 90, 100, 110, 120]:
    gap = gap_nm / 1000.0
    d = gap + (w1 + w2) / 2.0
    y1 = +d / 2.0
    y2 = -d / 2.0

    # Amorphous supermodes
    geom_am_coupler = [
        mp.Block(center=mp.Vector3(0, 0, -0.75), size=mp.Vector3(mp.inf, mp.inf, 1.5), material=clad),
        # WG1 (with Amorphous PCM)
        mp.Block(center=mp.Vector3(0, y1, h_wg/2.0), size=mp.Vector3(mp.inf, w1, h_wg), material=core),
        mp.Block(center=mp.Vector3(0, y1, h_wg + t_pcm/2.0), size=mp.Vector3(mp.inf, w1, t_pcm), material=mp.Medium(index=2.70)),
        # WG2 (Bare Si, NO PCM)
        mp.Block(center=mp.Vector3(0, y2, h_wg/2.0), size=mp.Vector3(mp.inf, w2, h_wg), material=core),
    ]

    ms_c = mpb.ModeSolver(geometry=geom_am_coupler, geometry_lattice=mp.Lattice(size=cell),
                          k_points=[mp.Vector3(2.5*fcen, 0, 0)], resolution=30,
                          default_material=clad, num_bands=2)
    k_am_even = ms_c.find_k(mp.NO_PARITY, fcen, 1, 1, mp.Vector3(1, 0, 0), 1e-4, 2.5*fcen, 0.5*fcen, 3.48*fcen)[0]
    k_am_odd  = ms_c.find_k(mp.NO_PARITY, fcen, 2, 2, mp.Vector3(1, 0, 0), 1e-4, 2.5*fcen, 0.5*fcen, 3.48*fcen)[0]
    dk_am = (k_am_even - k_am_odd) * 2 * math.pi
    Lc = math.pi / dk_am

    # Crystalline supermodes
    geom_cr_coupler = [
        mp.Block(center=mp.Vector3(0, 0, -0.75), size=mp.Vector3(mp.inf, mp.inf, 1.5), material=clad),
        # WG1 (with Crystalline PCM)
        mp.Block(center=mp.Vector3(0, y1, h_wg/2.0), size=mp.Vector3(mp.inf, w1, h_wg), material=core),
        mp.Block(center=mp.Vector3(0, y1, h_wg + t_pcm/2.0), size=mp.Vector3(mp.inf, w1, t_pcm), material=mp.Medium(index=3.30)),
        # WG2 (Bare Si, NO PCM)
        mp.Block(center=mp.Vector3(0, y2, h_wg/2.0), size=mp.Vector3(mp.inf, w2, h_wg), material=core),
    ]
    ms_c.geometry = geom_cr_coupler
    k_cr_even = ms_c.find_k(mp.NO_PARITY, fcen, 1, 1, mp.Vector3(1, 0, 0), 1e-4, 2.5*fcen, 0.5*fcen, 3.48*fcen)[0]
    k_cr_odd  = ms_c.find_k(mp.NO_PARITY, fcen, 2, 2, mp.Vector3(1, 0, 0), 1e-4, 2.5*fcen, 0.5*fcen, 3.48*fcen)[0]
    dk_cr = (k_cr_even - k_cr_odd) * 2 * math.pi

    ratio = dk_cr / dk_am
    print(f"{gap_nm:<10.0f} {dk_am:<16.4f} {dk_cr:<16.4f} {ratio:<20.4f} {Lc:<10.3f}")
