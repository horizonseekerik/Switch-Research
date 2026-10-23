import sys, os, math
import numpy as np
import meep as mp
from meep import mpb

def solve_mode_quick(platform="sin", w=0.450, h=0.300, t_pcm=0.050, n_pcm=2.70):
    lambda_0 = 1.064
    fcen = 1.0 / lambda_0
    clad = mp.Medium(index=1.449)
    if platform == "sin":
        core = mp.Medium(index=2.016)
        k_guess = 1.75 * fcen
    else:
        core = mp.Medium(index=3.480)
        k_guess = 2.45 * fcen
    pcm = mp.Medium(index=n_pcm)

    sy, sz = 3.0, 3.0
    cell = mp.Vector3(0, sy, sz)
    geom = [
        mp.Block(center=mp.Vector3(0, 0, -sz/4.0), size=mp.Vector3(mp.inf, mp.inf, sz/2.0), material=clad),
        mp.Block(center=mp.Vector3(0, 0, h/2.0), size=mp.Vector3(mp.inf, w, h), material=core),
        mp.Block(center=mp.Vector3(0, 0, h + t_pcm/2.0), size=mp.Vector3(mp.inf, w, t_pcm), material=pcm),
    ]
    ms = mpb.ModeSolver(geometry=geom, geometry_lattice=mp.Lattice(size=cell),
                        k_points=[mp.Vector3(k_guess, 0, 0)], resolution=30,
                        default_material=clad, num_bands=2)
    old_stdout = sys.stdout
    sys.stdout = open(os.devnull, 'w')
    try:
        k_res = ms.find_k(mp.NO_PARITY, fcen, 1, 1, mp.Vector3(1, 0, 0), 1e-4, k_guess, 0.5*fcen, 3.48*fcen)
    finally:
        sys.stdout.close()
        sys.stdout = old_stdout
    return k_res[0] / fcen

print("=== PLATFORM COMPARISON: Delta_n_eff vs t_pcm ===")
print("Target: L_c <= 3.50 um requires Delta_n_eff >= 0.263\n")

print("1. SILICON NITRIDE (Si3N4: w=450nm, h=300nm):")
for t in [0.030, 0.050, 0.080, 0.100, 0.150]:
    na = solve_mode_quick("sin", 0.450, 0.300, t, 2.70)
    nc = solve_mode_quick("sin", 0.450, 0.300, t, 3.30)
    dn = nc - na
    lc = math.sqrt(3.0) * 1.064 / (2.0 * dn)
    print(f"  t_pcm = {t*1000:3.0f} nm: n_am = {na:.4f}, n_cr = {nc:.4f}, Delta_n = {dn:.4f} -> L_c = {lc:.3f} um")

print("\n2. SILICON (c-Si: w=250nm, h=220nm):")
for t in [0.030, 0.050, 0.080, 0.100, 0.150]:
    na = solve_mode_quick("si", 0.250, 0.220, t, 2.70)
    nc = solve_mode_quick("si", 0.250, 0.220, t, 3.30)
    dn = nc - na
    lc = math.sqrt(3.0) * 1.064 / (2.0 * dn)
    print(f"  t_pcm = {t*1000:3.0f} nm: n_am = {na:.4f}, n_cr = {nc:.4f}, Delta_n = {dn:.4f} -> L_c = {lc:.3f} um")
