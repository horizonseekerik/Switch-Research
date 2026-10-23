import sys, os, math
import numpy as np
import meep as mp
from meep import mpb

lambda_0 = 1.064
fcen = 1.0 / lambda_0
h_wg = 0.220
clad_mat = mp.Medium(index=1.449)
core_mat = mp.Medium(index=3.480)

def solve_supermodes(w_wg, gap, t_pcm, resolution=32):
    y_sep = w_wg + gap
    y1 = y_sep / 2.0
    y2 = -y_sep / 2.0
    sy, sz = 3.0, 2.8
    cell = mp.Vector3(0, sy, sz)
    
    # State 0: Both Amorphous (n = 2.70)
    geom_am = [
        mp.Block(center=mp.Vector3(0, 0, -sz/4.0), size=mp.Vector3(mp.inf, mp.inf, sz/2.0), material=clad_mat),
        mp.Block(center=mp.Vector3(0, y1, h_wg/2.0), size=mp.Vector3(mp.inf, w_wg, h_wg), material=core_mat),
        mp.Block(center=mp.Vector3(0, y2, h_wg/2.0), size=mp.Vector3(mp.inf, w_wg, h_wg), material=core_mat),
        mp.Block(center=mp.Vector3(0, y1, h_wg + t_pcm/2.0), size=mp.Vector3(mp.inf, w_wg, t_pcm), material=mp.Medium(index=2.70)),
        mp.Block(center=mp.Vector3(0, y2, h_wg + t_pcm/2.0), size=mp.Vector3(mp.inf, w_wg, t_pcm), material=mp.Medium(index=2.70)),
    ]
    k_guess = 2.50 * fcen
    ms_am = mpb.ModeSolver(geometry=geom_am, geometry_lattice=mp.Lattice(size=cell),
                           k_points=[mp.Vector3(k_guess, 0, 0)], resolution=resolution,
                           default_material=clad_mat, num_bands=4)
    old_stdout = sys.stdout
    sys.stdout = open(os.devnull, 'w')
    try:
        k_am = ms_am.find_k(mp.NO_PARITY, fcen, 1, 2, mp.Vector3(1, 0, 0), 1e-4, k_guess, 0.5*fcen, 3.48*fcen)
    finally:
        sys.stdout.close()
        sys.stdout = old_stdout
    
    dn_am = abs(k_am[0] - k_am[1]) / fcen
    lc_0 = lambda_0 / (2.0 * dn_am) if dn_am > 0 else 1e9
    kappa_am = math.pi / (2.0 * lc_0)

    # State 1: WG1 Crystalline (n = 3.30), WG2 Amorphous (n = 2.70)
    geom_cr = [
        mp.Block(center=mp.Vector3(0, 0, -sz/4.0), size=mp.Vector3(mp.inf, mp.inf, sz/2.0), material=clad_mat),
        mp.Block(center=mp.Vector3(0, y1, h_wg/2.0), size=mp.Vector3(mp.inf, w_wg, h_wg), material=core_mat),
        mp.Block(center=mp.Vector3(0, y2, h_wg/2.0), size=mp.Vector3(mp.inf, w_wg, h_wg), material=core_mat),
        mp.Block(center=mp.Vector3(0, y1, h_wg + t_pcm/2.0), size=mp.Vector3(mp.inf, w_wg, t_pcm), material=mp.Medium(index=3.30)),
        mp.Block(center=mp.Vector3(0, y2, h_wg + t_pcm/2.0), size=mp.Vector3(mp.inf, w_wg, t_pcm), material=mp.Medium(index=2.70)),
    ]
    k_guess_cr = 2.65 * fcen
    ms_cr = mpb.ModeSolver(geometry=geom_cr, geometry_lattice=mp.Lattice(size=cell),
                           k_points=[mp.Vector3(k_guess_cr, 0, 0)], resolution=resolution,
                           default_material=clad_mat, num_bands=4)
    old_stdout = sys.stdout
    sys.stdout = open(os.devnull, 'w')
    try:
        k_cr = ms_cr.find_k(mp.NO_PARITY, fcen, 1, 2, mp.Vector3(1, 0, 0), 1e-4, k_guess_cr, 0.5*fcen, 3.48*fcen)
    finally:
        sys.stdout.close()
        sys.stdout = old_stdout
    
    dn_cr = abs(k_cr[0] - k_cr[1]) / fcen
    lc_1 = lambda_0 / dn_cr if dn_cr > 0 else 1e9
    s_cr = (math.pi / lambda_0) * dn_cr
    ratio = s_cr / kappa_am if kappa_am > 0 else 0.0
    offset = lc_1 - lc_0

    return {
        "gap": gap, "dn_am": dn_am, "lc_0": lc_0, "dn_cr": dn_cr, "lc_1": lc_1,
        "s_cr": s_cr, "kappa_am": kappa_am, "ratio": ratio, "offset": offset
    }

print("="*110)
print("SUPERMODE ANALYSIS FOR THICK PCM (t_pcm = 300 nm & 320 nm)")
print("="*110)
print(f"{'t_pcm(nm)':<10} {'Gap(nm)':<10} {'dn_am':<12} {'L_c(0)(um)':<14} {'dn_cr':<12} {'L_c(1)(um)':<14} {'S/kappa':<12} {'Offset(um)':<12}")
print("-" * 110)

w = 0.240
for t_nm in [300, 320]:
    t = t_nm * 1e-3
    for g_nm in [80, 85, 90, 95, 100]:
        g = g_nm * 1e-3
        r = solve_supermodes(w, g, t)
        print(f"{t_nm:<10} {g_nm:<10} {r['dn_am']:<12.6f} {r['lc_0']:<14.6f} {r['dn_cr']:<12.6f} {r['lc_1']:<14.6f} {r['ratio']:<12.6f} {r['offset']:<12.6f}")
