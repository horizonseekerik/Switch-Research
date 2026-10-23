import sys, os, math
import numpy as np
import meep as mp
from meep import mpb

def solve_candidate3_supermodes(gap, w1=0.240, w2=0.240, t_pcm=0.200, resolution=32):
    lambda_0 = 1.064
    fcen = 1.0 / lambda_0
    h = 0.220
    sy = 3.2
    sz = 2.4
    cell = mp.Vector3(0, sy, sz)
    clad = mp.Medium(index=1.449)
    core = mp.Medium(index=3.480)
    
    y_sep = (w1 + w2) / 2.0 + gap
    y1 = y_sep / 2.0
    y2 = -y_sep / 2.0
    
    # State 0: Both Amorphous (n = 2.70)
    geom_am = [
        mp.Block(center=mp.Vector3(0, 0, -sz/4.0), size=mp.Vector3(mp.inf, mp.inf, sz/2.0), material=clad),
        mp.Block(center=mp.Vector3(0, y1, h/2.0), size=mp.Vector3(mp.inf, w1, h), material=core),
        mp.Block(center=mp.Vector3(0, y2, h/2.0), size=mp.Vector3(mp.inf, w2, h), material=core),
        mp.Block(center=mp.Vector3(0, y1, h + t_pcm/2.0), size=mp.Vector3(mp.inf, w1, t_pcm), material=mp.Medium(index=2.70)),
        mp.Block(center=mp.Vector3(0, y2, h + t_pcm/2.0), size=mp.Vector3(mp.inf, w2, t_pcm), material=mp.Medium(index=2.70))
    ]
    
    k_guess = 2.20 * fcen
    ms_am = mpb.ModeSolver(
        geometry=geom_am,
        geometry_lattice=mp.Lattice(size=cell),
        k_points=[mp.Vector3(k_guess, 0, 0)],
        resolution=resolution,
        default_material=clad,
        num_bands=4
    )
    k1_am = ms_am.find_k(mp.NO_PARITY, fcen, 1, 1, mp.Vector3(1, 0, 0), 1e-4, k_guess, 0.5*fcen, 3.48*fcen)[0] / fcen
    k2_am = ms_am.find_k(mp.NO_PARITY, fcen, 2, 2, mp.Vector3(1, 0, 0), 1e-4, k_guess, 0.5*fcen, 3.48*fcen)[0] / fcen
    
    dn_am = k1_am - k2_am
    lc_0 = lambda_0 / (2.0 * dn_am) if dn_am > 0 else 1e9
    kappa = math.pi / (2.0 * lc_0)
    
    # State 1: WG1 Crystalline (n = 3.30), WG2 Amorphous (n = 2.70)
    geom_cr = [
        mp.Block(center=mp.Vector3(0, 0, -sz/4.0), size=mp.Vector3(mp.inf, mp.inf, sz/2.0), material=clad),
        mp.Block(center=mp.Vector3(0, y1, h/2.0), size=mp.Vector3(mp.inf, w1, h), material=core),
        mp.Block(center=mp.Vector3(0, y2, h/2.0), size=mp.Vector3(mp.inf, w2, h), material=core),
        mp.Block(center=mp.Vector3(0, y1, h + t_pcm/2.0), size=mp.Vector3(mp.inf, w1, t_pcm), material=mp.Medium(index=3.30)),
        mp.Block(center=mp.Vector3(0, y2, h + t_pcm/2.0), size=mp.Vector3(mp.inf, w2, t_pcm), material=mp.Medium(index=2.70))
    ]
    
    k_guess_cr = 2.28 * fcen
    ms_cr = mpb.ModeSolver(
        geometry=geom_cr,
        geometry_lattice=mp.Lattice(size=cell),
        k_points=[mp.Vector3(k_guess_cr, 0, 0)],
        resolution=resolution,
        default_material=clad,
        num_bands=4
    )
    k1_cr = ms_cr.find_k(mp.NO_PARITY, fcen, 1, 1, mp.Vector3(1, 0, 0), 1e-4, k_guess_cr, 0.5*fcen, 3.48*fcen)[0] / fcen
    k2_cr = ms_cr.find_k(mp.NO_PARITY, fcen, 2, 2, mp.Vector3(1, 0, 0), 1e-4, k_guess_cr, 0.5*fcen, 3.48*fcen)[0] / fcen
    
    dn_cr = k1_cr - k2_cr
    s_cr = (math.pi / lambda_0) * dn_cr
    lc_1 = lambda_0 / dn_cr if dn_cr > 0 else 1e9
    
    ratio = s_cr / kappa if kappa > 0 else 0.0
    offset = lc_1 - lc_0
    
    return {
        "gap": gap, "k1_am": k1_am, "k2_am": k2_am, "dn_am": dn_am, "lc_0": lc_0,
        "k1_cr": k1_cr, "k2_cr": k2_cr, "dn_cr": dn_cr, "lc_1": lc_1,
        "ratio_S_kappa": ratio, "offset": offset
    }

if __name__ == "__main__":
    print("=" * 115)
    print("  EXACT MPB RIGOROUS SUPERMODE ANALYSIS FOR CANDIDATE 3")
    print("  Evaluating gap dependence of L_c(0) vs L_c(1) and S/kappa")
    print("=" * 115)
    print(f"{'Gap (nm)':<10} {'dn_am':<12} {'L_c(0) (um)':<14} {'dn_cr':<12} {'L_c(1) (um)':<14} {'S/kappa':<12} {'Offset(um)':<12}")
    print("-" * 115)
    
    for g_nm in [80, 85, 89, 94, 98, 102, 106, 110, 115, 120]:
        g = g_nm * 1e-3
        res = solve_candidate3_supermodes(g, resolution=28)
        print(f"{g_nm:<10.1f} {res['dn_am']:<12.6f} {res['lc_0']:<14.6f} {res['dn_cr']:<12.6f} {res['lc_1']:<14.6f} {res['ratio_S_kappa']:<12.6f} {res['offset']:<12.6f}")
        sys.stdout.flush()
    print("=" * 115)
