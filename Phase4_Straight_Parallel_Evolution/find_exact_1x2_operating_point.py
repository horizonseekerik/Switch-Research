import sys, os, math
import numpy as np
import meep as mp
from meep import mpb

def solve_mode(w, h=0.220, t_pcm=0.150, n_pcm=2.70, resolution=40):
    lambda_0 = 1.064
    fcen = 1.0 / lambda_0
    clad = mp.Medium(index=1.449)
    core = mp.Medium(index=3.480)
    pcm = mp.Medium(index=n_pcm)

    sy, sz = 2.5, 2.5
    cell = mp.Vector3(0, sy, sz)
    geom = [
        mp.Block(center=mp.Vector3(0, 0, -sz/4.0), size=mp.Vector3(mp.inf, mp.inf, sz/2.0), material=clad),
        mp.Block(center=mp.Vector3(0, 0, h/2.0), size=mp.Vector3(mp.inf, w, h), material=core),
        mp.Block(center=mp.Vector3(0, 0, h + t_pcm/2.0), size=mp.Vector3(mp.inf, w, t_pcm), material=pcm),
    ]
    k_guess = 2.45 * fcen
    ms = mpb.ModeSolver(geometry=geom, geometry_lattice=mp.Lattice(size=cell),
                        k_points=[mp.Vector3(k_guess, 0, 0)], resolution=resolution,
                        default_material=clad, num_bands=2)
    old_stdout = sys.stdout
    sys.stdout = open(os.devnull, 'w')
    try:
        k_res = ms.find_k(mp.NO_PARITY, fcen, 1, 1, mp.Vector3(1, 0, 0), 1e-4, k_guess, 0.5*fcen, 3.48*fcen)
    finally:
        sys.stdout.close()
        sys.stdout = old_stdout
    return k_res[0] / fcen

def solve_supermodes(w, gap, h=0.220, t_pcm=0.150, n_pcm=2.70, resolution=40):
    lambda_0 = 1.064
    fcen = 1.0 / lambda_0
    clad = mp.Medium(index=1.449)
    core = mp.Medium(index=3.480)
    pcm = mp.Medium(index=n_pcm)

    y_sep = w + gap
    y1, y2 = y_sep / 2.0, -y_sep / 2.0
    sy, sz = 3.0, 2.5
    cell = mp.Vector3(0, sy, sz)
    geom = [
        mp.Block(center=mp.Vector3(0, 0, -sz/4.0), size=mp.Vector3(mp.inf, mp.inf, sz/2.0), material=clad),
        mp.Block(center=mp.Vector3(0, y1, h/2.0), size=mp.Vector3(mp.inf, w, h), material=core),
        mp.Block(center=mp.Vector3(0, y2, h/2.0), size=mp.Vector3(mp.inf, w, h), material=core),
        mp.Block(center=mp.Vector3(0, y1, h + t_pcm/2.0), size=mp.Vector3(mp.inf, w, t_pcm), material=pcm),
        mp.Block(center=mp.Vector3(0, y2, h + t_pcm/2.0), size=mp.Vector3(mp.inf, w, t_pcm), material=pcm),
    ]
    k_guess = 2.45 * fcen
    ms = mpb.ModeSolver(geometry=geom, geometry_lattice=mp.Lattice(size=cell),
                        k_points=[mp.Vector3(k_guess, 0, 0)], resolution=resolution,
                        default_material=clad, num_bands=4)
    old_stdout = sys.stdout
    sys.stdout = open(os.devnull, 'w')
    try:
        k_res = ms.find_k(mp.NO_PARITY, fcen, 1, 2, mp.Vector3(1, 0, 0), 1e-4, k_guess, 0.5*fcen, 3.48*fcen)
    finally:
        sys.stdout.close()
        sys.stdout = old_stdout
    ne = max(k_res[0], k_res[1]) / fcen
    no = min(k_res[0], k_res[1]) / fcen
    return ne, no

if __name__ == '__main__':
    w = 0.240
    t_pcm = 0.150
    lambda_0 = 1.064
    print("=" * 70)
    print("EXACT MPB 1x2 SWITCH OPERATING POINT FINDER")
    print("=" * 70)

    print("Step 1: Calculating Single Waveguide Modal Indices...")
    n_am = solve_mode(w, t_pcm=t_pcm, n_pcm=2.70)
    n_cr = solve_mode(w, t_pcm=t_pcm, n_pcm=3.30)
    dn_eff = n_cr - n_am
    print(f"  n_eff(Amorphous)   = {n_am:.5f}")
    print(f"  n_eff(Crystalline) = {n_cr:.5f}")
    print(f"  delta_n_eff        = {dn_eff:.5f}")

    target_splitting = dn_eff / math.sqrt(3.0)
    optimal_lc = math.sqrt(3.0) * lambda_0 / (2.0 * dn_eff)
    print(f"\nStep 2: Switching Condition (u = sqrt(3) = 1.73205):")
    print(f"  Required Supermode Splitting (n_even - n_odd) = {target_splitting:.5f}")
    print(f"  Analytical Coupling Length L_c                = {optimal_lc:.3f} um")

    print(f"\nStep 3: Sweeping Coupling Gap for Exact Splitting:")
    print(f"{'Gap (nm)':<10} {'n_even':<10} {'n_odd':<10} {'dn_split':<10} {'L_c (um)':<10} {'u':<8} {'P_bar(Cr)%':<12} {'XT(Cr) dB':<12}")
    print("-" * 82)
    gaps = np.arange(0.034, 0.046, 0.002)
    best_gap = None
    best_diff = 1e9

    for g in gaps:
        ne, no = solve_supermodes(w, g, t_pcm=t_pcm, n_pcm=2.70)
        split = ne - no
        lc = lambda_0 / (2.0 * max(split, 1e-6))
        u = dn_eff / max(split, 1e-6)
        # Power in bar port for crystalline state at length lc:
        # P_cross = 1/(1+u^2) * sin^2(pi/2 * sqrt(1+u^2))
        arg = 0.5 * math.pi * math.sqrt(1.0 + u**2)
        p_cross_cr = (1.0 / (1.0 + u**2)) * (math.sin(arg)**2)
        p_bar_cr = 1.0 - p_cross_cr
        xt_cr = 10.0 * math.log10(max(p_cross_cr / max(p_bar_cr, 1e-9), 1e-9))
        print(f"{g*1000:<10.1f} {ne:<10.5f} {no:<10.5f} {split:<10.5f} {lc:<10.3f} {u:<8.3f} {p_bar_cr*100:<11.2f}% {xt_cr:<12.2f}")
        if abs(split - target_splitting) < best_diff:
            best_diff = abs(split - target_splitting)
            best_gap = g

    print("=" * 82)
    print(f"[OPTIMAL OPERATING POINT]: Gap = {best_gap*1000:.1f} nm, L_c = {optimal_lc:.3f} um")
