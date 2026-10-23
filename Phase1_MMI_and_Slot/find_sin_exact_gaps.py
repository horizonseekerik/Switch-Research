import sys, os, math
import numpy as np
import meep as mp
from meep import mpb

def solve_sin_supermodes(w=0.450, gap=0.100, h=0.300, t_pcm=0.080, n_pcm=2.70, resolution=40):
    lambda_0 = 1.064
    fcen = 1.0 / lambda_0
    clad = mp.Medium(index=1.449)
    core = mp.Medium(index=2.016)
    pcm = mp.Medium(index=n_pcm)

    y_sep = w + gap
    y1, y2 = y_sep / 2.0, -y_sep / 2.0
    sy, sz = 4.0, 3.0
    cell = mp.Vector3(0, sy, sz)
    geom = [
        mp.Block(center=mp.Vector3(0, 0, -sz/4.0), size=mp.Vector3(mp.inf, mp.inf, sz/2.0), material=clad),
        mp.Block(center=mp.Vector3(0, y1, h/2.0), size=mp.Vector3(mp.inf, w, h), material=core),
        mp.Block(center=mp.Vector3(0, y2, h/2.0), size=mp.Vector3(mp.inf, w, h), material=core),
        mp.Block(center=mp.Vector3(0, y1, h + t_pcm/2.0), size=mp.Vector3(mp.inf, w, t_pcm), material=pcm),
        mp.Block(center=mp.Vector3(0, y2, h + t_pcm/2.0), size=mp.Vector3(mp.inf, w, t_pcm), material=pcm),
    ]

    k_guess = 1.80 * fcen
    ms = mpb.ModeSolver(geometry=geom, geometry_lattice=mp.Lattice(size=cell),
                        k_points=[mp.Vector3(k_guess, 0, 0)], resolution=resolution,
                        default_material=clad, num_bands=4)
    old_stdout = sys.stdout
    sys.stdout = open(os.devnull, 'w')
    try:
        k_res = ms.find_k(mp.NO_PARITY, fcen, 1, 2, mp.Vector3(1, 0, 0), 1e-4, k_guess, 0.5*fcen, 2.5*fcen)
    finally:
        sys.stdout.close()
        sys.stdout = old_stdout
    ne = max(k_res[0], k_res[1]) / fcen
    no = min(k_res[0], k_res[1]) / fcen
    return ne, no

lambda_0 = 1.064
for t_pcm, dn_eff in [(0.080, 0.2788), (0.100, 0.3557)]:
    target_split = dn_eff / math.sqrt(3.0)
    optimal_lc = math.sqrt(3.0) * lambda_0 / (2.0 * dn_eff)
    print(f"\n==========================================================================")
    print(f"Si3N4 with t_pcm = {t_pcm*1000:.0f} nm: dn_eff = {dn_eff:.4f}")
    print(f"Target Splitting = {target_split:.5f}, Target L_c = {optimal_lc:.3f} um")
    print(f"==========================================================================")
    print(f"{'Gap (nm)':<10} {'n_even':<10} {'n_odd':<10} {'dn_split':<10} {'L_c (um)':<10} {'u':<8} {'P_bar(Cr)%':<12} {'XT(Cr) dB':<12}")
    print("-" * 82)
    for g in [0.030, 0.040, 0.050, 0.060, 0.070, 0.080]:
        ne, no = solve_sin_supermodes(0.450, g, 0.300, t_pcm=t_pcm, n_pcm=2.70)
        split = ne - no
        lc = lambda_0 / (2.0 * max(split, 1e-6))
        u = dn_eff / max(split, 1e-6)
        arg = 0.5 * math.pi * math.sqrt(1.0 + u**2)
        p_cross_cr = (1.0 / (1.0 + u**2)) * (math.sin(arg)**2)
        p_bar_cr = 1.0 - p_cross_cr
        xt_cr = 10.0 * math.log10(max(p_cross_cr / max(p_bar_cr, 1e-9), 1e-9))
        print(f"{g*1000:<10.1f} {ne:<10.5f} {no:<10.5f} {split:<10.5f} {lc:<10.3f} {u:<8.3f} {p_bar_cr*100:<11.2f}% {xt_cr:<12.2f}")
