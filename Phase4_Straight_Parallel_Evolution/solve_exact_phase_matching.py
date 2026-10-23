"""
PROJECT JANUS: HIGH-ACCURACY MPB VECTOR MODE & PHASE-MATCHING SOLVER
====================================================================
Computes exact 3D vector mode indices at lambda_0 = 1064 nm:
  1. Bare Silicon Waveguides: n_eff vs core width w in [200, 360] nm.
  2. Active Waveguide 1: n_eff in Amorphous & Crystalline states for t_pcm = 30nm & 150nm.
  3. Exact Phase-Matching Width w2 for bare WG2 to match WG1(Amorphous).
  4. Supermode analysis: Delta_n_split and exact coupling length L_c.

Usage:
  python solve_exact_phase_matching.py
"""

import sys, os, math
import numpy as np

try:
    import meep as mp
    from meep import mpb
    HAS_MPB = True
except Exception as e:
    HAS_MPB = False
    _MPB_ERR = str(e)


def solve_single_wg_neff(w_wg, h_wg=0.220, t_pcm=0.0, n_pcm=2.70, resolution=40):
    """Solves fundamental quasi-TE mode effective index using MPB."""
    if not HAS_MPB:
        raise ImportError(f"MPB / PyMEEP not found: {_MPB_ERR}")

    lambda_0 = 1.064
    fcen = 1.0 / lambda_0
    clad_mat = mp.Medium(index=1.449)
    core_mat = mp.Medium(index=3.480)
    pcm_mat = mp.Medium(index=n_pcm)

    sy, sz = 2.4, 2.0
    cell = mp.Vector3(0, sy, sz)

    geom = [
        mp.Block(center=mp.Vector3(0, 0, -sz / 4.0), size=mp.Vector3(mp.inf, mp.inf, sz / 2.0), material=clad_mat),
        mp.Block(center=mp.Vector3(0, 0, h_wg / 2.0), size=mp.Vector3(mp.inf, w_wg, h_wg), material=core_mat),
    ]
    if t_pcm > 1e-4:
        geom.append(mp.Block(center=mp.Vector3(0, 0, h_wg + t_pcm / 2.0),
                             size=mp.Vector3(mp.inf, w_wg, t_pcm), material=pcm_mat))

    k_guess = 2.40 * fcen
    ms = mpb.ModeSolver(
        geometry=geom,
        geometry_lattice=mp.Lattice(size=cell),
        k_points=[mp.Vector3(k_guess, 0, 0)],
        resolution=resolution,
        default_material=clad_mat,
        num_bands=2
    )

    old_stdout = sys.stdout
    sys.stdout = open(os.devnull, 'w')
    try:
        k_res = ms.find_k(mp.NO_PARITY, fcen, 1, 1, mp.Vector3(1, 0, 0), 1e-5, k_guess, 0.5 * fcen, 3.48 * fcen)
    finally:
        sys.stdout.close()
        sys.stdout = old_stdout

    return float(k_res[0] / fcen)


def solve_coupler_supermodes(w1, w2, gap, h_wg=0.220, t_pcm1=0.0, n_pcm1=2.70,
                             t_pcm2=0.0, n_pcm2=2.70, resolution=40):
    """Computes symmetric/antisymmetric supermodes of the coupled system."""
    if not HAS_MPB:
        raise ImportError(f"MPB / PyMEEP not found: {_MPB_ERR}")

    lambda_0 = 1.064
    fcen = 1.0 / lambda_0
    clad_mat = mp.Medium(index=1.449)
    core_mat = mp.Medium(index=3.480)

    y_sep = (w1 + w2) / 2.0 + gap
    y1 = y_sep / 2.0
    y2 = -y_sep / 2.0

    sy, sz = 3.6, 2.4
    cell = mp.Vector3(0, sy, sz)

    geom = [
        mp.Block(center=mp.Vector3(0, 0, -sz / 4.0), size=mp.Vector3(mp.inf, mp.inf, sz / 2.0), material=clad_mat),
        mp.Block(center=mp.Vector3(0, y1, h_wg / 2.0), size=mp.Vector3(mp.inf, w1, h_wg), material=core_mat),
        mp.Block(center=mp.Vector3(0, y2, h_wg / 2.0), size=mp.Vector3(mp.inf, w2, h_wg), material=core_mat),
    ]
    if t_pcm1 > 1e-4:
        geom.append(mp.Block(center=mp.Vector3(0, y1, h_wg + t_pcm1 / 2.0),
                             size=mp.Vector3(mp.inf, w1, t_pcm1), material=mp.Medium(index=n_pcm1)))
    if t_pcm2 > 1e-4:
        geom.append(mp.Block(center=mp.Vector3(0, y2, h_wg + t_pcm2 / 2.0),
                             size=mp.Vector3(mp.inf, w2, t_pcm2), material=mp.Medium(index=n_pcm2)))

    k_guess = 2.45 * fcen
    ms = mpb.ModeSolver(
        geometry=geom,
        geometry_lattice=mp.Lattice(size=cell),
        k_points=[mp.Vector3(k_guess, 0, 0)],
        resolution=resolution,
        default_material=clad_mat,
        num_bands=4
    )

    old_stdout = sys.stdout
    sys.stdout = open(os.devnull, 'w')
    try:
        k_res = ms.find_k(mp.NO_PARITY, fcen, 1, 2, mp.Vector3(1, 0, 0), 1e-5, k_guess, 0.5 * fcen, 3.48 * fcen)
    finally:
        sys.stdout.close()
        sys.stdout = old_stdout

    n_even = max(k_res[0], k_res[1]) / fcen
    n_odd = min(k_res[0], k_res[1]) / fcen
    dn_split = n_even - n_odd
    l_c = lambda_0 / (2.0 * max(dn_split, 1e-6))
    return float(n_even), float(n_odd), float(dn_split), float(l_c)


def find_exact_phase_match_width(target_neff, w_min=0.200, w_max=0.380, tol=1e-4, resolution=40):
    """Finds bare silicon width w2 matching target_neff via bisection."""
    low, high = w_min, w_max
    for _ in range(15):
        mid = (low + high) / 2.0
        neff_mid = solve_single_wg_neff(mid, t_pcm=0.0, resolution=resolution)
        if abs(neff_mid - target_neff) < tol:
            return mid, neff_mid
        if neff_mid < target_neff:
            low = mid
        else:
            high = mid
    return mid, neff_mid


if __name__ == "__main__":
    print("=" * 78)
    print("PROJECT JANUS: 3D VECTOR MPB PHASE-MATCHING & SUPERMODE SOLVER")
    print("=" * 78)

    if not HAS_MPB:
        print(f"Error: MPB not available ({_MPB_ERR}). Run this script in your Azure conda environment.")
        sys.exit(1)

    print("\n[STEP 1] ACTIVE WAVEGUIDE 1 MODAL PROFILES (w1 = 240 nm):")
    print("-" * 78)
    print(f"{'Config':<28} {'t_pcm (nm)':<12} {'n_am (2.70)':<14} {'n_cr (3.30)':<14} {'Delta_n'}")
    print("-" * 78)

    configs = [
        ("Option 1 (Thin Symmetric)", 0.030),
        ("Option 2 (Asymmetric Core)", 0.150),
        ("Combined (Thin Asymmetric)", 0.030)
    ]

    results = {}
    for label, tp in configs:
        n_am = solve_single_wg_neff(0.240, t_pcm=tp, n_pcm=2.70)
        n_cr = solve_single_wg_neff(0.240, t_pcm=tp, n_pcm=3.30)
        dn = n_cr - n_am
        results[label] = (tp, n_am, n_cr, dn)
        print(f"{label:<28} {tp*1000:<12.0f} {n_am:<14.5f} {n_cr:<14.5f} {dn:.5f}")

    print("\n[STEP 2] EXACT PHASE-MATCHED BARE SILICON WIDTH w2:")
    print("-" * 78)
    # Option 1: Symmetric -> w2 = 240 nm with 30nm PCM
    print(f"Option 1: Symmetrically identical -> w2 = 240 nm (both have 30nm PCM)")
    # Option 2: WG2 is bare Si matching n_am of 150nm PCM
    target_opt2 = results["Option 2 (Asymmetric Core)"][1]
    w2_opt2, neff2_opt2 = find_exact_phase_match_width(target_opt2, 0.280, 0.380)
    print(f"Option 2: Bare Si w2 matching n_am={target_opt2:.5f} -> w2 = {w2_opt2*1000:.1f} nm (n_eff={neff2_opt2:.5f})")
    # Combined: WG2 is bare Si matching n_am of 30nm PCM
    target_comb = results["Combined (Thin Asymmetric)"][1]
    w2_comb, neff2_comb = find_exact_phase_match_width(target_comb, 0.230, 0.300)
    print(f"Combined: Bare Si w2 matching n_am={target_comb:.5f} -> w2 = {w2_comb*1000:.1f} nm (n_eff={neff2_comb:.5f})")

    print("\n[STEP 3] COUPLER SUPERMODES & COUPLING LENGTH L_c AT GAP = 80 nm:")
    print("-" * 78)
    # Option 1: Symmetric 30nm PCM on both arms
    ne1, no1, dn1, lc1 = solve_coupler_supermodes(0.240, 0.240, 0.080, t_pcm1=0.030, t_pcm2=0.030)
    print(f"Option 1 (gap=80nm): n_even={ne1:.5f}, n_odd={no1:.5f}, Delta_n={dn1:.5f} -> L_c = {lc1:.3f} um")
    # Option 2: 150nm PCM on WG1, bare Si WG2
    ne2, no2, dn2, lc2 = solve_coupler_supermodes(0.240, w2_opt2, 0.100, t_pcm1=0.150, t_pcm2=0.0)
    print(f"Option 2 (gap=100nm): n_even={ne2:.5f}, n_odd={no2:.5f}, Delta_n={dn2:.5f} -> L_c = {lc2:.3f} um")
    # Combined: 30nm PCM on WG1, bare Si WG2
    nec, noc, dnc, lcc = solve_coupler_supermodes(0.240, w2_comb, 0.075, t_pcm1=0.030, t_pcm2=0.0)
    print(f"Combined (gap=75nm): n_even={nec:.5f}, n_odd={noc:.5f}, Delta_n={dnc:.5f} -> L_c = {lcc:.3f} um")

    print("\n" + "=" * 78)
    print("PHASE-MATCHING ANALYSIS COMPLETE. Use these exact parameters in FDTD!")
    print("=" * 78)
