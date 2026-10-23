"""
Check supermode overlap integrals in MPB to understand directional coupler transfer efficiency
"""
import sys, os, math
import numpy as np
import meep as mp
from meep import mpb

lambda_0 = 1.064
fcen = 1.0 / lambda_0
h_wg = 0.220
clad_mat = mp.Medium(index=1.449)
core_mat = mp.Medium(index=3.480)

def test_overlap(w_wg=0.240, gap=0.100, t_pcm=0.150, n_pcm=2.70, resolution=40):
    pcm_mat = mp.Medium(index=n_pcm)
    y_sep = w_wg + gap
    y1 = y_sep / 2.0
    y2 = -y_sep / 2.0
    sy, sz = 3.0, 2.5
    cell = mp.Vector3(0, sy, sz)

    # 1. Isolated WG1
    geom_iso = [
        mp.Block(center=mp.Vector3(0, 0, -sz/4.0), size=mp.Vector3(mp.inf, mp.inf, sz/2.0), material=clad_mat),
        mp.Block(center=mp.Vector3(0, y1, h_wg/2.0), size=mp.Vector3(mp.inf, w_wg, h_wg), material=core_mat),
        mp.Block(center=mp.Vector3(0, y1, h_wg + t_pcm/2.0), size=mp.Vector3(mp.inf, w_wg, t_pcm), material=pcm_mat)
    ]
    k_guess = 2.40 * fcen
    ms_iso = mpb.ModeSolver(geometry=geom_iso, geometry_lattice=mp.Lattice(size=cell),
                            k_points=[mp.Vector3(k_guess, 0, 0)], resolution=resolution,
                            default_material=clad_mat, num_bands=2)
    old_stdout = sys.stdout
    sys.stdout = open(os.devnull, 'w')
    try:
        k_iso = ms_iso.find_k(mp.NO_PARITY, fcen, 1, 1, mp.Vector3(1, 0, 0), 1e-4, k_guess, 0.5*fcen, 3.48*fcen)
        efield_iso = ms_iso.get_efield(1) # array of [y, z, component]
    finally:
        sys.stdout.close()
        sys.stdout = old_stdout

    # 2. Coupled System (Both WGs)
    geom_cpl = [
        mp.Block(center=mp.Vector3(0, 0, -sz/4.0), size=mp.Vector3(mp.inf, mp.inf, sz/2.0), material=clad_mat),
        mp.Block(center=mp.Vector3(0, y1, h_wg/2.0), size=mp.Vector3(mp.inf, w_wg, h_wg), material=core_mat),
        mp.Block(center=mp.Vector3(0, y2, h_wg/2.0), size=mp.Vector3(mp.inf, w_wg, h_wg), material=core_mat),
        mp.Block(center=mp.Vector3(0, y1, h_wg + t_pcm/2.0), size=mp.Vector3(mp.inf, w_wg, t_pcm), material=pcm_mat),
        mp.Block(center=mp.Vector3(0, y2, h_wg + t_pcm/2.0), size=mp.Vector3(mp.inf, w_wg, t_pcm), material=pcm_mat),
    ]
    ms_cpl = mpb.ModeSolver(geometry=geom_cpl, geometry_lattice=mp.Lattice(size=cell),
                            k_points=[mp.Vector3(k_guess, 0, 0)], resolution=resolution,
                            default_material=clad_mat, num_bands=4)
    old_stdout = sys.stdout
    sys.stdout = open(os.devnull, 'w')
    try:
        k_cpl = ms_cpl.find_k(mp.NO_PARITY, fcen, 1, 2, mp.Vector3(1, 0, 0), 1e-4, k_guess, 0.5*fcen, 3.48*fcen)
        efield_even = ms_cpl.get_efield(1)
        efield_odd = ms_cpl.get_efield(2)
    finally:
        sys.stdout.close()
        sys.stdout = old_stdout

    # Calculate overlaps (using Ey which is dominant for TE)
    # E-field array shape: (Nx, Ny, Nz, 3) or (Ny, Nz, 3)
    Ey_iso = efield_iso[..., 1]
    Ey_even = efield_even[..., 1]
    Ey_odd = efield_odd[..., 1]

    # Normalize fields
    norm_iso = np.sqrt(np.sum(np.abs(Ey_iso)**2))
    norm_even = np.sqrt(np.sum(np.abs(Ey_even)**2))
    norm_odd = np.sqrt(np.sum(np.abs(Ey_odd)**2))

    c_even = np.abs(np.sum(np.conj(Ey_even) * Ey_iso)) / (norm_iso * norm_even)
    c_odd = np.abs(np.sum(np.conj(Ey_odd) * Ey_iso)) / (norm_iso * norm_odd)

    p_even = c_even**2
    p_odd = c_odd**2
    tot_p = p_even + p_odd
    p_even_rel = p_even / tot_p
    p_odd_rel = p_odd / tot_p

    max_p2 = 4.0 * p_even_rel * p_odd_rel
    min_p1 = (np.sqrt(p_even_rel) - np.sqrt(p_odd_rel))**2

    print(f"Gap = {gap*1000:3.0f} nm: c_even^2={p_even_rel:.4f}, c_odd^2={p_odd_rel:.4f} -> Max P2={max_p2*100:.2f}%, Min P1={min_p1*100:.2f}%")
    return max_p2, min_p1

if __name__ == "__main__":
    print("Testing supermode excitation asymmetry vs gap (w=240nm, tpcm=150nm, am):")
    for g in [0.060, 0.080, 0.100, 0.120, 0.150, 0.180, 0.200]:
        test_overlap(gap=g)
