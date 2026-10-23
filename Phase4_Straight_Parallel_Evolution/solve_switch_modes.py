"""
PROJECT JANUS: 3D VECTOR MAXWELL MODE SOLVER (MPB)
Finds exact modal effective indices, coupling coefficients kappa,
and confinement factors for various Sb2S3 switch cross-sections.
"""
import math
import numpy as np
import meep as mp
from meep import mpb

def solve_coupler_modes(w_wg=0.250, h_wg=0.220, gap=0.045, t_pcm=0.020, 
                        pcm_loc="top", pcm_state="amorphous", resolution=32):
    """
    Computes fundamental even and odd modes for a coupled waveguide system.
    Returns: n_even, n_odd, kappa, L_c (pi / (2 * kappa))
    """
    lambda_0 = 1.064
    k0 = 2 * math.pi / lambda_0
    fcen = 1.0 / lambda_0

    n_si = 3.480
    n_sio2 = 1.449
    n_pcm = 2.70 if pcm_state == "amorphous" else 3.30

    clad_mat = mp.Medium(index=n_sio2)
    core_mat = mp.Medium(index=n_si)
    pcm_mat = mp.Medium(index=n_pcm)

    y_sep = w_wg + gap
    y1 = y_sep / 2.0
    y2 = -y_sep / 2.0

    sy = 2.5
    sz = 2.0
    cell = mp.Vector3(0, sy, sz)

    geometry = [
        # Cladding background is handled by default_material
        # Lower SiO2 substrate
        mp.Block(center=mp.Vector3(0, 0, -sz/4.0), size=mp.Vector3(mp.inf, mp.inf, sz/2.0), material=clad_mat),
        # Waveguide 1
        mp.Block(center=mp.Vector3(0, y1, h_wg/2.0), size=mp.Vector3(mp.inf, w_wg, h_wg), material=core_mat),
        # Waveguide 2
        mp.Block(center=mp.Vector3(0, y2, h_wg/2.0), size=mp.Vector3(mp.inf, w_wg, h_wg), material=core_mat),
    ]

    if pcm_loc == "top" and t_pcm > 0:
        # PCM on WG1
        geometry.append(mp.Block(
            center=mp.Vector3(0, y1, h_wg + t_pcm/2.0),
            size=mp.Vector3(mp.inf, w_wg, t_pcm),
            material=pcm_mat
        ))
    elif pcm_loc == "slot" and gap > 0:
        # PCM in the gap
        geometry.append(mp.Block(
            center=mp.Vector3(0, 0, h_wg/2.0),
            size=mp.Vector3(mp.inf, gap, h_wg),
            material=pcm_mat
        ))

    # Guess effective index around 2.45
    k_guess = 2.45 * fcen
    k_points = [mp.Vector3(k_guess, 0, 0)]

    ms = mpb.ModeSolver(
        geometry=geometry,
        geometry_lattice=mp.Lattice(size=cell),
        k_points=k_points,
        resolution=resolution,
        default_material=clad_mat,
        num_bands=4
    )
    # Find k for specified frequency fcen
    k_res = ms.find_k(mp.NO_PARITY, fcen, 1, 2, mp.Vector3(1, 0, 0), 1e-4, k_guess, 0.5*fcen, 3.48*fcen)
    
    n_eff_1 = k_res[0] / fcen
    n_eff_2 = k_res[1] / fcen

    n_even = max(n_eff_1, n_eff_2)
    n_odd = min(n_eff_1, n_eff_2)

    delta_n = n_even - n_odd
    kappa = math.pi * delta_n / lambda_0
    l_c = math.pi / (2.0 * max(kappa, 1e-6))

    return {
        "pcm_loc": pcm_loc, "pcm_state": pcm_state, "t_pcm": t_pcm, "gap": gap,
        "n_even": round(n_even, 5), "n_odd": round(n_odd, 5), "delta_n": round(delta_n, 5),
        "kappa_um": round(kappa, 4), "l_c_um": round(l_c, 3)
    }

if __name__ == "__main__":
    print("=" * 70)
    print("MPB 3D VECTOR MAXWELL MODE SOLVER FOR SWITCH ARCHITECTURES")
    print("=" * 70)
    
    # 1. Baseline 20 nm top PCM
    print("\n[1] Top-Loaded 20 nm Sb2S3 (gap = 45 nm, w = 250 nm):")
    res_top20_am = solve_coupler_modes(pcm_loc="top", t_pcm=0.020, gap=0.045, pcm_state="amorphous")
    res_top20_cr = solve_coupler_modes(pcm_loc="top", t_pcm=0.020, gap=0.045, pcm_state="crystalline")
    print(f"  Amorphous  : n_even={res_top20_am['n_even']}, n_odd={res_top20_am['n_odd']}, L_c={res_top20_am['l_c_um']} um")
    print(f"  Crystalline: n_even={res_top20_cr['n_even']}, n_odd={res_top20_cr['n_odd']}, L_c={res_top20_cr['l_c_um']} um")

    # 2. Thicker 60 nm top PCM
    print("\n[2] Top-Loaded 60 nm Sb2S3 (gap = 45 nm, w = 250 nm):")
    res_top60_am = solve_coupler_modes(pcm_loc="top", t_pcm=0.060, gap=0.045, pcm_state="amorphous")
    res_top60_cr = solve_coupler_modes(pcm_loc="top", t_pcm=0.060, gap=0.045, pcm_state="crystalline")
    print(f"  Amorphous  : n_even={res_top60_am['n_even']}, n_odd={res_top60_am['n_odd']}, L_c={res_top60_am['l_c_um']} um")
    print(f"  Crystalline: n_even={res_top60_cr['n_even']}, n_odd={res_top60_cr['n_odd']}, L_c={res_top60_cr['l_c_um']} um")

    # 3. Slot-loaded PCM in gap
    print("\n[3] Slot-Loaded Sb2S3 in Gap (gap = 45 nm, w = 250 nm):")
    res_slot_am = solve_coupler_modes(pcm_loc="slot", gap=0.045, pcm_state="amorphous")
    res_slot_cr = solve_coupler_modes(pcm_loc="slot", gap=0.045, pcm_state="crystalline")
    print(f"  Amorphous  : n_even={res_slot_am['n_even']}, n_odd={res_slot_am['n_odd']}, L_c={res_slot_am['l_c_um']} um")
    print(f"  Crystalline: n_even={res_slot_cr['n_even']}, n_odd={res_slot_cr['n_odd']}, L_c={res_slot_cr['l_c_um']} um")
