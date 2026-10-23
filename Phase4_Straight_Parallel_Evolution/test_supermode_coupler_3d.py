"""
SUPERMODE RIGOROUS POWER EXTRACTION IN PURE PARALLEL SWITCH CELL
================================================================
Measures exact guided power in Port 1 and Port 2 via even/odd supermode projection:
  A1 = (alpha_even + alpha_odd) / sqrt(2)
  A2 = (alpha_even - alpha_odd) / sqrt(2)
  P1 = |A1|^2 / P_ref
  P2 = |A2|^2 / P_ref
  P_tot = P1 + P2
"""

import sys, math, time
import numpy as np
import meep as mp

def run_supermode_analysis(state="amorphous", l_coupler=6.0000000000, resolution=20):
    lambda_0 = 1.0640000000
    fcen = 1.0 / lambda_0
    df = 0.10 * fcen

    w_wg = 0.2400000000
    h_wg = 0.2200000000
    gap = 0.1055100000
    t_pcm = 0.1600000000

    y1 = 0.0000000000
    y2 = -(w_wg + gap)
    y_mid = (y1 + y2) / 2.0

    clad_mat = mp.Medium(index=1.4490000000)
    core_mat = mp.Medium(index=3.4800000000)
    pcm_am = mp.Medium(index=2.7000000000)
    pcm_cr = mp.Medium(index=3.3000000000)

    mat_pcm1 = pcm_cr if state.lower() == "crystalline" else pcm_am
    mat_pcm2 = pcm_am

    l_lead = 0.4000000000
    total_len = l_lead + l_coupler + l_lead
    dpml = 0.5000000000
    pad_y = 0.5000000000
    pad_z = 0.5000000000

    sx = total_len + 2 * dpml
    sy = abs(y2) + w_wg + 2 * (dpml + pad_y)
    sz = h_wg + t_pcm + 2 * (dpml + pad_z)

    cell = mp.Vector3(sx, sy, sz)
    pml_layers = [mp.PML(dpml)]

    x_left = -sx / 2.0
    x_c_start = x_left + dpml + l_lead
    x_c_end = x_c_start + l_coupler
    src_pos_x = x_left + dpml + 0.1500000000

    sources = [
        mp.EigenModeSource(
            src=mp.GaussianSource(fcen, fwidth=df),
            center=mp.Vector3(src_pos_x, y1, h_wg / 2.0),
            size=mp.Vector3(0, w_wg * 3.0, (h_wg + t_pcm) * 3.0),
            eig_band=1,
            direction=mp.X
        )
    ]

    # Reference run: isolated straight waveguide with PCM
    geom_ref = [
        mp.Block(center=mp.Vector3(0, 0, -sz / 4.0), size=mp.Vector3(mp.inf, mp.inf, sz / 2.0), material=clad_mat),
        mp.Block(center=mp.Vector3(0, y1, h_wg / 2.0), size=mp.Vector3(mp.inf, w_wg, h_wg), material=core_mat),
        mp.Block(center=mp.Vector3(0, y1, h_wg + t_pcm / 2.0), size=mp.Vector3(mp.inf, w_wg, t_pcm), material=mat_pcm1)
    ]
    sim_ref = mp.Simulation(cell_size=cell, boundary_layers=pml_layers, geometry=geom_ref,
                            sources=sources, resolution=resolution, default_material=clad_mat)
    f_ref = sim_ref.add_flux(fcen, 0, 1, mp.FluxRegion(center=mp.Vector3(x_c_start, y1, h_wg / 2.0), size=mp.Vector3(0, 1.2, 1.2)))
    sim_ref.run(until_after_sources=mp.stop_when_fields_decayed(20, mp.Hz, mp.Vector3(x_c_start, y1, h_wg / 2.0), 1e-4))
    res_ref = sim_ref.get_eigenmode_coefficients(f_ref, [1], eig_parity=mp.NO_PARITY)
    p_ref_mode = abs(res_ref.alpha[0, 0, 0])**2

    # Device run
    geom_dev = [
        mp.Block(center=mp.Vector3(0, 0, -sz / 4.0), size=mp.Vector3(mp.inf, mp.inf, sz / 2.0), material=clad_mat),
        # WG1 core and PCM
        mp.Block(center=mp.Vector3(0, y1, h_wg / 2.0), size=mp.Vector3(mp.inf, w_wg, h_wg), material=core_mat),
        mp.Block(center=mp.Vector3((x_c_start + x_c_end) / 2.0, y1, h_wg + t_pcm / 2.0),
                 size=mp.Vector3(l_coupler, w_wg, t_pcm), material=mat_pcm1),
        # WG2 core and PCM
        mp.Block(center=mp.Vector3(0, y2, h_wg / 2.0), size=mp.Vector3(mp.inf, w_wg, h_wg), material=core_mat),
        mp.Block(center=mp.Vector3((x_c_start + x_c_end) / 2.0, y2, h_wg + t_pcm / 2.0),
                 size=mp.Vector3(l_coupler, w_wg, t_pcm), material=mat_pcm2),
    ]

    sim_dev = mp.Simulation(cell_size=cell, boundary_layers=pml_layers, geometry=geom_dev,
                            sources=sources, resolution=resolution, default_material=clad_mat)

    # Monitor plane across BOTH waveguides to decompose into supermodes
    mon_w_both = abs(y2) + w_wg + 0.6000000000
    mon_h_both = (h_wg + t_pcm) * 3.0

    x_points = np.arange(5.0, 6.61, 0.10)
    flux_monitors = []
    for xp in x_points:
        f = sim_dev.add_flux(fcen, 0, 1, mp.FluxRegion(center=mp.Vector3(x_c_start + xp, y_mid, (h_wg + t_pcm) / 2.0),
                                                      size=mp.Vector3(0, mon_w_both, mon_h_both)))
        flux_monitors.append((xp, f))

    sim_dev.run(until_after_sources=mp.stop_when_fields_decayed(20, mp.Hz, mp.Vector3(x_c_end, y_mid, h_wg / 2.0), 1e-4))

    print(f"\n{'='*105}")
    print(f"[{state.upper()} STATE] SUPERMODE DECOMPOSITION (Resolution = {resolution})")
    print(f"{'Length x (um)':<14} {'P1 (WG1) %':<18} {'P2 (WG2) %':<18} {'Total Guided %':<18} {'P2 / P_tot %':<16}")
    print("-" * 105)

    results = []
    for xp, f in flux_monitors:
        res = sim_dev.get_eigenmode_coefficients(f, [1, 2], eig_parity=mp.NO_PARITY)
        a_even = res.alpha[0, 0, 0]  # Band 1 (forward)
        a_odd = res.alpha[0, 0, 1]   # Band 2 (forward)

        # In coupled waveguides, WG1 and WG2 are orthogonal symmetric/antisymmetric superpositions
        # A1 = (a_even + a_odd) / sqrt(2)
        # A2 = (a_even - a_odd) / sqrt(2)
        p1 = abs(a_even + a_odd)**2 / (2.0 * p_ref_mode) * 100.0
        p2 = abs(a_even - a_odd)**2 / (2.0 * p_ref_mode) * 100.0
        tot = p1 + p2
        frac2 = (p2 / tot) * 100.0 if tot > 0 else 0.0
        print(f"{xp:<14.2f} {p1:<18.10f} {p2:<18.10f} {tot:<18.10f} {frac2:<16.10f}")
        results.append((xp, p1, p2, tot, frac2))
    print(f"{'='*105}\n")
    sys.stdout.flush()
    return results

if __name__ == "__main__":
    t0 = time.time()
    r0 = run_supermode_analysis(state="amorphous")
    r1 = run_supermode_analysis(state="crystalline")
    print(f"Total supermode analysis runtime = {time.time() - t0:.2f} s")
