"""
EXACT 3D MAXWELL FDTD SPATIAL POWER EVOLUTION IN PURE PARALLEL SWITCH CELL
===========================================================================
- Two straight parallel waveguides: WG1 (y=0), WG2 (y=-(w+gap))
- Both covered by Sb2S3 (t_pcm = 160 nm)
- Measures power in WG1 and WG2 along x from 0 to 6.0 um
- State 0: Both amorphous (n = 2.70)
- State 1: WG1 crystalline (n = 3.30), WG2 amorphous (n = 2.70)
- NO BENDS! Pure fundamental spatial coupling physics.
"""

import sys, math, time
import numpy as np
import meep as mp

def run_spatial_evolution(state="amorphous", resolution=20):
    lambda_0 = 1.0640000000
    fcen = 1.0 / lambda_0
    df = 0.10 * fcen

    w_wg = 0.2400000000
    h_wg = 0.2200000000
    gap = 0.1055100000
    t_pcm = 0.1600000000
    l_coupler = 6.0000000000

    y1 = 0.0000000000
    y2 = -(w_wg + gap)
    y_mid = (y1 + y2) / 2.0

    clad_mat = mp.Medium(index=1.4490000000)
    core_mat = mp.Medium(index=3.4800000000)
    pcm_am = mp.Medium(index=2.7000000000)
    pcm_cr = mp.Medium(index=3.3000000000)

    mat_pcm1 = pcm_cr if state.lower() == "crystalline" else pcm_am
    mat_pcm2 = pcm_am

    l_lead = 0.5000000000
    total_len = l_lead + l_coupler + l_lead
    dpml = 0.5000000000
    pad_y = 0.6000000000
    pad_z = 0.6000000000

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

    # Reference run: single straight waveguide with PCM
    geom_ref = [
        mp.Block(center=mp.Vector3(0, 0, -sz / 4.0), size=mp.Vector3(mp.inf, mp.inf, sz / 2.0), material=clad_mat),
        mp.Block(center=mp.Vector3(0, y1, h_wg / 2.0), size=mp.Vector3(mp.inf, w_wg, h_wg), material=core_mat),
        mp.Block(center=mp.Vector3(0, y1, h_wg + t_pcm / 2.0), size=mp.Vector3(mp.inf, w_wg, t_pcm), material=mat_pcm1)
    ]
    sim_ref = mp.Simulation(cell_size=cell, boundary_layers=pml_layers, geometry=geom_ref,
                            sources=sources, resolution=resolution, default_material=clad_mat)
    f_ref_in = sim_ref.add_flux(fcen, 0, 1, mp.FluxRegion(center=mp.Vector3(x_c_start, y1, h_wg / 2.0), size=mp.Vector3(0, 1.2, 1.2)))
    sim_ref.run(until_after_sources=mp.stop_when_fields_decayed(20, mp.Hz, mp.Vector3(x_c_start, y1, h_wg / 2.0), 1e-4))
    p_ref = max(mp.get_fluxes(f_ref_in)[0], 1e-15)

    # Device run: coupled waveguides with PCM
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

    # Place monitor planes along x
    x_sample = np.linspace(0.5, 5.8, 15)
    flux_p1 = []
    flux_p2 = []
    flux_tot = []

    mon_size = mp.Vector3(0, abs(y2) * 1.5, (h_wg + t_pcm) * 3.0)
    mon_wg1_size = mp.Vector3(0, abs(y_mid - y1) * 2.0, (h_wg + t_pcm) * 3.0)
    mon_wg2_size = mp.Vector3(0, abs(y_mid - y2) * 2.0, (h_wg + t_pcm) * 3.0)

    for dist in x_sample:
        xp = x_c_start + dist
        # Total flux
        ft = sim_dev.add_flux(fcen, 0, 1, mp.FluxRegion(center=mp.Vector3(xp, y_mid, (h_wg + t_pcm) / 2.0), size=mon_size))
        # WG1 flux (upper half: y in [y_mid, +inf))
        f1 = sim_dev.add_flux(fcen, 0, 1, mp.FluxRegion(center=mp.Vector3(xp, y1, (h_wg + t_pcm) / 2.0), size=mon_wg1_size))
        # WG2 flux (lower half: y in (-inf, y_mid])
        f2 = sim_dev.add_flux(fcen, 0, 1, mp.FluxRegion(center=mp.Vector3(xp, y2, (h_wg + t_pcm) / 2.0), size=mon_wg2_size))
        flux_tot.append((dist, ft))
        flux_p1.append((dist, f1))
        flux_p2.append((dist, f2))

    sim_dev.run(until_after_sources=mp.stop_when_fields_decayed(20, mp.Hz, mp.Vector3(x_c_end, y_mid, h_wg / 2.0), 1e-4))

    print(f"\n{'='*95}")
    print(f"[{state.upper()} STATE] PURE PARALLEL COUPLER SPATIAL EVOLUTION (Resolution = {resolution})")
    print(f"{'Dist x (um)':<14} {'WG1 P1 Flux %':<18} {'WG2 P2 Flux %':<18} {'Total Guided %':<18} {'Rad Loss %':<16}")
    print("-" * 95)

    res_data = []
    for (d, ft), (_, f1), (_, f2) in zip(flux_tot, flux_p1, flux_p2):
        p_tot = (mp.get_fluxes(ft)[0] / p_ref) * 100.0
        p_1 = (mp.get_fluxes(f1)[0] / p_ref) * 100.0
        p_2 = (mp.get_fluxes(f2)[0] / p_ref) * 100.0
        rad_loss = max(100.0 - p_tot, 0.0)
        print(f"{d:<14.4f} {p_1:<18.10f} {p_2:<18.10f} {p_tot:<18.10f} {rad_loss:<16.10f}")
        res_data.append((d, p_1, p_2, p_tot, rad_loss))
    print(f"{'='*95}\n")
    sys.stdout.flush()
    return res_data

if __name__ == "__main__":
    t0 = time.time()
    print("=" * 95)
    print("3D FDTD RIGOROUS SPATIAL POWER TRACKING IN PURE PARALLEL SWITCH")
    print("=" * 95)
    d0 = run_spatial_evolution(state="amorphous")
    d1 = run_spatial_evolution(state="crystalline")
    print(f"Total scan runtime = {time.time() - t0:.2f} s")
