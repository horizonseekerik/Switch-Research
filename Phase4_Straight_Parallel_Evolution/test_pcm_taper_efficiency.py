"""
PROJECT JANUS: PCM TAPER ADIABATIC TRANSMISSION BENCHMARK (3D FDTD)
===================================================================
Tests pure TE0 mode transmission through:
  Bare Si -> PCM Taper (l_tip) -> Uniform PCM Body (1.0 um) -> Exit Taper (l_tip) -> Bare Si
Sweeps l_tip = [0.25, 0.40, 0.60, 0.80] um to find the adiabatic threshold for > 99% modal power.
"""

import sys, os, math, time
import numpy as np
import meep as mp

def test_taper_transmission(l_tip=0.250, t_pcm=0.190, pcm_state="amorphous", resolution=15):
    lambda_0 = 1.064
    fcen = 1.0 / lambda_0
    df = 0.10 * fcen

    w_wg = 0.240
    h_wg = 0.220
    core_index = 3.480
    clad_mat = mp.Medium(index=1.449)
    core_mat = mp.Medium(index=core_index)

    n_pcm = 3.30 if pcm_state == "crystalline" else 2.70
    cond = 2 * math.pi * fcen * (8.0e-5 if pcm_state == "crystalline" else 1.0e-5)
    pcm_mat = mp.Medium(index=n_pcm, D_conductivity=cond)

    l_body = 1.000
    l_lead = 0.600
    cell_len = 2 * l_lead + 2 * l_tip + l_body

    dpml = 0.50
    pad = 0.35
    sx = cell_len + 2 * (dpml + pad)
    sy = w_wg + 2 * (dpml + pad)
    sz = h_wg + t_pcm + 2 * (dpml + pad)

    cell = mp.Vector3(sx, sy, sz)
    pml_layers = [mp.PML(dpml)]

    x_left = -sx / 2.0
    x0 = -cell_len / 2.0
    x_taper_in = x0 + l_lead
    x_body_start = x_taper_in + l_tip
    x_body_end = x_body_start + l_body
    x_taper_out = x_body_end + l_tip
    x_mon = x_taper_out + 0.30

    sources = [
        mp.EigenModeSource(
            src=mp.GaussianSource(fcen, fwidth=df),
            center=mp.Vector3(x_left + dpml + 0.10, 0, h_wg / 2.0),
            size=mp.Vector3(0, w_wg * 3.0, (h_wg + t_pcm) * 3.0),
            eig_band=1,
            direction=mp.X
        )
    ]

    # Reference straight bare-Si waveguide
    geom_ref = [
        mp.Block(center=mp.Vector3(0, 0, -sz / 4.0), size=mp.Vector3(mp.inf, mp.inf, sz / 2.0), material=clad_mat),
        mp.Block(center=mp.Vector3(0, 0, h_wg / 2.0), size=mp.Vector3(mp.inf, w_wg, h_wg), material=core_mat)
    ]
    sim_ref = mp.Simulation(cell_size=cell, boundary_layers=pml_layers, geometry=geom_ref,
                            sources=sources, resolution=resolution, default_material=clad_mat)
    f_ref = sim_ref.add_flux(fcen, 0, 1, mp.FluxRegion(center=mp.Vector3(x_mon, 0, h_wg / 2.0), size=mp.Vector3(0, 0.40, 0.40)))
    sim_ref.run(until_after_sources=mp.stop_when_fields_decayed(15, mp.Hz, mp.Vector3(x_mon, 0, h_wg / 2.0), 1e-3))
    p_ref_flux = max(mp.get_fluxes(f_ref)[0], 1e-12)
    ref_mode = sim_ref.get_eigenmode_coefficients(f_ref, [1], eig_parity=mp.NO_PARITY)
    p_ref_mode = max(abs(ref_mode.alpha[0, 0, 0])**2, 1e-12)

    # Taper + Body + Taper waveguide
    def make_linear_taper(x_start, x_end, w_start, w_end, n_pts=20):
        u = np.linspace(0, 1, n_pts)
        x = x_start + u * (x_end - x_start)
        w = w_start + u * (w_end - w_start)
        bot = [mp.Vector3(xi, -wi / 2.0, 0) for xi, wi in zip(x, w)]
        top = [mp.Vector3(xi, wi / 2.0, 0) for xi, wi in reversed(list(zip(x, w)))]
        return bot + top

    poly_in = make_linear_taper(x_taper_in, x_body_start, 0.030, w_wg)
    cx_in = float(np.mean([pt.x for pt in poly_in]))
    prism_in = mp.Prism(vertices=poly_in, height=t_pcm, axis=mp.Vector3(0, 0, 1),
                        center=mp.Vector3(cx_in, 0, h_wg + t_pcm / 2.0), material=pcm_mat)

    poly_out = make_linear_taper(x_body_end, x_taper_out, w_wg, 0.030)
    cx_out = float(np.mean([pt.x for pt in poly_out]))
    prism_out = mp.Prism(vertices=poly_out, height=t_pcm, axis=mp.Vector3(0, 0, 1),
                         center=mp.Vector3(cx_out, 0, h_wg + t_pcm / 2.0), material=pcm_mat)

    geom_dev = [
        mp.Block(center=mp.Vector3(0, 0, -sz / 4.0), size=mp.Vector3(mp.inf, mp.inf, sz / 2.0), material=clad_mat),
        mp.Block(center=mp.Vector3(0, 0, h_wg / 2.0), size=mp.Vector3(mp.inf, w_wg, h_wg), material=core_mat),
        prism_in,
        mp.Block(center=mp.Vector3((x_body_start + x_body_end) / 2.0, 0, h_wg + t_pcm / 2.0),
                 size=mp.Vector3(l_body, w_wg, t_pcm), material=pcm_mat),
        prism_out
    ]

    sim_dev = mp.Simulation(cell_size=cell, boundary_layers=pml_layers, geometry=geom_dev,
                            sources=sources, resolution=resolution, default_material=clad_mat)
    f_dev = sim_dev.add_flux(fcen, 0, 1, mp.FluxRegion(center=mp.Vector3(x_mon, 0, h_wg / 2.0), size=mp.Vector3(0, 0.40, 0.40)))
    sim_dev.run(until_after_sources=mp.stop_when_fields_decayed(15, mp.Hz, mp.Vector3(x_mon, 0, h_wg / 2.0), 1e-3))
    p_dev_flux = max(mp.get_fluxes(f_dev)[0], 1e-12)
    dev_mode = sim_dev.get_eigenmode_coefficients(f_dev, [1], eig_parity=mp.NO_PARITY)
    p_dev_mode = max(abs(dev_mode.alpha[0, 0, 0])**2, 1e-12)

    flux_trans = (p_dev_flux / p_ref_flux) * 100.0
    mode_trans = (p_dev_mode / p_ref_mode) * 100.0
    flux_loss = 100.0 - flux_trans
    mode_loss = 100.0 - mode_trans

    print(f"[{pcm_state.upper()}] l_tip = {l_tip*1e3:5.1f} nm: Modal Trans = {mode_trans:12.8f}% (Loss = {mode_loss:10.8f}%), Flux Trans = {flux_trans:12.8f}% (Loss = {flux_loss:10.8f}%)")
    sys.stdout.flush()
    return mode_trans, flux_trans

if __name__ == '__main__':
    print("=" * 100)
    print("  PCM TAPER EFFICIENCY vs TAPER LENGTH (3D FDTD)")
    print("=" * 100)
    for state in ["amorphous", "crystalline"]:
        print(f"\n--- State: {state.upper()} ---")
        for tip in [0.250, 0.400, 0.600]:
            test_taper_transmission(l_tip=tip, t_pcm=0.190, pcm_state=state, resolution=15)
    print("=" * 100)
