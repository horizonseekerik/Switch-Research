"""
PROJECT JANUS: GATE 1 ISOLATED C^2 HERMITE S-BEND COLD TRANSMISSION BENCHMARK
=============================================================================
Physics Principle:
  Evaluates 3D Maxwell FDTD transmission through the isolated bare-silicon C^2 
  Hermite S-bend to confirm that bend radiation loss is strictly < 0.1500000000%.
  
  C^2 polynomial formulation:
    y(u) = y_start + dy * (10*u^3 - 15*u^4 + 6*u^5)
    Guarantees d^2y/dx^2 = 0 at both endpoints (zero curvature transition).

Reports all values to >= 10 decimal places.
"""

import sys, os, math, time
import numpy as np
import meep as mp

def run_gate1_sbend(l_bend=1.200, dy=0.108, resolution=20):
    lambda_0 = 1.064
    fcen = 1.0 / lambda_0
    df = 0.10 * fcen

    w_wg = 0.240
    h_wg = 0.220
    core_index = 3.480
    clad_index = 1.449

    clad_mat = mp.Medium(index=clad_index)
    core_mat = mp.Medium(index=core_index)

    l_in = 0.500
    l_out = 0.600
    cell_len = l_in + l_bend + l_out
    cell_w = dy + w_wg + 1.200
    cell_h = h_wg + 1.000

    dpml = 0.500
    pad = 0.350

    sx = cell_len + 2 * (dpml + pad)
    sy = cell_w + 2 * (dpml + pad)
    sz = cell_h + 2 * (dpml + pad)

    cell = mp.Vector3(sx, sy, sz)
    pml_layers = [mp.PML(dpml)]

    x_left = -sx / 2.0
    x0 = -cell_len / 2.0
    x_in_end = x0 + l_in
    x_bend_end = x_in_end + l_bend
    x_mon = x_bend_end + 0.300

    mon_size = mp.Vector3(0, 0.600, 0.600)

    # C^2 Hermite polynomial vertices
    n_pts = 45
    u = np.linspace(0, 1, n_pts)
    x_arr = x_in_end + u * l_bend
    y_arr = -dy * (10.0 * u**3 - 15.0 * u**4 + 6.0 * u**5)

    bot = [mp.Vector3(xi, yi - w_wg / 2.0, 0) for xi, yi in zip(x_arr, y_arr)]
    top = [mp.Vector3(xi, yi + w_wg / 2.0, 0) for xi, yi in reversed(list(zip(x_arr, y_arr)))]
    poly = bot + top

    cx = float(np.mean([pt.x for pt in poly]))
    cy = float(np.mean([pt.y for pt in poly]))
    cz = float(h_wg / 2.0)
    bend_prism = mp.Prism(vertices=poly, height=h_wg, axis=mp.Vector3(0, 0, 1),
                          center=mp.Vector3(cx, cy, cz), material=core_mat)

    src_pos_x = x_left + dpml + 0.100
    sources = [
        mp.EigenModeSource(
            src=mp.GaussianSource(fcen, fwidth=df),
            center=mp.Vector3(src_pos_x, 0, h_wg / 2.0),
            size=mp.Vector3(0, w_wg * 3.0, h_wg * 3.0),
            eig_band=1,
            direction=mp.X
        )
    ]

    # --- 1. Straight Reference Waveguide ---
    geom_ref = [
        mp.Block(center=mp.Vector3(0, 0, -sz / 4.0), size=mp.Vector3(mp.inf, mp.inf, sz / 2.0), material=clad_mat),
        mp.Block(center=mp.Vector3(0, 0, h_wg / 2.0), size=mp.Vector3(mp.inf, w_wg, h_wg), material=core_mat)
    ]
    sim_ref = mp.Simulation(cell_size=cell, boundary_layers=pml_layers, geometry=geom_ref,
                            sources=sources, resolution=resolution, default_material=clad_mat)
    f_ref = sim_ref.add_flux(fcen, 0, 1, mp.FluxRegion(center=mp.Vector3(x_mon, 0, h_wg / 2.0), size=mon_size))
    sim_ref.run(until_after_sources=mp.stop_when_fields_decayed(20, mp.Hz, mp.Vector3(x_mon, 0, h_wg / 2.0), 1e-4))
    p_ref = max(mp.get_fluxes(f_ref)[0], 1e-15)

    # Modal decomposition for reference
    mode_ref = sim_ref.get_eigenmode_coefficients(f_ref, [1])
    p_mode_ref = abs(mode_ref.alpha[0, 0, 0])**2

    # --- 2. C^2 Hermite S-Bend Waveguide ---
    geom_bend = [
        mp.Block(center=mp.Vector3(0, 0, -sz / 4.0), size=mp.Vector3(mp.inf, mp.inf, sz / 2.0), material=clad_mat),
        mp.Block(center=mp.Vector3((x_left + x_in_end) / 2.0, 0, h_wg / 2.0), size=mp.Vector3(x_in_end - x_left, w_wg, h_wg), material=core_mat),
        bend_prism,
        mp.Block(center=mp.Vector3((x_bend_end + sx / 2.0) / 2.0, -dy, h_wg / 2.0), size=mp.Vector3(sx / 2.0 - x_bend_end, w_wg, h_wg), material=core_mat)
    ]
    sim_bend = mp.Simulation(cell_size=cell, boundary_layers=pml_layers, geometry=geom_bend,
                             sources=sources, resolution=resolution, default_material=clad_mat)
    f_bend = sim_bend.add_flux(fcen, 0, 1, mp.FluxRegion(center=mp.Vector3(x_mon, -dy, h_wg / 2.0), size=mon_size))
    sim_bend.run(until_after_sources=mp.stop_when_fields_decayed(20, mp.Hz, mp.Vector3(x_mon, -dy, h_wg / 2.0), 1e-4))
    p_bend = mp.get_fluxes(f_bend)[0]

    # Modal decomposition for bend
    mode_bend = sim_bend.get_eigenmode_coefficients(f_bend, [1])
    p_mode_bend = abs(mode_bend.alpha[0, 0, 0])**2

    flux_trans = (p_bend / p_ref) * 100.0
    flux_loss = 100.0 - flux_trans
    modal_trans = (p_mode_bend / p_mode_ref) * 100.0
    modal_loss = 100.0 - modal_trans
    il_db = -10.0 * math.log10(max(modal_trans / 100.0, 1e-12))

    print("=" * 95)
    print("GATE 1: C^2 HERMITE S-BEND 3D MAXWELL FDTD BENCHMARK RESULTS")
    print(f"Parameters: L_bend = {l_bend:.4f} um, dy = {dy*1e3:.2f} nm, Resolution = {resolution}")
    print("=" * 95)
    print(f"  Modal Transmission (TE0) = {modal_trans:.10f}%")
    print(f"  Modal Radiation Loss    = {modal_loss:.10f}%")
    print(f"  Spatial Flux Trans      = {flux_trans:.10f}%")
    print(f"  Spatial Flux Loss       = {flux_loss:.10f}%")
    print(f"  Insertion Loss (IL)     = {il_db:.10f} dB")
    print("=" * 95)
    if modal_loss < 0.1500000000:
        print(">> GATE 1 PASS: Modal bend loss strictly < 0.1500000000%!")
    else:
        print(f">> GATE 1 WARNING: Modal bend loss = {modal_loss:.10f}%")
    print("=" * 95)

if __name__ == '__main__':
    run_gate1_sbend()
