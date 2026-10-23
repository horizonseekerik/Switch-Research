"""
DIAGNOSTIC: TRACK FLUX ALONG THE SWITCH TO LOCATE EXACT LOSS MECHANISM
======================================================================
Places flux planes at:
- x1: Input waveguide (bare Si)
- x2: End of input taper / Start of coupler
- x3: End of coupler / Start of S-bend
- x4: End of S-bend
- x5: Output monitors
"""

import sys, math
import numpy as np
import meep as mp

def make_ccw_sbend(x_start, x_end, y_start, y_end, w1, w2, n_pts=30):
    dx = x_end - x_start
    dy = y_end - y_start
    u = np.linspace(0, 1, n_pts)
    x = x_start + u * dx
    y = y_start + dy * (3 * u**2 - 2 * u**3)
    w = w1 + (w2 - w1) * (3 * u**2 - 2 * u**3)
    bot = [mp.Vector3(xi, yi - wi / 2.0, 0) for xi, yi, wi in zip(x, y, w)]
    top = [mp.Vector3(xi, yi + wi / 2.0, 0) for xi, yi, wi in reversed(list(zip(x, y, w)))]
    return bot + top

def make_ccw_linear_taper(x_start, x_end, y, w_start, w_end, n_pts=25):
    u = np.linspace(0, 1, n_pts)
    x = x_start + u * (x_end - x_start)
    w = w_start + u * (w_end - w_start)
    bot = [mp.Vector3(xi, y - wi / 2.0, 0) for xi, wi in zip(x, w)]
    top = [mp.Vector3(xi, y + wi / 2.0, 0) for xi, wi in reversed(list(zip(x, w)))]
    return bot + top

def create_prism_exact(vertices, height, z_base, material):
    cx = float(np.mean([pt.x for pt in vertices]))
    cy = float(np.mean([pt.y for pt in vertices]))
    cz = float(z_base + height / 2.0)
    return mp.Prism(vertices=vertices, height=height, axis=mp.Vector3(0, 0, 1),
                    center=mp.Vector3(cx, cy, cz), material=material)

def diagnose_loss(state="amorphous", l_c=4.65, gap=0.10551, t_pcm=0.160, l_tip=0.800, l_bend=0.850, y_pitch=0.450):
    lambda_0 = 1.064
    fcen = 1.0 / lambda_0
    df = 0.10 * fcen

    w_wg = 0.240
    h_wg = 0.220
    core_index = 3.480

    clad_mat = mp.Medium(index=1.449)
    core_mat = mp.Medium(index=core_index)
    pcm_am = mp.Medium(index=2.70)
    pcm_cr = mp.Medium(index=3.30)

    l_in = 0.30
    l_out = 0.30
    iso_buf = 0.08

    y1 = 0.0
    y2_c = -(w_wg + gap)
    y_out2 = -y_pitch

    cell_len = l_in + l_tip + l_c + l_bend + l_tip + l_out
    cell_w = y_pitch + w_wg + 2 * iso_buf
    dpml = 0.50
    pad = 0.35
    sx = cell_len + 2 * (dpml + pad)
    sy = cell_w + 2 * (dpml + pad)
    sz = h_wg + t_pcm + 2 * (dpml + pad)

    cell = mp.Vector3(sx, sy, sz)
    pml_layers = [mp.PML(dpml)]

    x_left = -sx / 2.0
    x0 = -cell_len / 2.0

    x_in_end = x0 + l_in
    x_c_start = x_in_end + l_tip
    x_c_end = x_c_start + l_c
    x_bend_end = x_c_end + l_bend
    x_out_tip_end = x_bend_end + l_tip
    x_mon = x_out_tip_end + 0.15

    src_pos_x = x_left + dpml + 0.10
    sources = [
        mp.EigenModeSource(
            src=mp.GaussianSource(fcen, fwidth=df),
            center=mp.Vector3(src_pos_x, y1, h_wg / 2.0),
            size=mp.Vector3(0, w_wg * 3.0, (h_wg + t_pcm) * 3.0),
            eig_band=1,
            direction=mp.X
        )
    ]

    # Reference
    geom_ref = [
        mp.Block(center=mp.Vector3(0, 0, -sz / 4.0), size=mp.Vector3(mp.inf, mp.inf, sz / 2.0), material=clad_mat),
        mp.Block(center=mp.Vector3(0, y1, h_wg / 2.0), size=mp.Vector3(mp.inf, w_wg, h_wg), material=core_mat)
    ]
    sim_ref = mp.Simulation(cell_size=cell, boundary_layers=pml_layers, geometry=geom_ref,
                            sources=sources, resolution=15, default_material=clad_mat)
    f_ref = sim_ref.add_flux(fcen, 0, 1, mp.FluxRegion(center=mp.Vector3(src_pos_x + 0.2, 0, h_wg / 2.0), size=mp.Vector3(0, 1.0, 1.0)))
    sim_ref.run(until_after_sources=mp.stop_when_fields_decayed(20, mp.Hz, mp.Vector3(src_pos_x + 0.2, 0, h_wg / 2.0), 1e-4))
    p_ref = max(mp.get_fluxes(f_ref)[0], 1e-15)

    # Device
    geom_dev = [
        mp.Block(center=mp.Vector3(0, 0, -sz / 4.0), size=mp.Vector3(mp.inf, mp.inf, sz / 2.0), material=clad_mat),
        mp.Block(center=mp.Vector3(0, y1, h_wg / 2.0), size=mp.Vector3(mp.inf, w_wg, h_wg), material=core_mat),
    ]

    poly_wg2_in = make_ccw_linear_taper(x_in_end, x_c_start, y2_c, 0.02, w_wg)
    geom_dev.append(create_prism_exact(poly_wg2_in, h_wg, 0.0, core_mat))
    geom_dev.append(mp.Block(center=mp.Vector3((x_c_start + x_c_end) / 2.0, y2_c, h_wg / 2.0),
                             size=mp.Vector3(l_c, w_wg, h_wg), material=core_mat))
    poly_wg2_bend = make_ccw_sbend(x_c_end, x_bend_end, y2_c, y_out2, w_wg, w_wg)
    geom_dev.append(create_prism_exact(poly_wg2_bend, h_wg, 0.0, core_mat))
    geom_dev.append(mp.Block(center=mp.Vector3((x_bend_end + sx/2.0) / 2.0, y_out2, h_wg / 2.0),
                             size=mp.Vector3(sx/2.0 - x_bend_end, w_wg, h_wg), material=core_mat))

    mat_pcm1 = pcm_cr if state == "crystalline" else pcm_am
    mat_pcm2 = pcm_am

    poly_pcm1_in = make_ccw_linear_taper(x_in_end, x_c_start, y1, 0.02, w_wg)
    geom_dev.append(create_prism_exact(poly_pcm1_in, t_pcm, h_wg, mat_pcm1))
    geom_dev.append(mp.Block(center=mp.Vector3((x_c_start + x_c_end) / 2.0, y1, h_wg + t_pcm / 2.0),
                             size=mp.Vector3(l_c, w_wg, t_pcm), material=mat_pcm1))
    geom_dev.append(mp.Block(center=mp.Vector3((x_c_end + x_bend_end) / 2.0, y1, h_wg + t_pcm / 2.0),
                             size=mp.Vector3(l_bend, w_wg, t_pcm), material=mat_pcm1))
    poly_pcm1_out = make_ccw_linear_taper(x_bend_end, x_out_tip_end, y1, w_wg, 0.02)
    geom_dev.append(create_prism_exact(poly_pcm1_out, t_pcm, h_wg, mat_pcm1))

    poly_pcm2_in = make_ccw_linear_taper(x_in_end, x_c_start, y2_c, 0.02, w_wg)
    geom_dev.append(create_prism_exact(poly_pcm2_in, t_pcm, h_wg, mat_pcm2))
    geom_dev.append(mp.Block(center=mp.Vector3((x_c_start + x_c_end) / 2.0, y2_c, h_wg + t_pcm / 2.0),
                             size=mp.Vector3(l_c, w_wg, t_pcm), material=mat_pcm2))
    poly_pcm2_sbend = make_ccw_sbend(x_c_end, x_bend_end, y2_c, y_out2, w_wg, 0.02)
    geom_dev.append(create_prism_exact(poly_pcm2_sbend, t_pcm, h_wg, mat_pcm2))

    sim_dev = mp.Simulation(cell_size=cell, boundary_layers=pml_layers, geometry=geom_dev,
                            sources=sources, resolution=15, default_material=clad_mat)

    # Monitor planes along the entire device
    planes = [
        ("Input bare Si", x_in_end),
        ("Coupler Start", x_c_start),
        ("Coupler Mid  ", (x_c_start + x_c_end) / 2.0),
        ("Coupler End  ", x_c_end),
        ("S-Bend End   ", x_bend_end),
        ("Output Lead  ", x_mon)
    ]
    f_planes = []
    for name, xp in planes:
        fp = sim_dev.add_flux(fcen, 0, 1, mp.FluxRegion(center=mp.Vector3(xp, -y_pitch/2.0, h_wg/2.0), size=mp.Vector3(0, sy*0.8, sz*0.8)))
        f_planes.append((name, fp))

    sim_dev.run(until_after_sources=mp.stop_when_fields_decayed(20, mp.Hz, mp.Vector3(x_mon, 0, h_wg/2.0), 1e-4))

    print(f"\n{'='*70}")
    print(f"[{state.upper()}] SPATIAL FLUX TRACKING (Total Input Power = 100.0%):")
    print(f"{'Plane Location':<22} {'Cross-Section x':<16} {'Total Flux %':<16} {'Loss Delta %':<16}")
    print("-" * 70)
    prev = 100.0
    for name, fp in f_planes:
        flx = (mp.get_fluxes(fp)[0] / p_ref) * 100.0
        delta = prev - flx
        print(f"{name:<22} {flx:<16.2f}% {delta:<16.2f}%")
        prev = flx
    print(f"{'='*70}\n")

if __name__ == "__main__":
    diagnose_loss(state="amorphous")
    diagnose_loss(state="crystalline")
