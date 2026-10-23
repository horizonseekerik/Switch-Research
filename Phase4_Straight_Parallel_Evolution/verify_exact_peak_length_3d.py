"""
PINPOINT EXACT COUPLING LENGTH FOR DUAL >98% IN 3D FDTD
=========================================================
Sweeps L_c in [6.10, 6.30, 6.50] um to find the peak power transfer.
Resolution = 20 pixels/um
"""

import sys, math, time
import numpy as np
import meep as mp

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

def test_length(l_c=6.3000000000, resolution=20):
    lambda_0 = 1.0640000000
    fcen = 1.0 / lambda_0
    df = 0.10 * fcen

    w_wg = 0.2400000000
    h_wg = 0.2200000000
    gap = 0.1030000000
    t_pcm = 0.1600000000
    l_tip = 1.5000000000
    l_body = l_c - 2.0 * l_tip

    y1 = 0.0000000000
    y2 = -(w_wg + gap)

    clad_mat = mp.Medium(index=1.4490000000)
    core_mat = mp.Medium(index=3.4800000000)
    pcm_am = mp.Medium(index=2.7000000000)
    pcm_cr = mp.Medium(index=3.3000000000)

    l_in = 0.3000000000
    l_out = 0.3000000000
    iso_buf = 0.0800000000

    cell_len = l_in + l_c + l_out
    cell_w = abs(y2) + w_wg + 2.0 * iso_buf
    cell_area = cell_len * cell_w

    dpml = 0.5000000000
    pad_y = 0.4000000000
    pad_z = 0.4500000000

    sx = cell_len + 2.0 * dpml
    sy = abs(y2) + w_wg + 2.0 * (dpml + pad_y)
    sz = h_wg + t_pcm + 2.0 * (dpml + pad_z)

    cell = mp.Vector3(sx, sy, sz)
    pml_layers = [mp.PML(dpml)]

    x_left = -sx / 2.0
    x0 = -cell_len / 2.0
    x_c_start = x0 + l_in
    x_body_start = x_c_start + l_tip
    x_body_end = x_body_start + l_body
    x_c_end = x_body_end + l_tip
    x_mon = x_c_end + 0.1500000000

    src_pos_x = x_left + dpml + 0.1200000000
    sources = [
        mp.EigenModeSource(
            src=mp.GaussianSource(fcen, fwidth=df),
            center=mp.Vector3(src_pos_x, y1, h_wg / 2.0),
            size=mp.Vector3(0, w_wg * 3.0, (h_wg + t_pcm) * 3.0),
            eig_band=1,
            direction=mp.X
        )
    ]

    mon_w = min(abs(y2) * 0.90, 0.3200000000)
    mon_h = (h_wg + t_pcm) * 2.5

    # Reference
    geom_ref = [
        mp.Block(center=mp.Vector3(0, 0, -sz / 4.0), size=mp.Vector3(mp.inf, mp.inf, sz / 2.0), material=clad_mat),
        mp.Block(center=mp.Vector3(0, y1, h_wg / 2.0), size=mp.Vector3(mp.inf, w_wg, h_wg), material=core_mat)
    ]
    sim_ref = mp.Simulation(cell_size=cell, boundary_layers=pml_layers, geometry=geom_ref,
                            sources=sources, resolution=resolution, default_material=clad_mat)
    f_ref = sim_ref.add_flux(fcen, 0, 1, mp.FluxRegion(center=mp.Vector3(x_mon, y1, h_wg / 2.0), size=mp.Vector3(0, mon_w, mon_h)))
    sim_ref.run(until_after_sources=mp.stop_when_fields_decayed(20, mp.Hz, mp.Vector3(x_mon, y1, h_wg / 2.0), 1e-4))
    res_ref = sim_ref.get_eigenmode_coefficients(f_ref, [1], eig_parity=mp.NO_PARITY)
    p_ref = max(abs(res_ref.alpha[0, 0, 0])**2, 1e-15)

    # State 0
    geom_s0 = [
        mp.Block(center=mp.Vector3(0, 0, -sz / 4.0), size=mp.Vector3(mp.inf, mp.inf, sz / 2.0), material=clad_mat),
        mp.Block(center=mp.Vector3(0, y1, h_wg / 2.0), size=mp.Vector3(mp.inf, w_wg, h_wg), material=core_mat),
        mp.Block(center=mp.Vector3(0, y2, h_wg / 2.0), size=mp.Vector3(mp.inf, w_wg, h_wg), material=core_mat),
    ]
    poly1_in = make_ccw_linear_taper(x_c_start, x_body_start, y1, 0.015, w_wg)
    geom_s0.append(create_prism_exact(poly1_in, t_pcm, h_wg, pcm_am))
    geom_s0.append(mp.Block(center=mp.Vector3((x_body_start + x_body_end) / 2.0, y1, h_wg + t_pcm / 2.0),
                            size=mp.Vector3(l_body, w_wg, t_pcm), material=pcm_am))
    poly1_out = make_ccw_linear_taper(x_body_end, x_c_end, y1, w_wg, 0.015)
    geom_s0.append(create_prism_exact(poly1_out, t_pcm, h_wg, pcm_am))

    poly2_in = make_ccw_linear_taper(x_c_start, x_body_start, y2, 0.015, w_wg)
    geom_s0.append(create_prism_exact(poly2_in, t_pcm, h_wg, pcm_am))
    geom_s0.append(mp.Block(center=mp.Vector3((x_body_start + x_body_end) / 2.0, y2, h_wg + t_pcm / 2.0),
                            size=mp.Vector3(l_body, w_wg, t_pcm), material=pcm_am))
    poly2_out = make_ccw_linear_taper(x_body_end, x_c_end, y2, w_wg, 0.015)
    geom_s0.append(create_prism_exact(poly2_out, t_pcm, h_wg, pcm_am))

    sim_s0 = mp.Simulation(cell_size=cell, boundary_layers=pml_layers, geometry=geom_s0,
                           sources=sources, resolution=resolution, default_material=clad_mat)
    f0_p1 = sim_s0.add_flux(fcen, 0, 1, mp.FluxRegion(center=mp.Vector3(x_mon, y1, h_wg / 2.0), size=mp.Vector3(0, mon_w, mon_h)))
    f0_p2 = sim_s0.add_flux(fcen, 0, 1, mp.FluxRegion(center=mp.Vector3(x_mon, y2, h_wg / 2.0), size=mp.Vector3(0, mon_w, mon_h)))
    sim_s0.run(until_after_sources=mp.stop_when_fields_decayed(20, mp.Hz, mp.Vector3(x_mon, y2, h_wg / 2.0), 1e-4))
    r0_p1 = sim_s0.get_eigenmode_coefficients(f0_p1, [1], eig_parity=mp.NO_PARITY)
    r0_p2 = sim_s0.get_eigenmode_coefficients(f0_p2, [1], eig_parity=mp.NO_PARITY)
    p0_1 = abs(r0_p1.alpha[0, 0, 0])**2 / p_ref
    p0_2 = abs(r0_p2.alpha[0, 0, 0])**2 / p_ref

    # State 1
    geom_s1 = [
        mp.Block(center=mp.Vector3(0, 0, -sz / 4.0), size=mp.Vector3(mp.inf, mp.inf, sz / 2.0), material=clad_mat),
        mp.Block(center=mp.Vector3(0, y1, h_wg / 2.0), size=mp.Vector3(mp.inf, w_wg, h_wg), material=core_mat),
        mp.Block(center=mp.Vector3(0, y2, h_wg / 2.0), size=mp.Vector3(mp.inf, w_wg, h_wg), material=core_mat),
    ]
    poly1_in_cr = make_ccw_linear_taper(x_c_start, x_body_start, y1, 0.015, w_wg)
    geom_s1.append(create_prism_exact(poly1_in_cr, t_pcm, h_wg, pcm_cr))
    geom_s1.append(mp.Block(center=mp.Vector3((x_body_start + x_body_end) / 2.0, y1, h_wg + t_pcm / 2.0),
                            size=mp.Vector3(l_body, w_wg, t_pcm), material=pcm_cr))
    poly1_out_cr = make_ccw_linear_taper(x_body_end, x_c_end, y1, w_wg, 0.015)
    geom_s1.append(create_prism_exact(poly1_out_cr, t_pcm, h_wg, pcm_cr))

    poly2_in_am = make_ccw_linear_taper(x_c_start, x_body_start, y2, 0.015, w_wg)
    geom_s1.append(create_prism_exact(poly2_in_am, t_pcm, h_wg, pcm_am))
    geom_s1.append(mp.Block(center=mp.Vector3((x_body_start + x_body_end) / 2.0, y2, h_wg + t_pcm / 2.0),
                            size=mp.Vector3(l_body, w_wg, t_pcm), material=pcm_am))
    poly2_out_am = make_ccw_linear_taper(x_body_end, x_c_end, y2, w_wg, 0.015)
    geom_s1.append(create_prism_exact(poly2_out_am, t_pcm, h_wg, pcm_am))

    sim_s1 = mp.Simulation(cell_size=cell, boundary_layers=pml_layers, geometry=geom_s1,
                           sources=sources, resolution=resolution, default_material=clad_mat)
    f1_p1 = sim_s1.add_flux(fcen, 0, 1, mp.FluxRegion(center=mp.Vector3(x_mon, y1, h_wg / 2.0), size=mp.Vector3(0, mon_w, mon_h)))
    f1_p2 = sim_s1.add_flux(fcen, 0, 1, mp.FluxRegion(center=mp.Vector3(x_mon, y2, h_wg / 2.0), size=mp.Vector3(0, mon_w, mon_h)))
    sim_s1.run(until_after_sources=mp.stop_when_fields_decayed(20, mp.Hz, mp.Vector3(x_mon, y1, h_wg / 2.0), 1e-4))
    r1_p1 = sim_s1.get_eigenmode_coefficients(f1_p1, [1], eig_parity=mp.NO_PARITY)
    r1_p2 = sim_s1.get_eigenmode_coefficients(f1_p2, [1], eig_parity=mp.NO_PARITY)
    p1_1 = abs(r1_p1.alpha[0, 0, 0])**2 / p_ref
    p1_2 = abs(r1_p2.alpha[0, 0, 0])**2 / p_ref

    print(f"L_c = {l_c:.3f} um | S0 P2 = {p0_2*100:10.6f}% (P1 = {p0_1*100:10.6f}%) | S1 P1 = {p1_1*100:10.6f}% (P2 = {p1_2*100:10.6f}%) | Footprint = {cell_area:.3f} um^2")
    sys.stdout.flush()
    return l_c, p0_2*100, p0_1*100, p1_1*100, p1_2*100, cell_area

if __name__ == "__main__":
    t0 = time.time()
    print("=" * 105)
    print("SWEEPING COUPLING LENGTH TO HIT DUAL >98% PERFECTION")
    print("=" * 105)
    for lc in [6.10, 6.30, 6.50, 6.70]:
        test_length(l_c=lc)
    print(f"Total scan runtime = {time.time() - t0:.2f} s")
