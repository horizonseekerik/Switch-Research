import math
import time
import numpy as np
import meep as mp

def make_ccw_hermite_c2_sbend(x_start, x_end, y_start, y_end, w_start, w_end, n_pts=35):
    u = np.linspace(0, 1, n_pts)
    s = 10.0 * u**3 - 15.0 * u**4 + 6.0 * u**5
    x = x_start + u * (x_end - x_start)
    y_center = y_start + s * (y_end - y_start)
    w = w_start + u * (w_end - w_start)
    bot = [mp.Vector3(xi, yi - wi / 2.0, 0) for xi, yi, wi in zip(x, y_center, w)]
    top = [mp.Vector3(xi, yi + wi / 2.0, 0) for xi, yi, wi in reversed(list(zip(x, y_center, w)))]
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

def run_point(l_c, state="amorphous", gap=0.1037165000, resolution=15):
    lambda_0 = 1.064
    fcen = 1.0 / lambda_0
    df = 0.10 * fcen

    w_wg = 0.240
    h_wg = 0.220
    t_pcm = 0.190
    l_bend = 1.000
    l_tip = 0.600
    y_pitch = 0.450
    l_in = 0.350
    l_out = 0.250
    iso_buf = 0.080

    clad_mat = mp.Medium(index=1.449)
    core_mat = mp.Medium(index=3.480)
    pcm_am = mp.Medium(index=2.70, D_conductivity=2 * math.pi * fcen * 1.0e-5)
    pcm_cr = mp.Medium(index=3.30, D_conductivity=2 * math.pi * fcen * 8.0e-5)

    y1 = 0.0
    y2_c = -(w_wg + gap)
    y_out2 = -y_pitch

    cell_len = l_in + l_c + l_bend + l_tip + l_out
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
    x_c_start = x_in_end
    x_c_end = x_c_start + l_c
    x_bend_end = x_c_end + l_bend
    x_taper_end = x_bend_end + l_tip
    x_mon = x_taper_end + 0.12

    mon_w = min(y_pitch * 0.85, 0.40)
    mon_h = (h_wg + t_pcm) * 2.2

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
                            sources=sources, resolution=resolution, default_material=clad_mat)
    f_ref = sim_ref.add_flux(fcen, 0, 1, mp.FluxRegion(center=mp.Vector3(x_mon, y1, h_wg / 2.0), size=mp.Vector3(0, mon_w, mon_h)))
    sim_ref.run(until_after_sources=mp.stop_when_fields_decayed(20, mp.Hz, mp.Vector3(x_mon, y1, h_wg / 2.0), 1e-4))
    ref_mode = max(abs(sim_ref.get_eigenmode_coefficients(f_ref, [1], eig_parity=mp.NO_PARITY).alpha[0, 0, 0])**2, 1e-15)

    # Device
    geom_dev = [
        mp.Block(center=mp.Vector3(0, 0, -sz / 4.0), size=mp.Vector3(mp.inf, mp.inf, sz / 2.0), material=clad_mat),
        mp.Block(center=mp.Vector3(0, y1, h_wg / 2.0), size=mp.Vector3(mp.inf, w_wg, h_wg), material=core_mat),
    ]
    poly_wg2_in = make_ccw_linear_taper(x_c_start - l_in, x_c_start, y2_c, 0.02, w_wg)
    geom_dev.append(create_prism_exact(poly_wg2_in, h_wg, 0.0, core_mat))
    geom_dev.append(mp.Block(center=mp.Vector3((x_c_start + x_c_end) / 2.0, y2_c, h_wg / 2.0), size=mp.Vector3(l_c, w_wg, h_wg), material=core_mat))
    poly_wg2_bend = make_ccw_hermite_c2_sbend(x_c_end, x_bend_end, y2_c, y_out2, w_wg, w_wg)
    geom_dev.append(create_prism_exact(poly_wg2_bend, h_wg, 0.0, core_mat))
    geom_dev.append(mp.Block(center=mp.Vector3((x_bend_end + sx / 2.0) / 2.0, y_out2, h_wg / 2.0), size=mp.Vector3(sx / 2.0 - x_bend_end, w_wg, h_wg), material=core_mat))

    mat_pcm1 = pcm_cr if state.lower() == "crystalline" else pcm_am
    mat_pcm2 = pcm_am

    poly_pcm1_in = make_ccw_linear_taper(x_c_start, x_c_start + l_tip, y1, 0.03, w_wg)
    geom_dev.append(create_prism_exact(poly_pcm1_in, t_pcm, h_wg, mat_pcm1))
    geom_dev.append(mp.Block(center=mp.Vector3((x_c_start + l_tip + x_bend_end) / 2.0, y1, h_wg + t_pcm / 2.0), size=mp.Vector3(x_bend_end - (x_c_start + l_tip), w_wg, t_pcm), material=mat_pcm1))
    poly_pcm1_out = make_ccw_linear_taper(x_bend_end, x_taper_end, y1, w_wg, 0.03)
    geom_dev.append(create_prism_exact(poly_pcm1_out, t_pcm, h_wg, mat_pcm1))

    poly_pcm2_in = make_ccw_linear_taper(x_c_start, x_c_start + l_tip, y2_c, 0.03, w_wg)
    geom_dev.append(create_prism_exact(poly_pcm2_in, t_pcm, h_wg, mat_pcm2))
    if l_c > l_tip:
        geom_dev.append(mp.Block(center=mp.Vector3((x_c_start + l_tip + x_c_end) / 2.0, y2_c, h_wg + t_pcm / 2.0), size=mp.Vector3(l_c - l_tip, w_wg, t_pcm), material=mat_pcm2))
    poly_pcm2_bend = make_ccw_hermite_c2_sbend(x_c_end, x_bend_end, y2_c, y_out2, w_wg, w_wg)
    geom_dev.append(create_prism_exact(poly_pcm2_bend, t_pcm, h_wg, mat_pcm2))
    poly_pcm2_out = make_ccw_linear_taper(x_bend_end, x_taper_end, y_out2, w_wg, 0.03)
    geom_dev.append(create_prism_exact(poly_pcm2_out, t_pcm, h_wg, mat_pcm2))

    sim_dev = mp.Simulation(cell_size=cell, boundary_layers=pml_layers, geometry=geom_dev,
                            sources=sources, resolution=resolution, default_material=clad_mat)
    f_p1 = sim_dev.add_flux(fcen, 0, 1, mp.FluxRegion(center=mp.Vector3(x_mon, y1, h_wg / 2.0), size=mp.Vector3(0, mon_w, mon_h)))
    f_p2 = sim_dev.add_flux(fcen, 0, 1, mp.FluxRegion(center=mp.Vector3(x_mon, y_out2, h_wg / 2.0), size=mp.Vector3(0, mon_w, mon_h)))

    y_dest = y_out2 if state.lower() == "amorphous" else y1
    sim_dev.run(until_after_sources=mp.stop_when_fields_decayed(20, mp.Hz, mp.Vector3(x_mon, y_dest, h_wg / 2.0), 1e-4))

    c_p1 = abs(sim_dev.get_eigenmode_coefficients(f_p1, [1], eig_parity=mp.NO_PARITY).alpha[0, 0, 0])**2
    c_p2 = abs(sim_dev.get_eigenmode_coefficients(f_p2, [1], eig_parity=mp.NO_PARITY).alpha[0, 0, 0])**2

    t1 = c_p1 / ref_mode
    t2 = c_p2 / ref_mode
    return t1, t2

if __name__ == "__main__":
    print("=" * 75)
    print("3-POINT L_c SCAN (Gap = 103.7165 nm, Res = 15)")
    print("=" * 75)
    print(f"{'L_c (um)':<10} {'State 0 P2 (%)':<18} {'State 0 P1 (%)':<18} {'State 1 P1 (%)':<18} {'State 1 P2 (%)':<18}")
    print("-" * 75)
    for lc in [4.30, 4.45, 4.60]:
        t1_am, t2_am = run_point(lc, "amorphous")
        t1_cr, t2_cr = run_point(lc, "crystalline")
        print(f"{lc:<10.2f} {t2_am*100:<18.6f} {t1_am*100:<18.6f} {t1_cr*100:<18.6f} {t2_cr*100:<18.6f}")
