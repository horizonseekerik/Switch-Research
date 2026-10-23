"""
Directly query epsilon in MEEP simulation
"""
import math
import numpy as np
import meep as mp

def check_eps(state="crystalline"):
    lambda_0 = 1.064
    fcen = 1.0 / lambda_0
    w_wg = 0.240
    h_wg = 0.220
    t_pcm = 0.150
    gap = 0.040
    l_c = 4.40
    y_pitch_out = 0.480
    iso_buf = 0.080
    l_in, l_out, l_bend = 0.350, 0.350, 0.950

    y_sep = w_wg + gap
    y_wg1 = y_sep / 2.0
    y_wg2 = -y_sep / 2.0
    y_out1 = y_pitch_out / 2.0
    y_out2 = -y_pitch_out / 2.0

    cell_len = l_in + l_c + l_bend + l_out
    cell_w = y_pitch_out + 2 * iso_buf
    dpml, pad_x, pad_y, pad_z = 0.50, 0.35, 0.35, 0.35
    sx = cell_len + 2 * (dpml + pad_x)
    sy = cell_w + 2 * (dpml + pad_y)
    sz = h_wg + t_pcm + 2 * (dpml + pad_z)

    cell = mp.Vector3(sx, sy, sz)
    clad_mat = mp.Medium(index=1.449)
    core_mat = mp.Medium(index=3.480)

    pcm_am = mp.Medium(index=2.70)
    pcm_cr = mp.Medium(index=3.30)
    mat_wg1_pcm = pcm_am if state == "amorphous" else pcm_cr
    mat_wg2_pcm = pcm_am

    x_min_cell = -cell_len / 2.0
    x_in_end = x_min_cell + l_in
    x_coupler_end = x_in_end + l_c
    x_bend_end = x_coupler_end + l_bend
    x_c_center = (x_in_end + x_coupler_end) / 2.0

    def make_parabolic_pcm(x_c, y_c, l_tot, w_p, l_tip=0.60, w_tip=0.030, n_pts=20):
        x_lt = x_c - l_tot / 2.0
        x_lb = x_lt + l_tip
        x_rb = x_c + l_tot / 2.0 - l_tip
        u_l = np.linspace(0, 1, n_pts)[:-1]
        xl = x_lt + u_l * l_tip
        wl = w_tip + (w_p - w_tip) * (2 * u_l - u_l**2)
        u_r = np.linspace(0, 1, n_pts)[1:]
        xr = x_rb + u_r * l_tip
        wr = w_tip + (w_p - w_tip) * (1.0 - u_r**2)
        top = [mp.Vector3(xi, y_c + wi/2.0, 0) for xi, wi in zip(xl, wl)]
        top.append(mp.Vector3(x_lb, y_c + w_p/2.0, 0))
        top.append(mp.Vector3(x_rb, y_c + w_p/2.0, 0))
        top.extend([mp.Vector3(xi, y_c + wi/2.0, 0) for xi, wi in zip(xr, wr)])
        bot = [mp.Vector3(xi, y_c - wi/2.0, 0) for xi, wi in reversed(list(zip(xr, wr)))]
        bot.append(mp.Vector3(x_rb, y_c - w_p/2.0, 0))
        bot.append(mp.Vector3(x_lb, y_c - w_p/2.0, 0))
        bot.extend([mp.Vector3(xi, y_c - wi/2.0, 0) for xi, wi in reversed(list(zip(xl, wl)))])
        return top + bot

    pcm_len = l_c + 2 * 0.40
    patch_pts1 = make_parabolic_pcm(x_c_center, y_wg1, pcm_len, w_wg)
    cx_p1 = float(np.mean([pt.x for pt in patch_pts1]))
    cy_p1 = float(np.mean([pt.y for pt in patch_pts1]))

    patch_pts2 = make_parabolic_pcm(x_c_center, y_wg2, pcm_len, w_wg)
    cx_p2 = float(np.mean([pt.x for pt in patch_pts2]))
    cy_p2 = float(np.mean([pt.y for pt in patch_pts2]))

    geom = [
        mp.Block(center=mp.Vector3(0, 0, -sz/4.0), size=mp.Vector3(mp.inf, mp.inf, sz/2.0), material=clad_mat),
        mp.Block(center=mp.Vector3(x_c_center, y_wg1, h_wg/2.0), size=mp.Vector3(l_c, w_wg, h_wg), material=core_mat),
        mp.Block(center=mp.Vector3(x_c_center, y_wg2, h_wg/2.0), size=mp.Vector3(l_c, w_wg, h_wg), material=core_mat),
        mp.Prism(vertices=patch_pts1, height=t_pcm, axis=mp.Vector3(0, 0, 1),
                 center=mp.Vector3(cx_p1, cy_p1, h_wg + t_pcm/2.0), material=mat_wg1_pcm),
        mp.Prism(vertices=patch_pts2, height=t_pcm, axis=mp.Vector3(0, 0, 1),
                 center=mp.Vector3(cx_p2, cy_p2, h_wg + t_pcm/2.0), material=mat_wg2_pcm),
    ]

    sim = mp.Simulation(cell_size=cell, geometry=geom, resolution=20, default_material=clad_mat)
    sim.init_sim()

    xc = x_c_center
    pts = {
        "WG1 Core (Si)": mp.Vector3(xc, y_wg1, h_wg/2.0),
        "WG1 PCM Layer": mp.Vector3(xc, y_wg1, h_wg + t_pcm/2.0),
        "WG2 Core (Si)": mp.Vector3(xc, y_wg2, h_wg/2.0),
        "WG2 PCM Layer": mp.Vector3(xc, y_wg2, h_wg + t_pcm/2.0),
        "Coupling Gap ": mp.Vector3(xc, 0, h_wg/2.0),
    }
    print(f"\nEpsilon values for {state.upper()} state at coupler center x = {xc:.2f}:")
    for name, pt in pts.items():
        eps = sim.get_epsilon_point(pt)
        n = math.sqrt(eps)
        print(f"  {name}: eps = {eps:.4f} -> n = {n:.4f}")

if __name__ == "__main__":
    check_eps("amorphous")
    check_eps("crystalline")
