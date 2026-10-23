"""
PROJECT JANUS: HIGH-EFFICIENCY 1x2 OPTICAL SWITCH (TARGET: T > 97%, XT < -20 dB)
===============================================================================
Architecture: 1x2 Directional Coupler Switch (1 Input -> 2 Outputs)
  - Port 0: Single Input Waveguide (WG1 at y = +y_sep/2)
  - WG1 (Active): With apodized Sb2S3 PCM (L_tip = 0.50 um, w_tip = 30 nm)
  - WG2 (Passive): Bare silicon, receives 100% in Cross state
  - Output Ports: Port 1 (Bar) and Port 2 (Cross) with gentle L_bend = 0.85 um fanout
"""
import sys
import os
import math
import time
import numpy as np
import meep as mp

def simulate_1x2_switch(w1=0.240, w2=0.240, gap=0.035, l_c=4.50, t_pcm=0.150, 
                        state="amorphous", resolution=20):
    lambda_0 = 1.064
    fcen = 1.0 / lambda_0
    df = 0.1 * fcen

    h_wg = 0.220
    l_in = 0.350
    l_ext = 0.400
    l_bend = 0.850
    y_pitch_out = 0.500
    l_out = 0.350
    iso_buf = 0.080

    y_sep = (w1 + w2) / 2.0 + gap
    y_wg1 = y_sep / 2.0
    y_wg2 = -y_sep / 2.0
    y_out1 = y_pitch_out / 2.0
    y_out2 = -y_pitch_out / 2.0

    cell_len = l_in + l_c + l_bend + l_out
    cell_w = y_pitch_out + 2 * iso_buf
    cell_area = cell_len * cell_w

    dpml = 0.50
    pad_x = 0.35
    pad_y = 0.35
    pad_z = 0.35
    sx = cell_len + 2 * (dpml + pad_x)
    sy = cell_w + 2 * (dpml + pad_y)
    sz = h_wg + t_pcm + 2 * (dpml + pad_z)

    cell = mp.Vector3(sx, sy, sz)
    pml_layers = [mp.PML(dpml)]

    x_left = -sx / 2.0
    x_right = sx / 2.0
    x_min_cell = -cell_len / 2.0
    x_in_end = x_min_cell + l_in
    x_coupler_end = x_in_end + l_c
    x_bend_end = x_coupler_end + l_bend

    clad_mat = mp.Medium(index=1.449)
    core_mat = mp.Medium(index=3.480)

    n_pcm = 2.70 if state == "amorphous" else 3.30
    k_pcm = 1.0e-5 if state == "amorphous" else 8.0e-5
    patch_mat = mp.Medium(index=n_pcm, D_conductivity=2*math.pi*fcen*k_pcm)

    # Reference straight waveguide
    x_mon = x_bend_end + 0.15
    mon_w = max(w1, w2) * 2.5
    mon_h = h_wg * 2.5

    sources = [
        mp.EigenModeSource(
            src=mp.GaussianSource(fcen, fwidth=df),
            center=mp.Vector3(x_left + dpml + 0.10, y_wg1, h_wg / 2.0),
            size=mp.Vector3(0, w1 * 3.0, h_wg * 3.0),
            eig_band=1,
            direction=mp.X
        )
    ]

    geom_ref = [
        mp.Block(center=mp.Vector3(0, 0, -sz/4.0), size=mp.Vector3(mp.inf, mp.inf, sz/2.0), material=clad_mat),
        mp.Block(center=mp.Vector3(0, y_wg1, h_wg/2.0), size=mp.Vector3(mp.inf, w1, h_wg), material=core_mat)
    ]
    sim_ref = mp.Simulation(
        cell_size=cell, boundary_layers=pml_layers, geometry=geom_ref,
        sources=sources, resolution=resolution, default_material=clad_mat
    )
    flux_ref = sim_ref.add_flux(fcen, 0, 1, mp.FluxRegion(
        center=mp.Vector3(x_mon, y_wg1, h_wg/2.0), size=mp.Vector3(0, mon_w, mon_h)
    ))
    sim_ref.run(until_after_sources=mp.stop_when_fields_decayed(20, mp.Hz, mp.Vector3(x_mon, y_wg1, h_wg/2.0), 1e-3))
    p_ref = max(mp.get_fluxes(flux_ref)[0], 1e-12)

    def make_sbend(x_start, x_end, y_start, y_end, w, n_pts=25):
        dx = x_end - x_start
        dy = y_end - y_start
        u = np.linspace(0, 1, n_pts)
        x = x_start + u * dx
        y = y_start + dy * (3 * u**2 - 2 * u**3)
        top = [mp.Vector3(xi, yi + w / 2.0, 0) for xi, yi in zip(x, y)]
        bot = [mp.Vector3(xi, yi - w / 2.0, 0) for xi, yi in reversed(list(zip(x, y)))]
        return top + bot

    def make_apodized_pcm(x_c, y_c, l_tot, w_p, l_tip=0.50, w_tip=0.030, n_pts=15):
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

    x_c_center = (x_in_end + x_coupler_end) / 2.0
    patch_pts = make_apodized_pcm(x_c_center, y_wg1, l_c + 2*l_ext, w1)
    cx_p = float(np.mean([pt.x for pt in patch_pts]))
    cy_p = float(np.mean([pt.y for pt in patch_pts]))

    geom = [
        mp.Block(center=mp.Vector3(0, 0, -sz/4.0), size=mp.Vector3(mp.inf, mp.inf, sz/2.0), material=clad_mat),
        # 1. Single Input Lead on WG1
        mp.Block(center=mp.Vector3((x_left + x_in_end)/2.0, y_wg1, h_wg/2.0), size=mp.Vector3(x_in_end - x_left, w1, h_wg), material=core_mat),
        # 2. Coupler Section
        mp.Block(center=mp.Vector3(x_c_center, y_wg1, h_wg/2.0), size=mp.Vector3(l_c, w1, h_wg), material=core_mat),
        mp.Block(center=mp.Vector3(x_c_center, y_wg2, h_wg/2.0), size=mp.Vector3(l_c, w2, h_wg), material=core_mat),
        # Active PCM on WG1
        mp.Prism(vertices=patch_pts, height=t_pcm, axis=mp.Vector3(0, 0, 1),
                 center=mp.Vector3(cx_p, cy_p, h_wg + t_pcm/2.0), material=patch_mat),
        # 3. Gentle S-bends to Output Ports
        mp.Prism(vertices=make_sbend(x_coupler_end, x_bend_end, y_wg1, y_out1, w1), height=h_wg, axis=mp.Vector3(0,0,1), material=core_mat),
        mp.Prism(vertices=make_sbend(x_coupler_end, x_bend_end, y_wg2, y_out2, w2), height=h_wg, axis=mp.Vector3(0,0,1), material=core_mat),
        # 4. Output Leads
        mp.Block(center=mp.Vector3((x_bend_end + x_right)/2.0, y_out1, h_wg/2.0), size=mp.Vector3(x_right - x_bend_end, w1, h_wg), material=core_mat),
        mp.Block(center=mp.Vector3((x_bend_end + x_right)/2.0, y_out2, h_wg/2.0), size=mp.Vector3(x_right - x_bend_end, w2, h_wg), material=core_mat),
    ]

    sim = mp.Simulation(
        cell_size=cell, boundary_layers=pml_layers, geometry=geom,
        sources=sources, resolution=resolution, default_material=clad_mat
    )
    flux_p1 = sim.add_flux(fcen, 0, 1, mp.FluxRegion(center=mp.Vector3(x_mon, y_out1, h_wg/2.0), size=mp.Vector3(0, mon_w, mon_h)))
    flux_p2 = sim.add_flux(fcen, 0, 1, mp.FluxRegion(center=mp.Vector3(x_mon, y_out2, h_wg/2.0), size=mp.Vector3(0, mon_w, mon_h)))

    y_dest = y_out2 if state == "amorphous" else y_out1
    sim.run(until_after_sources=mp.stop_when_fields_decayed(20, mp.Hz, mp.Vector3(x_mon, y_dest, h_wg/2.0), 1e-3))

    p1 = max(mp.get_fluxes(flux_p1)[0] / p_ref, 1e-12)
    p2 = max(mp.get_fluxes(flux_p2)[0] / p_ref, 1e-12)
    tot = p1 + p2

    sig = p2 if state == "amorphous" else p1
    leak = p1 if state == "amorphous" else p2
    il_db = -10.0 * math.log10(max(sig, 1e-6))
    xt_db = 10.0 * math.log10(max(leak / max(sig, 1e-6), 1e-6))

    return {
        "state": state, "cell_len": cell_len, "cell_w": cell_w, "area_um2": cell_area,
        "p1_bar_pct": p1 * 100.0, "p2_cross_pct": p2 * 100.0, "tot_trans_pct": tot * 100.0,
        "il_db": il_db, "xt_db": xt_db
    }

if __name__ == "__main__":
    st = sys.argv[1] if len(sys.argv) > 1 else "amorphous"
    res = simulate_1x2_switch(state=st)
    print("\n" + "=" * 65)
    print(f"HIGH-EFFICIENCY 1x2 SWITCH RESULTS ({st.upper()}):")
    print(f"  Footprint Area : {res['area_um2']:.3f} um^2 ({res['cell_len']:.2f} um x {res['cell_w']:.2f} um)")
    print(f"  Port 1 (Bar)   : {res['p1_bar_pct']:.2f}%")
    print(f"  Port 2 (Cross) : {res['p2_cross_pct']:.2f}%")
    print(f"  Total Trans    : {res['tot_trans_pct']:.2f}%")
    print(f"  Insertion Loss : {res['il_db']:.3f} dB")
    print(f"  Crosstalk      : {res['xt_db']:.2f} dB")
    print("=" * 65)
