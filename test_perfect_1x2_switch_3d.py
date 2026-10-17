"""
PROJECT JANUS: FULL-WAVE 3D VECTOR MAXWELL MEEP FDTD SIMULATION
ULTRA-LOW-LOSS 1x2 NON-VOLATILE OPTICAL SWITCH (TARGET: IL <= 0.08 dB, XT <= -20 dB)
===================================================================================
Architecture: 1x2 Directional Coupler Switch (1 Input -> 2 Outputs)
  - Input: Port 0 on WG1 (y = +y_sep/2)
  - Output 1: Port 1 (Through / Bar on WG1 at y = +y_out/2)
  - Output 2: Port 2 (Cross on WG2 at y = -y_out/2)
Silicon Core: h = 220 nm, w1 = w2 = 240 nm, n = 3.480
Cladding: SiO2, n = 1.449
PCM: Sb2S3 (Amorphous n=2.70, Crystalline n=3.30, k_cr=8e-5)
  - WG1: Active PCM patch (switched: amorphous -> crystalline)
  - WG2: Matching amorphous PCM patch (provides exact phase matching in amorphous state)
Apodization: Parabolic adiabatic tapers (L_tip = 0.60 um, w_tip = 30 nm)
Fanout: Cosine S-bends (L_bend = 0.95 um, output pitch = 480 nm)
"""
import sys
import os
import math
import time
import json
import numpy as np
import meep as mp

def run_1x2_switch_3d(l_c=4.40, gap=0.040, state="amorphous", resolution=20):
    lambda_0 = 1.064
    fcen = 1.0 / lambda_0
    df = 0.1 * fcen

    w_wg = 0.240       # Silicon width
    h_wg = 0.220       # Silicon height
    t_pcm = 0.150      # Sb2S3 PCM thickness
    l_in = 0.350       # Input lead inside cell
    l_tip = 0.600      # Adiabatic taper length
    w_tip = 0.030      # Taper tip width
    l_bend = 0.950     # Cosine S-bend fanout length
    y_pitch_out = 0.480# Output port pitch
    l_out = 0.350      # Output lead inside cell
    iso_buf = 0.080    # Buffer to cell boundary

    y_sep = w_wg + gap
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

    # PCM Materials:
    # In Amorphous state: WG1 is amorphous (n=2.70), WG2 is amorphous (n=2.70)
    # In Crystalline state: WG1 is crystalline (n=3.30, k=8e-5), WG2 is amorphous (n=2.70)
    pcm_am = mp.Medium(index=2.70, D_conductivity=2*math.pi*fcen*1.0e-5)
    pcm_cr = mp.Medium(index=3.30, D_conductivity=2*math.pi*fcen*8.0e-5)

    mat_wg1_pcm = pcm_am if state == "amorphous" else pcm_cr
    mat_wg2_pcm = pcm_am  # WG2 is always amorphous matching layer

    # Reference straight waveguide monitor
    x_mon = x_bend_end + 0.15
    mon_w = w_wg * 2.5
    mon_h = h_wg * 2.5

    sources = [
        mp.EigenModeSource(
            src=mp.GaussianSource(fcen, fwidth=df),
            center=mp.Vector3(x_left + dpml + 0.10, y_wg1, h_wg / 2.0),
            size=mp.Vector3(0, w_wg * 3.0, h_wg * 3.0),
            eig_band=1,
            direction=mp.X
        )
    ]

    print(f"[{state.upper()}] 1. Simulating Reference Straight Waveguide for Flux Normalization...")
    geom_ref = [
        mp.Block(center=mp.Vector3(0, 0, -sz/4.0), size=mp.Vector3(mp.inf, mp.inf, sz/2.0), material=clad_mat),
        mp.Block(center=mp.Vector3(0, y_wg1, h_wg/2.0), size=mp.Vector3(mp.inf, w_wg, h_wg), material=core_mat)
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
    print(f"[{state.upper()}] Reference Poynting Flux = {p_ref:.6e}")

    # C2-Continuous Cosine S-bend generator
    def make_cosine_sbend(x_start, x_end, y_start, y_end, w, n_pts=30):
        dx = x_end - x_start
        dy = y_end - y_start
        u = np.linspace(0, 1, n_pts)
        x = x_start + u * dx
        y = y_start + dy * (0.5 * (1.0 - np.cos(np.pi * u)))
        top = [mp.Vector3(xi, yi + w / 2.0, 0) for xi, yi in zip(x, y)]
        bot = [mp.Vector3(xi, yi - w / 2.0, 0) for xi, yi in reversed(list(zip(x, y)))]
        return top + bot

    # Parabolic adiabatic PCM taper generator
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

    x_c_center = (x_in_end + x_coupler_end) / 2.0
    pcm_len = l_c + 2 * 0.40  # Tapers extend slightly to ensure adiabatic transition

    patch_pts1 = make_parabolic_pcm(x_c_center, y_wg1, pcm_len, w_wg, l_tip=l_tip, w_tip=w_tip)
    cx_p1 = float(np.mean([pt.x for pt in patch_pts1]))
    cy_p1 = float(np.mean([pt.y for pt in patch_pts1]))

    patch_pts2 = make_parabolic_pcm(x_c_center, y_wg2, pcm_len, w_wg, l_tip=l_tip, w_tip=w_tip)
    cx_p2 = float(np.mean([pt.x for pt in patch_pts2]))
    cy_p2 = float(np.mean([pt.y for pt in patch_pts2]))

    geom_switch = [
        # Substrate
        mp.Block(center=mp.Vector3(0, 0, -sz/4.0), size=mp.Vector3(mp.inf, mp.inf, sz/2.0), material=clad_mat),
        # 1. Single Input Waveguide (Port 0 on WG1)
        mp.Block(center=mp.Vector3((x_left + x_in_end)/2.0, y_wg1, h_wg/2.0), size=mp.Vector3(x_in_end - x_left, w_wg, h_wg), material=core_mat),
        # 2. Coupler Section (Identical silicon cores)
        mp.Block(center=mp.Vector3(x_c_center, y_wg1, h_wg/2.0), size=mp.Vector3(l_c, w_wg, h_wg), material=core_mat),
        mp.Block(center=mp.Vector3(x_c_center, y_wg2, h_wg/2.0), size=mp.Vector3(l_c, w_wg, h_wg), material=core_mat),
        # 3. Active PCM on WG1 and Matching Layer on WG2
        mp.Prism(vertices=patch_pts1, height=t_pcm, axis=mp.Vector3(0, 0, 1),
                 center=mp.Vector3(cx_p1, cy_p1, h_wg + t_pcm/2.0), material=mat_wg1_pcm),
        mp.Prism(vertices=patch_pts2, height=t_pcm, axis=mp.Vector3(0, 0, 1),
                 center=mp.Vector3(cx_p2, cy_p2, h_wg + t_pcm/2.0), material=mat_wg2_pcm),
        # 4. Ultra-Smooth Cosine S-bends to Output Ports
        mp.Prism(vertices=make_cosine_sbend(x_coupler_end, x_bend_end, y_wg1, y_out1, w_wg), height=h_wg, axis=mp.Vector3(0,0,1), material=core_mat),
        mp.Prism(vertices=make_cosine_sbend(x_coupler_end, x_bend_end, y_wg2, y_out2, w_wg), height=h_wg, axis=mp.Vector3(0,0,1), material=core_mat),
        # 5. Output Leads (Port 1 and Port 2)
        mp.Block(center=mp.Vector3((x_bend_end + x_right)/2.0, y_out1, h_wg/2.0), size=mp.Vector3(x_right - x_bend_end, w_wg, h_wg), material=core_mat),
        mp.Block(center=mp.Vector3((x_bend_end + x_right)/2.0, y_out2, h_wg/2.0), size=mp.Vector3(x_right - x_bend_end, w_wg, h_wg), material=core_mat),
    ]

    print(f"[{state.upper()}] 2. Simulating 3D Switch Device...")
    sim_switch = mp.Simulation(
        cell_size=cell, boundary_layers=pml_layers, geometry=geom_switch,
        sources=sources, resolution=resolution, default_material=clad_mat
    )
    flux_p1 = sim_switch.add_flux(fcen, 0, 1, mp.FluxRegion(center=mp.Vector3(x_mon, y_out1, h_wg/2.0), size=mp.Vector3(0, mon_w, mon_h)))
    flux_p2 = sim_switch.add_flux(fcen, 0, 1, mp.FluxRegion(center=mp.Vector3(x_mon, y_out2, h_wg/2.0), size=mp.Vector3(0, mon_w, mon_h)))

    y_dest = y_out2 if state == "amorphous" else y_out1
    sim_switch.run(until_after_sources=mp.stop_when_fields_decayed(20, mp.Hz, mp.Vector3(x_mon, y_dest, h_wg/2.0), 1e-3))

    p1 = max(mp.get_fluxes(flux_p1)[0] / p_ref, 1e-12)
    p2 = max(mp.get_fluxes(flux_p2)[0] / p_ref, 1e-12)
    tot = p1 + p2

    # In 1x2 switch:
    # Amorphous State: Signal is Port 2 (Cross), Leakage is Port 1 (Bar)
    # Crystalline State: Signal is Port 1 (Bar), Leakage is Port 2 (Cross)
    sig = p2 if state == "amorphous" else p1
    leak = p1 if state == "amorphous" else p2
    il_db = -10.0 * math.log10(max(sig, 1e-6))
    xt_db = 10.0 * math.log10(max(leak / max(sig, 1e-6), 1e-6))

    return {
        "state": state, "l_c": l_c, "gap_nm": gap * 1000,
        "cell_len": cell_len, "cell_w": cell_w, "area_um2": cell_area,
        "p1_bar_pct": p1 * 100.0, "p2_cross_pct": p2 * 100.0, "tot_trans_pct": tot * 100.0,
        "il_db": il_db, "xt_db": xt_db
    }

if __name__ == "__main__":
    st = sys.argv[1] if len(sys.argv) > 1 else "amorphous"
    lc = float(sys.argv[2]) if len(sys.argv) > 2 else 4.40
    gp = float(sys.argv[3]) if len(sys.argv) > 3 else 0.040

    res = run_1x2_switch_3d(l_c=lc, gap=gp, state=st)
    print("\n" + "=" * 70)
    print(f"PROJECT JANUS: HIGH-EFFICIENCY 1x2 SWITCH 3D FDTD RESULT ({st.upper()}):")
    print(f"  Footprint Area  : {res['area_um2']:.3f} um^2 ({res['cell_len']:.2f} um x {res['cell_w']:.2f} um)")
    print(f"  Port 1 (Bar)    : {res['p1_bar_pct']:.2f}%")
    print(f"  Port 2 (Cross)  : {res['p2_cross_pct']:.2f}%")
    print(f"  Total Trans     : {res['tot_trans_pct']:.2f}%")
    print(f"  Insertion Loss  : {res['il_db']:.3f} dB (Target: <= 0.08 dB)")
    print(f"  Crosstalk       : {res['xt_db']:.2f} dB (Target: <= -20.0 dB)")
    print("=" * 70)
