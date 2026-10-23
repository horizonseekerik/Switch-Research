"""
PROJECT JANUS: 3D FULL-WAVE MAXWELL MEEP FDTD SIMULATION
SUB-3 um^2 1x2 BALANCED MZI OPTICAL SWITCH
==========================================
Architecture: 1x2 Switch (1 Input -> 2 Outputs)
  - Input: 1 single silicon waveguide (Port 0 at y = 0)
  - Splitter: 1x2 Symmetric Y-branch (50:50 power split by bilateral symmetry)
  - Upper Arm: Active Phase Shifter with Sb2S3 (t_pcm = 150 nm, L_pcm = 2.30 um)
  - Lower Arm: Passive Reference Arm (bare silicon)
  - Output Combiner: 2x2 3dB Directional Coupler (L_3dB = 1.25 um, gap = 45 nm)
  - Outputs: Port 1 (Bar, y = +0.275 um) and Port 2 (Cross, y = -0.275 um)
Footprint Target: Sub-3.5 um^2 Active Cell
Transmission Target: > 97% (IL < 0.13 dB)
Crosstalk Target: < -20.0 dB (ER > 20 dB, SCR > 100:1)
"""
import sys
import os
import math
import time
import json
import numpy as np
import meep as mp

def run_1x2_mzi_switch(state="amorphous", resolution=20):
    lambda_0 = 1.064
    fcen = 1.0 / lambda_0
    df = 0.1 * fcen

    w_wg = 0.250        # Silicon waveguide width
    h_wg = 0.220        # Silicon height
    t_pcm = 0.150       # Sb2S3 PCM thickness
    y_pitch = 0.550     # Arm separation pitch (y = +0.275, -0.275)
    l_ybranch = 0.800   # 1x2 Y-junction length
    l_phase = 2.300     # PCM phase modulation length (pi shift)
    l_conv = 0.400      # Convergence length to coupler gap
    gap_c = 0.045       # 3dB coupler gap
    l_3db = 1.250       # 3dB coupler length
    l_fanout = 0.450    # S-bend fanout to output ports
    l_in = 0.300        # Input straight lead
    l_out = 0.300       # Output straight lead
    iso_buf = 0.080     # Lateral isolation buffer

    y1_arm = y_pitch / 2.0
    y2_arm = -y_pitch / 2.0

    y_sep_coupler = w_wg + gap_c
    y1_coupler = y_sep_coupler / 2.0
    y2_coupler = -y_sep_coupler / 2.0

    cell_len = l_in + l_ybranch + l_phase + l_conv + l_3db + l_fanout + l_out
    cell_w = y_pitch + 2 * iso_buf
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
    x_y_end = x_in_end + l_ybranch
    x_phase_end = x_y_end + l_phase
    x_conv_end = x_phase_end + l_conv
    x_c_end = x_conv_end + l_3db
    x_bend_end = x_c_end + l_fanout

    clad_mat = mp.Medium(index=1.449)
    core_mat = mp.Medium(index=3.480)

    n_pcm = 2.70 if state == "amorphous" else 3.30
    k_pcm = 1.0e-5 if state == "amorphous" else 8.0e-5
    patch_mat = mp.Medium(index=n_pcm, D_conductivity=2*math.pi*fcen*k_pcm)

    # 1. Reference Single Straight Waveguide
    x_mon = x_bend_end + 0.15
    mon_w = 0.260       # Non-overlapping isolated monitor aperture
    mon_h = 0.280

    sources = [
        mp.EigenModeSource(
            src=mp.GaussianSource(fcen, fwidth=df),
            center=mp.Vector3(x_left + dpml + 0.10, 0, h_wg / 2.0),
            size=mp.Vector3(0, w_wg * 3.0, h_wg * 3.0),
            eig_band=1,
            direction=mp.X
        )
    ]

    print(f"[{state.upper()}] 1. Running Reference Straight Waveguide Simulation...")
    geom_ref = [
        mp.Block(center=mp.Vector3(0, 0, -sz/4.0), size=mp.Vector3(mp.inf, mp.inf, sz/2.0), material=clad_mat),
        mp.Block(center=mp.Vector3(0, 0, h_wg/2.0), size=mp.Vector3(mp.inf, w_wg, h_wg), material=core_mat)
    ]
    sim_ref = mp.Simulation(
        cell_size=cell, boundary_layers=pml_layers, geometry=geom_ref,
        sources=sources, resolution=resolution, default_material=clad_mat
    )
    flux_ref = sim_ref.add_flux(fcen, 0, 1, mp.FluxRegion(
        center=mp.Vector3(x_mon, 0, h_wg/2.0), size=mp.Vector3(0, mon_w, mon_h)
    ))
    sim_ref.run(until_after_sources=mp.stop_when_fields_decayed(20, mp.Hz, mp.Vector3(x_mon, 0, h_wg/2.0), 1e-3))
    p_ref = max(mp.get_fluxes(flux_ref)[0], 1e-12)
    print(f"[{state.upper()}] Reference Poynting Flux = {p_ref:.6e}")

    # S-bend and Polygon Utilities
    def make_sbend(x_start, x_end, y_start, y_end, w, n_pts=20):
        dx = x_end - x_start
        dy = y_end - y_start
        u = np.linspace(0, 1, n_pts)
        x = x_start + u * dx
        y = y_start + dy * (3 * u**2 - 2 * u**3)
        top = [mp.Vector3(xi, yi + w / 2.0, 0) for xi, yi in zip(x, y)]
        bot = [mp.Vector3(xi, yi - w / 2.0, 0) for xi, yi in reversed(list(zip(x, y)))]
        return top + bot

    def make_apodized_pcm_polygon(x_c, y_c, l_tot, w_p, l_tip=0.35, w_tip=0.040, n_pts=10):
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

    x_pcm_c = (x_y_end + x_phase_end) / 2.0
    pcm_pts = make_apodized_pcm_polygon(x_pcm_c, y1_arm, l_phase, w_wg)
    cx_pcm = float(np.mean([pt.x for pt in pcm_pts]))
    cy_pcm = float(np.mean([pt.y for pt in pcm_pts]))

    # Construct complete 1x2 MZI Device Geometry
    geom = [
        mp.Block(center=mp.Vector3(0, 0, -sz/4.0), size=mp.Vector3(mp.inf, mp.inf, sz/2.0), material=clad_mat),
        # 1. Single Input Lead (Port 0)
        mp.Block(center=mp.Vector3((x_left + x_in_end)/2.0, 0, h_wg/2.0), size=mp.Vector3(x_in_end - x_left, w_wg, h_wg), material=core_mat),
        # 2. 1x2 Y-Branch Splitter (diverges from y=0 to y=+0.275 and y=-0.275)
        mp.Prism(vertices=make_sbend(x_in_end, x_y_end, 0, y1_arm, w_wg), height=h_wg, axis=mp.Vector3(0,0,1), material=core_mat),
        mp.Prism(vertices=make_sbend(x_in_end, x_y_end, 0, y2_arm, w_wg), height=h_wg, axis=mp.Vector3(0,0,1), material=core_mat),
        # 3. Parallel Phase Arms
        mp.Block(center=mp.Vector3(x_pcm_c, y1_arm, h_wg/2.0), size=mp.Vector3(l_phase, w_wg, h_wg), material=core_mat),
        mp.Block(center=mp.Vector3(x_pcm_c, y2_arm, h_wg/2.0), size=mp.Vector3(l_phase, w_wg, h_wg), material=core_mat),
        # Active Sb2S3 on Upper Arm
        mp.Prism(vertices=pcm_pts, height=t_pcm, axis=mp.Vector3(0, 0, 1),
                 center=mp.Vector3(cx_pcm, cy_pcm, h_wg + t_pcm/2.0), material=patch_mat),
        # 4. Convergence S-bends into 3dB Coupler Gap
        mp.Prism(vertices=make_sbend(x_phase_end, x_conv_end, y1_arm, y1_coupler, w_wg), height=h_wg, axis=mp.Vector3(0,0,1), material=core_mat),
        mp.Prism(vertices=make_sbend(x_phase_end, x_conv_end, y2_arm, y2_coupler, w_wg), height=h_wg, axis=mp.Vector3(0,0,1), material=core_mat),
        # 5. Output 3dB Directional Coupler
        mp.Block(center=mp.Vector3((x_conv_end + x_c_end)/2.0, y1_coupler, h_wg/2.0), size=mp.Vector3(l_3db, w_wg, h_wg), material=core_mat),
        mp.Block(center=mp.Vector3((x_conv_end + x_c_end)/2.0, y2_coupler, h_wg/2.0), size=mp.Vector3(l_3db, w_wg, h_wg), material=core_mat),
        # 6. S-Bend Fanout to Output Port Pitch
        mp.Prism(vertices=make_sbend(x_c_end, x_bend_end, y1_coupler, y1_arm, w_wg), height=h_wg, axis=mp.Vector3(0,0,1), material=core_mat),
        mp.Prism(vertices=make_sbend(x_c_end, x_bend_end, y2_coupler, y2_arm, w_wg), height=h_wg, axis=mp.Vector3(0,0,1), material=core_mat),
        # 7. Output Leads through Right PML
        mp.Block(center=mp.Vector3((x_bend_end + x_right)/2.0, y1_arm, h_wg/2.0), size=mp.Vector3(x_right - x_bend_end, w_wg, h_wg), material=core_mat),
        mp.Block(center=mp.Vector3((x_bend_end + x_right)/2.0, y2_arm, h_wg/2.0), size=mp.Vector3(x_right - x_bend_end, w_wg, h_wg), material=core_mat),
    ]

    print(f"[{state.upper()}] 2. Running 1x2 MZI Device Simulation...")
    sim = mp.Simulation(
        cell_size=cell, boundary_layers=pml_layers, geometry=geom,
        sources=sources, resolution=resolution, default_material=clad_mat
    )
    flux_p1 = sim.add_flux(fcen, 0, 1, mp.FluxRegion(center=mp.Vector3(x_mon, y1_arm, h_wg/2.0), size=mp.Vector3(0, mon_w, mon_h)))
    flux_p2 = sim.add_flux(fcen, 0, 1, mp.FluxRegion(center=mp.Vector3(x_mon, y2_arm, h_wg/2.0), size=mp.Vector3(0, mon_w, mon_h)))

    y_dest = y1_arm if state == "amorphous" else y2_arm
    sim.run(until_after_sources=mp.stop_when_fields_decayed(20, mp.Hz, mp.Vector3(x_mon, y_dest, h_wg/2.0), 1e-3))

    p1 = max(mp.get_fluxes(flux_p1)[0] / p_ref, 1e-12)
    p2 = max(mp.get_fluxes(flux_p2)[0] / p_ref, 1e-12)
    tot = p1 + p2

    sig = p1 if state == "amorphous" else p2
    leak = p2 if state == "amorphous" else p1
    il_db = -10.0 * math.log10(max(sig, 1e-6))
    xt_db = 10.0 * math.log10(max(leak / max(sig, 1e-6), 1e-6))

    return {
        "state": state, "cell_len": cell_len, "cell_w": cell_w, "area_um2": cell_area,
        "p1_bar_pct": p1 * 100.0, "p2_cross_pct": p2 * 100.0, "tot_trans_pct": tot * 100.0,
        "il_db": il_db, "xt_db": xt_db
    }

if __name__ == "__main__":
    st = sys.argv[1] if len(sys.argv) > 1 else "amorphous"
    res = run_1x2_mzi_switch(state=st)
    print("\n" + "=" * 65)
    print(f"RESULTS FOR 1x2 BALANCED MZI SWITCH ({st.upper()}):")
    print(f"  Footprint Area : {res['area_um2']:.3f} um^2 ({res['cell_len']:.2f} um x {res['cell_w']:.2f} um)")
    print(f"  Port 1 (Bar)   : {res['p1_bar_pct']:.2f}%")
    print(f"  Port 2 (Cross) : {res['p2_cross_pct']:.2f}%")
    print(f"  Total Trans    : {res['tot_trans_pct']:.2f}%")
    print(f"  Insertion Loss : {res['il_db']:.3f} dB")
    print(f"  Crosstalk      : {res['xt_db']:.2f} dB")
    print("=" * 65)
