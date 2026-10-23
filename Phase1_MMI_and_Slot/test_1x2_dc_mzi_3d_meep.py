"""
PROJECT JANUS: 3D FULL-WAVE MAXWELL MEEP FDTD SIMULATION
1x2 BALANCED DIRECTIONAL-COUPLER MZI SWITCH (DC-MZI)
====================================================
Architecture: 1x2 Switch (1 Input -> 2 Outputs)
  - Single Input on WG1 (Port 0)
  - Stage 1: Input 3dB Directional Coupler (50:50 power split, pi/2 phase difference)
  - Stage 2: Phase Shifter Arms with gap expanded to prevent coupling
             Arm 1: Sb2S3 PCM (t_pcm = 150 nm, L_pcm = 2.30 um, apodized tapers)
             Arm 2: Bare silicon reference arm
  - Stage 3: Output 3dB Directional Coupler (Interferometric Combiner)
  - Stage 4: Output Leads to Port 1 and Port 2
Formulas:
  P1 = sin^2(Delta_phi / 2)
  P2 = cos^2(Delta_phi / 2)
State 0 (Amorphous, Delta_phi = 0):  P2 = 100%, P1 = 0%  (Cross)
State 1 (Crystalline, Delta_phi = pi): P1 = 100%, P2 = 0%  (Bar)
Target: Transmission > 97%, Crosstalk < -20.0 dB
"""
import sys
import os
import math
import time
import json
import numpy as np
import meep as mp

def run_1x2_dc_mzi(state="amorphous", resolution=20):
    lambda_0 = 1.064
    fcen = 1.0 / lambda_0
    df = 0.1 * fcen

    w_wg = 0.250        # Waveguide width
    h_wg = 0.220        # Waveguide height
    t_pcm = 0.150       # Sb2S3 thickness
    gap_c = 0.045       # 3dB coupler gap
    l_3db = 1.250       # 3dB coupler length
    gap_phase = 0.300   # Expanded gap in phase region (zero coupling)
    l_taper = 0.450     # Smooth transition from coupler gap to phase gap
    l_phase = 2.300     # PCM phase modulation length (pi shift)
    l_fanout = 0.550    # S-bend fanout to output pitch (500 nm)
    y_pitch_out = 0.500
    l_in = 0.350
    l_out = 0.350
    iso_buf = 0.080

    y_sep_c = w_wg + gap_c
    y1_c = y_sep_c / 2.0
    y2_c = -y_sep_c / 2.0

    y_sep_phase = w_wg + gap_phase
    y1_phase = y_sep_phase / 2.0
    y2_phase = -y_sep_phase / 2.0

    y_out1 = y_pitch_out / 2.0
    y_out2 = -y_pitch_out / 2.0

    cell_len = l_in + l_3db + l_taper + l_phase + l_taper + l_3db + l_fanout + l_out
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
    x_c1_end = x_in_end + l_3db
    x_tap1_end = x_c1_end + l_taper
    x_phase_end = x_tap1_end + l_phase
    x_tap2_end = x_phase_end + l_taper
    x_c2_end = x_tap2_end + l_3db
    x_bend_end = x_c2_end + l_fanout

    clad_mat = mp.Medium(index=1.449)
    core_mat = mp.Medium(index=3.480)

    n_pcm = 2.70 if state == "amorphous" else 3.30
    k_pcm = 1.0e-5 if state == "amorphous" else 8.0e-5
    patch_mat = mp.Medium(index=n_pcm, D_conductivity=2*math.pi*fcen*k_pcm)

    # 1. Reference Single Straight Waveguide
    x_mon = x_bend_end + 0.15
    mon_w = w_wg * 2.5
    mon_h = h_wg * 2.5

    sources = [
        mp.EigenModeSource(
            src=mp.GaussianSource(fcen, fwidth=df),
            center=mp.Vector3(x_left + dpml + 0.10, y1_c, h_wg / 2.0),
            size=mp.Vector3(0, w_wg * 3.0, h_wg * 3.0),
            eig_band=1,
            direction=mp.X
        )
    ]

    print(f"[{state.upper()}] 1. Running Reference Straight Waveguide Simulation...")
    geom_ref = [
        mp.Block(center=mp.Vector3(0, 0, -sz/4.0), size=mp.Vector3(mp.inf, mp.inf, sz/2.0), material=clad_mat),
        mp.Block(center=mp.Vector3(0, y1_c, h_wg/2.0), size=mp.Vector3(mp.inf, w_wg, h_wg), material=core_mat)
    ]
    sim_ref = mp.Simulation(
        cell_size=cell, boundary_layers=pml_layers, geometry=geom_ref,
        sources=sources, resolution=resolution, default_material=clad_mat
    )
    flux_ref = sim_ref.add_flux(fcen, 0, 1, mp.FluxRegion(
        center=mp.Vector3(x_mon, y1_c, h_wg/2.0), size=mp.Vector3(0, mon_w, mon_h)
    ))
    sim_ref.run(until_after_sources=mp.stop_when_fields_decayed(20, mp.Hz, mp.Vector3(x_mon, y1_c, h_wg/2.0), 1e-3))
    p_ref = max(mp.get_fluxes(flux_ref)[0], 1e-12)
    print(f"[{state.upper()}] Reference Poynting Flux = {p_ref:.6e}")

    # S-bend and Polygon Utilities
    def make_sbend(x_start, x_end, y_start, y_end, w, n_pts=25):
        dx = x_end - x_start
        dy = y_end - y_start
        u = np.linspace(0, 1, n_pts)
        x = x_start + u * dx
        y = y_start + dy * (3 * u**2 - 2 * u**3)
        top = [mp.Vector3(xi, yi + w / 2.0, 0) for xi, yi in zip(x, y)]
        bot = [mp.Vector3(xi, yi - w / 2.0, 0) for xi, yi in reversed(list(zip(x, y)))]
        return top + bot

    def make_apodized_pcm(x_c, y_c, l_tot, w_p, l_tip=0.45, w_tip=0.030, n_pts=15):
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

    x_pcm_c = (x_tap1_end + x_phase_end) / 2.0
    pcm_pts = make_apodized_pcm(x_pcm_c, y1_phase, l_phase, w_wg)
    cx_pcm = float(np.mean([pt.x for pt in pcm_pts]))
    cy_pcm = float(np.mean([pt.y for pt in pcm_pts]))

    # Construct Complete 1x2 DC-MZI Geometry
    geom = [
        mp.Block(center=mp.Vector3(0, 0, -sz/4.0), size=mp.Vector3(mp.inf, mp.inf, sz/2.0), material=clad_mat),
        # 1. Single Input Lead on WG1 (Port 0)
        mp.Block(center=mp.Vector3((x_left + x_in_end)/2.0, y1_c, h_wg/2.0), size=mp.Vector3(x_in_end - x_left, w_wg, h_wg), material=core_mat),
        # 2. Stage 1: Input 3dB Directional Coupler
        mp.Block(center=mp.Vector3((x_in_end + x_c1_end)/2.0, y1_c, h_wg/2.0), size=mp.Vector3(l_3db, w_wg, h_wg), material=core_mat),
        mp.Block(center=mp.Vector3((x_in_end + x_c1_end)/2.0, y2_c, h_wg/2.0), size=mp.Vector3(l_3db, w_wg, h_wg), material=core_mat),
        # 3. Expansion S-bends to Phase Shifter Pitch (gap = 300 nm)
        mp.Prism(vertices=make_sbend(x_c1_end, x_tap1_end, y1_c, y1_phase, w_wg), height=h_wg, axis=mp.Vector3(0,0,1), material=core_mat),
        mp.Prism(vertices=make_sbend(x_c1_end, x_tap1_end, y2_c, y2_phase, w_wg), height=h_wg, axis=mp.Vector3(0,0,1), material=core_mat),
        # 4. Stage 2: Phase Modulation Region
        mp.Block(center=mp.Vector3(x_pcm_c, y1_phase, h_wg/2.0), size=mp.Vector3(l_phase, w_wg, h_wg), material=core_mat),
        mp.Block(center=mp.Vector3(x_pcm_c, y2_phase, h_wg/2.0), size=mp.Vector3(l_phase, w_wg, h_wg), material=core_mat),
        # Active Sb2S3 on Upper Arm
        mp.Prism(vertices=pcm_pts, height=t_pcm, axis=mp.Vector3(0, 0, 1),
                 center=mp.Vector3(cx_pcm, cy_pcm, h_wg + t_pcm/2.0), material=patch_mat),
        # 5. Convergence S-bends back to Coupler Gap (45 nm)
        mp.Prism(vertices=make_sbend(x_phase_end, x_tap2_end, y1_phase, y1_c, w_wg), height=h_wg, axis=mp.Vector3(0,0,1), material=core_mat),
        mp.Prism(vertices=make_sbend(x_phase_end, x_tap2_end, y2_phase, y2_c, w_wg), height=h_wg, axis=mp.Vector3(0,0,1), material=core_mat),
        # 6. Stage 3: Output 3dB Directional Coupler (Interferometric Combiner)
        mp.Block(center=mp.Vector3((x_tap2_end + x_c2_end)/2.0, y1_c, h_wg/2.0), size=mp.Vector3(l_3db, w_wg, h_wg), material=core_mat),
        mp.Block(center=mp.Vector3((x_tap2_end + x_c2_end)/2.0, y2_c, h_wg/2.0), size=mp.Vector3(l_3db, w_wg, h_wg), material=core_mat),
        # 7. S-Bend Fanout to Output Port Pitch (500 nm)
        mp.Prism(vertices=make_sbend(x_c2_end, x_bend_end, y1_c, y_out1, w_wg), height=h_wg, axis=mp.Vector3(0,0,1), material=core_mat),
        mp.Prism(vertices=make_sbend(x_c2_end, x_bend_end, y2_c, y_out2, w_wg), height=h_wg, axis=mp.Vector3(0,0,1), material=core_mat),
        # 8. Output Leads through Right PML
        mp.Block(center=mp.Vector3((x_bend_end + x_right)/2.0, y_out1, h_wg/2.0), size=mp.Vector3(x_right - x_bend_end, w_wg, h_wg), material=core_mat),
        mp.Block(center=mp.Vector3((x_bend_end + x_right)/2.0, y_out2, h_wg/2.0), size=mp.Vector3(x_right - x_bend_end, w_wg, h_wg), material=core_mat),
    ]

    print(f"[{state.upper()}] 2. Running 1x2 DC-MZI Device Simulation...")
    sim = mp.Simulation(
        cell_size=cell, boundary_layers=pml_layers, geometry=geom,
        sources=sources, resolution=resolution, default_material=clad_mat
    )
    flux_p1 = sim.add_flux(fcen, 0, 1, mp.FluxRegion(center=mp.Vector3(x_mon, y_out1, h_wg/2.0), size=mp.Vector3(0, mon_w, mon_h)))
    flux_p2 = sim.add_flux(fcen, 0, 1, mp.FluxRegion(center=mp.Vector3(x_mon, y_out2, h_wg/2.0), size=mp.Vector3(0, mon_w, mon_h)))

    # In DC-MZI: Amorphous (Delta_phi = 0) -> Port 2 (Cross); Crystalline (Delta_phi = pi) -> Port 1 (Bar)
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
    res = run_1x2_dc_mzi(state=st)
    print("\n" + "=" * 65)
    print(f"RESULTS FOR 1x2 DC-MZI SWITCH ({st.upper()}):")
    print(f"  Footprint Area : {res['area_um2']:.3f} um^2 ({res['cell_len']:.2f} um x {res['cell_w']:.2f} um)")
    print(f"  Port 1 (Bar)   : {res['p1_bar_pct']:.2f}%")
    print(f"  Port 2 (Cross) : {res['p2_cross_pct']:.2f}%")
    print(f"  Total Trans    : {res['tot_trans_pct']:.2f}%")
    print(f"  Insertion Loss : {res['il_db']:.3f} dB")
    print(f"  Crosstalk      : {res['xt_db']:.2f} dB")
    print("=" * 65)
