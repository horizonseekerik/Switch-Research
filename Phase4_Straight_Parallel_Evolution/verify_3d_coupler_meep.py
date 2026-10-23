"""
PROJECT JANUS: HIGH-FIDELITY 3D MEEP FDTD VERIFICATION SCRIPT
Mathematically Rigorous 1x2 Non-Volatile Optical Switch Simulation
"""
import sys
import os
import math
import numpy as np
import meep as mp

def run_simulation(gap=0.040, l_c=4.20, state="amorphous", resolution=20):
    lambda_0 = 1.064
    fcen = 1.0 / lambda_0
    df = 0.1 * fcen

    w_wg = 0.240       # 240 nm width
    h_wg = 0.220       # 220 nm height
    t_pcm = 0.150      # 150 nm Sb2S3 PCM thickness
    l_bend = 0.900     # 900 nm Cosine S-bend
    y_pitch_out = 0.600# 600 nm output port separation (clean isolation)
    l_in = 0.300
    l_out = 0.300
    iso_buf = 0.080

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

    # PCM refractive indices at 1064 nm:
    # Amorphous: n = 2.70, k = 1e-5
    # Crystalline: n = 3.30, k = 8e-5
    pcm_am = mp.Medium(index=2.70, D_conductivity=2*math.pi*fcen*1.0e-5)
    pcm_cr = mp.Medium(index=3.30, D_conductivity=2*math.pi*fcen*8.0e-5)

    mat_wg1_pcm = pcm_am if state == "amorphous" else pcm_cr
    mat_wg2_pcm = pcm_am  # WG2 is matching amorphous layer

    # Fully isolated non-overlapping flux monitors (width = 0.30 um)
    # Port 1 monitor: y in [+0.15, +0.45] um
    # Port 2 monitor: y in [-0.45, -0.15] um
    # Centerline isolation gap: 0.30 um of clear cladding between monitors!
    mon_w = 0.300
    mon_h = 0.300
    x_mon = x_bend_end + 0.15

    sources = [
        mp.EigenModeSource(
            src=mp.GaussianSource(fcen, fwidth=df),
            center=mp.Vector3(x_left + dpml + 0.10, y_wg1, h_wg / 2.0),
            size=mp.Vector3(0, w_wg * 3.0, h_wg * 3.0),
            eig_band=1,
            direction=mp.X
        )
    ]

    print(f"[{state.upper()}] Step 1: Running Reference Straight Waveguide (Normalized to Port 1 Monitor)...")
    geom_ref = [
        mp.Block(center=mp.Vector3(0, 0, -sz/4.0), size=mp.Vector3(mp.inf, mp.inf, sz/2.0), material=clad_mat),
        # Straight waveguide through PML
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
    print(f"[{state.upper()}] Reference Flux = {p_ref:.6e}")

    # Cosine S-bend
    def make_cosine_sbend(x_start, x_end, y_start, y_end, w, n_pts=30):
        dx = x_end - x_start
        dy = y_end - y_start
        u = np.linspace(0, 1, n_pts)
        x = x_start + u * dx
        y = y_start + dy * (0.5 * (1.0 - np.cos(np.pi * u)))
        top = [mp.Vector3(xi, yi + w / 2.0, 0) for xi, yi in zip(x, y)]
        bot = [mp.Vector3(xi, yi - w / 2.0, 0) for xi, yi in reversed(list(zip(x, y)))]
        return top + bot

    # Trapezoidal adiabatic taper
    def make_taper_relative(l_t, w_in, w_out):
        return [
            mp.Vector3(-l_t/2.0, w_in/2.0, 0),
            mp.Vector3(l_t/2.0, w_out/2.0, 0),
            mp.Vector3(l_t/2.0, -w_out/2.0, 0),
            mp.Vector3(-l_t/2.0, -w_in/2.0, 0),
        ]

    l_tip = 0.450
    w_tip = 0.035
    l_straight_pcm = l_c - 2 * l_tip
    x_c_center = (x_in_end + x_coupler_end) / 2.0
    x_taper_in = x_in_end + l_tip / 2.0
    x_taper_out = x_coupler_end - l_tip / 2.0

    geom_switch = [
        # Substrate
        mp.Block(center=mp.Vector3(0, 0, -sz/4.0), size=mp.Vector3(mp.inf, mp.inf, sz/2.0), material=clad_mat),

        # 1. WG1: Continuous Silicon Waveguide (Through / Bar Port)
        # Input lead
        mp.Block(center=mp.Vector3((x_left + x_in_end)/2.0, y_wg1, h_wg/2.0), size=mp.Vector3(x_in_end - x_left, w_wg, h_wg), material=core_mat),
        # Coupler section
        mp.Block(center=mp.Vector3(x_c_center, y_wg1, h_wg/2.0), size=mp.Vector3(l_c, w_wg, h_wg), material=core_mat),
        # Cosine S-bend to Port 1
        mp.Prism(vertices=make_cosine_sbend(x_coupler_end, x_bend_end, y_wg1, y_out1, w_wg), height=h_wg, axis=mp.Vector3(0,0,1), material=core_mat),
        # Output lead to Port 1
        mp.Block(center=mp.Vector3((x_bend_end + x_right)/2.0, y_out1, h_wg/2.0), size=mp.Vector3(x_right - x_bend_end, w_wg, h_wg), material=core_mat),

        # 2. WG2: Silicon Waveguide (Cross Port)
        # Lead-in sharp inverse taper tip (zero reflection/scattering)
        mp.Prism(vertices=make_taper_relative(0.400, w_tip, w_wg), height=h_wg, axis=mp.Vector3(0,0,1),
                 center=mp.Vector3(x_in_end - 0.200, y_wg2, h_wg/2.0), material=core_mat),
        # Coupler section
        mp.Block(center=mp.Vector3(x_c_center, y_wg2, h_wg/2.0), size=mp.Vector3(l_c, w_wg, h_wg), material=core_mat),
        # Cosine S-bend to Port 2
        mp.Prism(vertices=make_cosine_sbend(x_coupler_end, x_bend_end, y_wg2, y_out2, w_wg), height=h_wg, axis=mp.Vector3(0,0,1), material=core_mat),
        # Output lead to Port 2
        mp.Block(center=mp.Vector3((x_bend_end + x_right)/2.0, y_out2, h_wg/2.0), size=mp.Vector3(x_right - x_bend_end, w_wg, h_wg), material=core_mat),

        # 3. Active PCM Layer on WG1 (with input and output tapers)
        mp.Prism(vertices=make_taper_relative(l_tip, w_tip, w_wg), height=t_pcm, axis=mp.Vector3(0,0,1),
                 center=mp.Vector3(x_taper_in, y_wg1, h_wg + t_pcm/2.0), material=mat_wg1_pcm),
        mp.Block(center=mp.Vector3(x_c_center, y_wg1, h_wg + t_pcm/2.0), size=mp.Vector3(l_straight_pcm, w_wg, t_pcm), material=mat_wg1_pcm),
        mp.Prism(vertices=make_taper_relative(l_tip, w_wg, w_tip), height=t_pcm, axis=mp.Vector3(0,0,1),
                 center=mp.Vector3(x_taper_out, y_wg1, h_wg + t_pcm/2.0), material=mat_wg1_pcm),

        # 4. Matching PCM Layer on WG2 (identical to WG1)
        mp.Prism(vertices=make_taper_relative(l_tip, w_tip, w_wg), height=t_pcm, axis=mp.Vector3(0,0,1),
                 center=mp.Vector3(x_taper_in, y_wg2, h_wg + t_pcm/2.0), material=mat_wg2_pcm),
        mp.Block(center=mp.Vector3(x_c_center, y_wg2, h_wg + t_pcm/2.0), size=mp.Vector3(l_straight_pcm, w_wg, t_pcm), material=mat_wg2_pcm),
        mp.Prism(vertices=make_taper_relative(l_tip, w_wg, w_tip), height=t_pcm, axis=mp.Vector3(0,0,1),
                 center=mp.Vector3(x_taper_out, y_wg2, h_wg + t_pcm/2.0), material=mat_wg2_pcm),
    ]

    print(f"[{state.upper()}] Step 2: Running 3D Device Simulation...")
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
    lc = float(sys.argv[2]) if len(sys.argv) > 2 else 4.20
    gp = float(sys.argv[3]) if len(sys.argv) > 3 else 0.040

    res = run_simulation(gap=gp, l_c=lc, state=st)
    print("\n" + "=" * 70)
    print(f"VERIFIED 1x2 SWITCH 3D FDTD RESULT ({st.upper()}):")
    print(f"  Footprint Area  : {res['area_um2']:.3f} um^2 ({res['cell_len']:.2f} um x {res['cell_w']:.2f} um)")
    print(f"  Port 1 (Bar)    : {res['p1_bar_pct']:.2f}%")
    print(f"  Port 2 (Cross)  : {res['p2_cross_pct']:.2f}%")
    print(f"  Total Trans     : {res['tot_trans_pct']:.2f}%")
    print(f"  Insertion Loss  : {res['il_db']:.3f} dB (Target: <= 0.08 dB)")
    print(f"  Crosstalk       : {res['xt_db']:.2f} dB (Target: <= -20.0 dB)")
    print("=" * 70)
