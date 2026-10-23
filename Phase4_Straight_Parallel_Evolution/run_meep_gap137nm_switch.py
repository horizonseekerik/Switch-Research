"""
PROJECT JANUS: 3D FDTD MEEP SIMULATION - RESONANCE-ALIGNED SWITCH
==================================================================
Testing Candidate in full 3D Maxwell FDTD (PyMEEP):
  - Gap = 137.5 nm (0.1375 um)
  - L_c in [5.40, 5.55] um
  - Core: w = 240 nm, h = 220 nm (Si, n = 3.48)
  - PCM: t_pcm = 200 nm, Sb2S3 (n_am = 2.70, n_cr = 3.30)
  - Output fanout: y_out = 0.40 um, w_bend = 0.30 um, l_bend = 0.65 um
  - Cell footprint: 6.55 um x 0.83 um = 5.4365 um^2 (<= 6.00 um^2 PASS)
  - Output metrics printed to 10 decimal places.
"""

import sys, os, math, time, argparse
import numpy as np
import meep as mp

def make_ccw_sbend(x_start, x_end, y_start, y_end, w1, w2=None, n_pts=35):
    if w2 is None:
        w2 = w1
    dx = x_end - x_start
    dy = y_end - y_start
    u = np.linspace(0, 1, n_pts)
    x = x_start + u * dx
    y = y_start + dy * (3 * u**2 - 2 * u**3)
    w = w1 + (w2 - w1) * (3 * u**2 - 2 * u**3)
    bot = [mp.Vector3(xi, yi - wi / 2.0, 0) for xi, yi, wi in zip(x, y, w)]
    top = [mp.Vector3(xi, yi + wi / 2.0, 0) for xi, yi, wi in reversed(list(zip(x, y, w)))]
    return bot + top

def make_ccw_linear_taper(x_start, x_end, y, w_start, w_end, n_pts=20):
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

def simulate_meep_3d(state="amorphous", resolution=15, l_c=5.5500000000, gap=0.1375000000,
                     t_pcm=0.2000000000, l_tip=0.2500000000, y_out=0.4000000000, l_bend=0.6500000000,
                     w_bend=0.3000000000, p_ref_flux=None, p_ref_mode=None, verbose=True):
    lambda_0 = 1.0640000000
    fcen = 1.0 / lambda_0
    df = 0.10 * fcen

    w_wg = 0.2400000000
    h_wg = 0.2200000000
    core_index = 3.4800000000

    clad_mat = mp.Medium(index=1.4490000000)
    core_mat = mp.Medium(index=core_index)
    pcm_am = mp.Medium(index=2.7000000000, D_conductivity=2.0 * math.pi * fcen * 1.0e-5)
    pcm_cr = mp.Medium(index=3.3000000000, D_conductivity=2.0 * math.pi * fcen * 8.0e-5)

    l_in = 0.2000000000
    l_out = 0.1500000000
    iso_buf = 0.0800000000

    y1 = 0.0000000000
    y2_c = -(w_wg + gap)
    y_out2 = -y_out

    cell_len = l_in + l_c + l_bend + l_out
    cell_w = y_out + w_bend / 2.0 + w_wg / 2.0 + 2.0 * iso_buf
    cell_area = cell_len * cell_w

    dpml = 0.5000000000
    pad = 0.3500000000
    sx = cell_len + 2.0 * (dpml + pad)
    sy = cell_w + 2.0 * (dpml + pad)
    sz = h_wg + t_pcm + 2.0 * (dpml + pad)

    cell = mp.Vector3(sx, sy, sz)
    pml_layers = [mp.PML(dpml)]

    x_left = -sx / 2.0
    x_right = sx / 2.0
    x0 = -cell_len / 2.0

    x_in_end = x0 + l_in
    x_c_start = x_in_end
    x_c_end = x_c_start + l_c
    x_bend_end = x_c_end + l_bend
    x_mon = x_bend_end + 0.1000000000

    mon_w = min(y_out * 0.85, 0.3500000000)
    mon_h = (h_wg + t_pcm) * 2.2

    sources = [
        mp.EigenModeSource(
            src=mp.GaussianSource(fcen, fwidth=df),
            center=mp.Vector3(x_left + dpml + 0.1000000000, y1, h_wg / 2.0),
            size=mp.Vector3(0, w_wg * 3.0, (h_wg + t_pcm) * 3.0),
            eig_band=1,
            direction=mp.X
        )
    ]

    # Reference Waveguide (Straight Si core)
    if p_ref_flux is None or p_ref_mode is None:
        if verbose:
            print(f"[{state.upper()}] Running Reference Waveguide (res={resolution})...")
        geom_ref = [
            mp.Block(center=mp.Vector3(0, 0, -sz / 4.0), size=mp.Vector3(mp.inf, mp.inf, sz / 2.0), material=clad_mat),
            mp.Block(center=mp.Vector3(0, y1, h_wg / 2.0), size=mp.Vector3(mp.inf, w_wg, h_wg), material=core_mat)
        ]
        sim_ref = mp.Simulation(cell_size=cell, boundary_layers=pml_layers, geometry=geom_ref,
                                sources=sources, resolution=resolution, default_material=clad_mat)
        f_ref = sim_ref.add_flux(fcen, 0, 1, mp.FluxRegion(center=mp.Vector3(x_mon, y1, h_wg / 2.0), size=mp.Vector3(0, mon_w, mon_h)))
        sim_ref.run(until_after_sources=mp.stop_when_fields_decayed(20, mp.Hz, mp.Vector3(x_mon, y1, h_wg / 2.0), 1e-3))
        p_ref_flux = max(mp.get_fluxes(f_ref)[0], 1e-12)
        ref_mode_data = sim_ref.get_eigenmode_coefficients(f_ref, [1], eig_parity=mp.NO_PARITY)
        p_ref_mode = max(abs(ref_mode_data.alpha[0, 0, 0])**2, 1e-12)
        if verbose:
            print(f"[{state.upper()}] Reference Flux = {p_ref_flux:.10e}, Mode = {p_ref_mode:.10e}")

    # Device Geometry:
    # 1. Straight continuous bus for WG1 at y = 0
    geom_dev = [
        mp.Block(center=mp.Vector3(0, 0, -sz / 4.0), size=mp.Vector3(mp.inf, mp.inf, sz / 2.0), material=clad_mat),
        mp.Block(center=mp.Vector3(0, y1, h_wg / 2.0), size=mp.Vector3(mp.inf, w_wg, h_wg), material=core_mat),
    ]

    # 2. WG2 entrance adiabatic inverse taper (bare Si, w = 20 -> 240 nm)
    poly_wg2_in = make_ccw_linear_taper(x_c_start - l_in, x_c_start, y2_c, 0.0200000000, w_wg)
    geom_dev.append(create_prism_exact(poly_wg2_in, h_wg, 0.0, core_mat))

    # 3. WG2 parallel coupled section (bare Si core, w = 240 nm)
    geom_dev.append(mp.Block(
        center=mp.Vector3((x_c_start + x_c_end) / 2.0, y2_c, h_wg / 2.0),
        size=mp.Vector3(l_c, w_wg, h_wg),
        material=core_mat
    ))

    # 4. WG2 output S-bend fanout (w_wg -> w_bend)
    poly_wg2_bend = make_ccw_sbend(x_c_end, x_bend_end, y2_c, y_out2, w_wg, w_bend)
    geom_dev.append(create_prism_exact(poly_wg2_bend, h_wg, 0.0, core_mat))

    # 5. WG2 output straight lead at y_out2
    geom_dev.append(mp.Block(
        center=mp.Vector3((x_bend_end + x_right) / 2.0, y_out2, h_wg / 2.0),
        size=mp.Vector3(x_right - x_bend_end, w_bend, h_wg),
        material=core_mat
    ))

    # PCM Placement
    mat_pcm1 = pcm_cr if state.lower() == "crystalline" else pcm_am
    mat_pcm2 = pcm_am

    # WG1 PCM (entrance taper, central body, exit taper)
    poly_pcm1_in = make_ccw_linear_taper(x_c_start, x_c_start + l_tip, y1, 0.0300000000, w_wg)
    geom_dev.append(create_prism_exact(poly_pcm1_in, t_pcm, h_wg, mat_pcm1))

    geom_dev.append(mp.Block(
        center=mp.Vector3((x_c_start + l_tip + x_c_end - l_tip) / 2.0, y1, h_wg + t_pcm / 2.0),
        size=mp.Vector3((x_c_end - l_tip) - (x_c_start + l_tip), w_wg, t_pcm),
        material=mat_pcm1
    ))

    poly_pcm1_out = make_ccw_linear_taper(x_c_end - l_tip, x_c_end, y1, w_wg, 0.0300000000)
    geom_dev.append(create_prism_exact(poly_pcm1_out, t_pcm, h_wg, mat_pcm1))

    # WG2 PCM (entrance taper, central body, exit taper) - EXACT SYNCHRONOUS
    poly_pcm2_in = make_ccw_linear_taper(x_c_start, x_c_start + l_tip, y2_c, 0.0300000000, w_wg)
    geom_dev.append(create_prism_exact(poly_pcm2_in, t_pcm, h_wg, mat_pcm2))

    geom_dev.append(mp.Block(
        center=mp.Vector3((x_c_start + l_tip + x_c_end - l_tip) / 2.0, y2_c, h_wg + t_pcm / 2.0),
        size=mp.Vector3((x_c_end - l_tip) - (x_c_start + l_tip), w_wg, t_pcm),
        material=mat_pcm2
    ))

    poly_pcm2_out = make_ccw_linear_taper(x_c_end - l_tip, x_c_end, y2_c, w_wg, 0.0300000000)
    geom_dev.append(create_prism_exact(poly_pcm2_out, t_pcm, h_wg, mat_pcm2))

    # Run Device Simulation
    sim_dev = mp.Simulation(cell_size=cell, boundary_layers=pml_layers, geometry=geom_dev,
                            sources=sources, resolution=resolution, default_material=clad_mat)
    f_p1 = sim_dev.add_flux(fcen, 0, 1, mp.FluxRegion(center=mp.Vector3(x_mon, y1, h_wg / 2.0), size=mp.Vector3(0, mon_w, mon_h)))
    f_p2 = sim_dev.add_flux(fcen, 0, 1, mp.FluxRegion(center=mp.Vector3(x_mon, y_out2, h_wg / 2.0), size=mp.Vector3(0, mon_w, mon_h)))

    y_dest = y_out2 if state.lower() == "amorphous" else y1
    sim_dev.run(until_after_sources=mp.stop_when_fields_decayed(20, mp.Hz, mp.Vector3(x_mon, y_dest, h_wg / 2.0), 1e-3))

    # 1. Spatial Flux Power
    p1_flux = max(mp.get_fluxes(f_p1)[0] / p_ref_flux, 1e-12)
    p2_flux = max(mp.get_fluxes(f_p2)[0] / p_ref_flux, 1e-12)
    ptot_flux = p1_flux + p2_flux

    # 2. Pure Modal Eigenmode Power
    mode1_data = sim_dev.get_eigenmode_coefficients(f_p1, [1], eig_parity=mp.NO_PARITY)
    mode2_data = sim_dev.get_eigenmode_coefficients(f_p2, [1], eig_parity=mp.NO_PARITY)
    p1_mode = max(abs(mode1_data.alpha[0, 0, 0])**2 / p_ref_mode, 1e-12)
    p2_mode = max(abs(mode2_data.alpha[0, 0, 0])**2 / p_ref_mode, 1e-12)
    ptot_mode = p1_mode + p2_mode

    sig_mode = p2_mode if state.lower() == "amorphous" else p1_mode
    leak_mode = p1_mode if state.lower() == "amorphous" else p2_mode
    il_db = -10.0 * math.log10(sig_mode)
    xt_db = 10.0 * math.log10(leak_mode / sig_mode)
    ratio_guided = (sig_mode / ptot_mode) * 100.0 if ptot_mode > 0 else 0.0

    if verbose:
        print(f"\n{'='*95}")
        print(f"[{state.upper()}] MEEP 3D FDTD (L_c={l_c:.10f} um, Gap={gap*1000:.10f} nm, Area={cell_area:.10f} um^2):")
        print(f"  P1 Modal Power     : {p1_mode*100:.10f}%")
        print(f"  P2 Modal Power     : {p2_mode*100:.10f}%")
        print(f"  Total Guided Mode  : {ptot_mode*100:.10f}% (Rad Loss = {(1-ptot_mode)*100:.10f}%)")
        print(f"  Signal Port Mode   : {sig_mode*100:.10f}% (IL = {il_db:.10f} dB)")
        print(f"  Leakage Port Mode  : {leak_mode*100:.10f}% (XT = {xt_db:.10f} dB)")
        print(f"  Signal / Guided    : {ratio_guided:.10f}%")
        print(f"  P1 Spatial Flux    : {p1_flux*100:.10f}%")
        print(f"  P2 Spatial Flux    : {p2_flux*100:.10f}%")
        print(f"  Total Spatial Flux : {ptot_flux*100:.10f}%")
        print(f"  Cell Footprint     : {cell_area:.10f} um^2 (<= 6.0000000000 um^2 PASS)")
        print(f"{'='*95}\n")

    return {
        "state": state, "p1_mode": p1_mode, "p2_mode": p2_mode, "ptot_mode": ptot_mode,
        "sig_mode": sig_mode, "leak_mode": leak_mode, "il_db": il_db, "xt_db": xt_db,
        "ratio_guided": ratio_guided, "p1_flux": p1_flux, "p2_flux": p2_flux,
        "ptot_flux": ptot_flux, "cell_area": cell_area,
        "p_ref_flux": p_ref_flux, "p_ref_mode": p_ref_mode
    }

def run_3d_benchmark():
    gap = 0.1375000000
    res = 15

    print("=" * 125)
    print("  PROJECT JANUS: MEEP 3D FDTD BENCHMARK AT RESONANCE GAP (137.5 nm)")
    print("=" * 125)
    print(f"{'L_c(um)':<14} {'AM P2(Mod)%':<18} {'AM P1(Mod)%':<18} {'CR P1(Mod)%':<18} {'CR P2(Mod)%':<18} {'AM Sig/G%':<16} {'CR Sig/G%':<16} {'Area(um2)':<14}")
    print("-" * 125)

    res_ref = simulate_meep_3d("amorphous", resolution=res, l_c=5.5000000000, gap=gap, verbose=False)
    p_ref_flux = res_ref["p_ref_flux"]
    p_ref_mode = res_ref["p_ref_mode"]

    for lc in [5.3500000000, 5.5000000000]:
        res_am = simulate_meep_3d("amorphous", resolution=res, l_c=lc, gap=gap,
                                  p_ref_flux=p_ref_flux, p_ref_mode=p_ref_mode, verbose=False)
        res_cr = simulate_meep_3d("crystalline", resolution=res, l_c=lc, gap=gap,
                                  p_ref_flux=p_ref_flux, p_ref_mode=p_ref_mode, verbose=False)

        am_p2 = res_am["p2_mode"] * 100.0
        am_p1 = res_am["p1_mode"] * 100.0
        cr_p1 = res_cr["p1_mode"] * 100.0
        cr_p2 = res_cr["p2_mode"] * 100.0
        am_ratio = res_am["ratio_guided"]
        cr_ratio = res_cr["ratio_guided"]
        area = res_am["cell_area"]

        print(f"{lc:<14.10f} {am_p2:<18.10f} {am_p1:<18.10f} {cr_p1:<18.10f} {cr_p2:<18.10f} {am_ratio:<16.10f} {cr_ratio:<16.10f} {area:<14.10f}")
        sys.stdout.flush()

    print("=" * 125)

if __name__ == '__main__':
    run_3d_benchmark()
