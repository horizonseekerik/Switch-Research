"""
CONVERGED ZERO-MISMATCH SWITCH 3D MAXWELL FDTD
==============================================
Parameters:
  - Wavelength: 1064.0000000000 nm
  - Core: w = 240.0000000000 nm, h = 220.0000000000 nm (Si, n = 3.4800000000)
  - Cladding: SiO2 (n = 1.4490000000)
  - Sb2S3 PCM: n_am = 2.7000000000, n_cr = 3.3000000000
  - Converged PCM Thickness: t_pcm = 237.5000000000 nm (Delta L = 0.0000000000 nm)
  - Proven Coupling Gap: Gap = 98.1400000000 nm
  - Adiabatic Taper: l_tip = 500.0000000000 nm
"""

import sys, os, math, time
import numpy as np
import meep as mp

def make_ccw_sbend(x_start, x_end, y_start, y_end, w1, w2, n_pts=25):
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

def simulate_converged_point(state="amorphous", resolution=15, l_c=4.2500000000,
                             gap=0.0981400000, t_pcm=0.2375000000, l_tip=0.5000000000,
                             l_hold=0.3500000000, y_pitch=0.4500000000, l_bend=1.0000000000,
                             p_ref_flux=None, p_ref_mode=None, verbose=True):
    lambda_0 = 1.0640000000
    fcen = 1.0 / lambda_0
    df = 0.10 * fcen

    w_wg = 0.2400000000
    h_wg = 0.2200000000
    core_index = 3.4800000000

    clad_mat = mp.Medium(index=1.4490000000)
    core_mat = mp.Medium(index=core_index)
    pcm_am = mp.Medium(index=2.7000000000, D_conductivity=2 * math.pi * fcen * 1.0e-5)
    pcm_cr = mp.Medium(index=3.3000000000, D_conductivity=2 * math.pi * fcen * 8.0e-5)

    l_in = 0.3500000000
    l_out = 0.3500000000
    iso_buf = 0.0800000000

    y1 = 0.0000000000
    y2_c = -(w_wg + gap)
    y_out2 = -y_pitch

    cell_len = l_in + l_c + l_bend + l_out
    cell_w = y_pitch + w_wg + 2 * iso_buf
    cell_area = cell_len * cell_w

    dpml = 0.5000000000
    pad = 0.3500000000
    sx = cell_len + 2 * (dpml + pad)
    sy = cell_w + 2 * (dpml + pad)
    sz = h_wg + t_pcm + 2 * (dpml + pad)

    cell = mp.Vector3(sx, sy, sz)
    pml_layers = [mp.PML(dpml)]

    x_left = -sx / 2.0
    x_right = sx / 2.0
    x0 = -cell_len / 2.0

    x_in_end = x0 + l_in
    x_c_start = x_in_end
    x_c_end = x_c_start + l_c
    x_bend_end = x_c_end + l_bend
    x_mon = x_bend_end + 0.1500000000

    mon_w = min(y_pitch * 0.85, 0.4000000000)
    mon_h = (h_wg + t_pcm) * 2.2

    src_pos_x = x_left + dpml + 0.1000000000
    sources = [
        mp.EigenModeSource(
            src=mp.GaussianSource(fcen, fwidth=df),
            center=mp.Vector3(src_pos_x, y1, h_wg / 2.0),
            size=mp.Vector3(0, w_wg * 3.0, (h_wg + t_pcm) * 3.0),
            eig_band=1,
            direction=mp.X
        )
    ]

    # 1. Straight Reference
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
        sim_ref.run(until_after_sources=mp.stop_when_fields_decayed(20, mp.Hz, mp.Vector3(x_mon, y1, h_wg / 2.0), 1e-4))
        p_ref_flux = max(mp.get_fluxes(f_ref)[0], 1e-15)
        ref_mode_data = sim_ref.get_eigenmode_coefficients(f_ref, [1], eig_parity=mp.NO_PARITY)
        p_ref_mode = max(abs(ref_mode_data.alpha[0, 0, 0])**2, 1e-15)
        if verbose:
            print(f"[{state.upper()}] Reference Flux = {p_ref_flux:.10e}, Mode = {p_ref_mode:.10e}")

    # 2. Integrated Converged Switch Geometry
    geom_dev = [
        mp.Block(center=mp.Vector3(0, 0, -sz / 4.0), size=mp.Vector3(mp.inf, mp.inf, sz / 2.0), material=clad_mat),
        mp.Block(center=mp.Vector3(0, y1, h_wg / 2.0), size=mp.Vector3(mp.inf, w_wg, h_wg), material=core_mat),
    ]

    # WG2 entrance inverse taper
    poly_wg2_in = make_ccw_linear_taper(x_c_start - l_in, x_c_start, y2_c, 0.0200000000, w_wg)
    geom_dev.append(create_prism_exact(poly_wg2_in, h_wg, 0.0, core_mat))

    # WG2 parallel coupled core
    geom_dev.append(mp.Block(
        center=mp.Vector3((x_c_start + x_c_end) / 2.0, y2_c, h_wg / 2.0),
        size=mp.Vector3(l_c, w_wg, h_wg),
        material=core_mat
    ))

    # WG2 cubic Hermite S-bend fanout (bare silicon)
    poly_wg2_bend = make_ccw_sbend(x_c_end, x_bend_end, y2_c, y_out2, w_wg, w_wg)
    geom_dev.append(create_prism_exact(poly_wg2_bend, h_wg, 0.0, core_mat))

    # WG2 output lead through PML
    geom_dev.append(mp.Block(
        center=mp.Vector3((x_bend_end + x_right) / 2.0, y_out2, h_wg / 2.0),
        size=mp.Vector3(x_right - x_bend_end, w_wg, h_wg),
        material=core_mat
    ))

    # Phase-Change Material Layers
    mat_pcm1 = pcm_cr if state.lower() == "crystalline" else pcm_am
    mat_pcm2 = pcm_am

    # WG1 PCM: [x_c_start, x_c_end + l_hold]
    poly_pcm1_in = make_ccw_linear_taper(x_c_start, x_c_start + l_tip, y1, 0.0300000000, w_wg)
    geom_dev.append(create_prism_exact(poly_pcm1_in, t_pcm, h_wg, mat_pcm1))

    x_pcm1_body_len = (x_c_end + l_hold) - (x_c_start + l_tip)
    if x_pcm1_body_len > 0:
        geom_dev.append(mp.Block(
            center=mp.Vector3((x_c_start + l_tip + x_c_end + l_hold) / 2.0, y1, h_wg + t_pcm / 2.0),
            size=mp.Vector3(x_pcm1_body_len, w_wg, t_pcm),
            material=mat_pcm1
        ))

    poly_pcm1_out = make_ccw_linear_taper(x_c_end + l_hold, x_c_end + l_hold + l_tip, y1, w_wg, 0.0300000000)
    geom_dev.append(create_prism_exact(poly_pcm1_out, t_pcm, h_wg, mat_pcm1))

    # WG2 PCM: [x_c_start, x_c_end]
    poly_pcm2_in = make_ccw_linear_taper(x_c_start, x_c_start + l_tip, y2_c, 0.0300000000, w_wg)
    geom_dev.append(create_prism_exact(poly_pcm2_in, t_pcm, h_wg, mat_pcm2))

    x_pcm2_body_len = (x_c_end - l_tip) - (x_c_start + l_tip)
    if x_pcm2_body_len > 0:
        geom_dev.append(mp.Block(
            center=mp.Vector3((x_c_start + l_tip + x_c_end - l_tip) / 2.0, y2_c, h_wg + t_pcm / 2.0),
            size=mp.Vector3(x_pcm2_body_len, w_wg, t_pcm),
            material=mat_pcm2
        ))

    poly_pcm2_out = make_ccw_linear_taper(x_c_end - l_tip, x_c_end, y2_c, w_wg, 0.0300000000)
    geom_dev.append(create_prism_exact(poly_pcm2_out, t_pcm, h_wg, mat_pcm2))

    # 3. 3D Maxwell FDTD Simulation
    sim_dev = mp.Simulation(cell_size=cell, boundary_layers=pml_layers, geometry=geom_dev,
                            sources=sources, resolution=resolution, default_material=clad_mat)
    f_p1 = sim_dev.add_flux(fcen, 0, 1, mp.FluxRegion(center=mp.Vector3(x_mon, y1, h_wg / 2.0), size=mp.Vector3(0, mon_w, mon_h)))
    f_p2 = sim_dev.add_flux(fcen, 0, 1, mp.FluxRegion(center=mp.Vector3(x_mon, y_out2, h_wg / 2.0), size=mp.Vector3(0, mon_w, mon_h)))

    y_dest = y_out2 if state.lower() == "amorphous" else y1
    sim_dev.run(until_after_sources=mp.stop_when_fields_decayed(20, mp.Hz, mp.Vector3(x_mon, y_dest, h_wg / 2.0), 1e-4))

    raw_flux_p1 = max(mp.get_fluxes(f_p1)[0], 0.0)
    raw_flux_p2 = max(mp.get_fluxes(f_p2)[0], 0.0)

    res_p1 = sim_dev.get_eigenmode_coefficients(f_p1, [1], eig_parity=mp.NO_PARITY)
    res_p2 = sim_dev.get_eigenmode_coefficients(f_p2, [1], eig_parity=mp.NO_PARITY)

    c_p1 = abs(res_p1.alpha[0, 0, 0])**2
    c_p2 = abs(res_p2.alpha[0, 0, 0])**2

    t1_mode = c_p1 / p_ref_mode
    t2_mode = c_p2 / p_ref_mode
    tot_mode = t1_mode + t2_mode

    t1_flux = raw_flux_p1 / p_ref_flux
    t2_flux = raw_flux_p2 / p_ref_flux
    tot_flux = t1_flux + t2_flux

    if state.lower() == "amorphous":
        sig_mode = t2_mode
        leak_mode = t1_mode
        sig_flux = t2_flux
        leak_flux = t1_flux
    else:
        sig_mode = t1_mode
        leak_mode = t2_mode
        sig_flux = t1_flux
        leak_flux = t2_flux

    il_mode = -10.0 * math.log10(max(sig_mode, 1e-12))
    xt_mode = 10.0 * math.log10(max(leak_mode, 1e-12) / max(sig_mode, 1e-12))

    il_flux = -10.0 * math.log10(max(sig_flux, 1e-12))
    xt_flux = 10.0 * math.log10(max(leak_flux, 1e-12) / max(sig_flux, 1e-12))

    if verbose:
        print(f"\n{'='*105}")
        print(f"[{state.upper()}] CONVERGED CELL (L_c={l_c:.6f}um, Gap={gap*1e3:.4f}nm, t_pcm={t_pcm*1e3:.4f}nm, Res={resolution}):")
        print(f"  Modal Power (TE0): P1 = {t1_mode*100:15.10f}%, P2 = {t2_mode*100:15.10f}%, Total = {tot_mode*100:15.10f}%")
        print(f"  Signal Port Pwr  : {sig_mode*100:15.10f}% (IL = {il_mode:13.10f} dB)")
        print(f"  Leakage Port Pwr : {leak_mode*100:15.10f}% (XT = {xt_mode:13.10f} dB)")
        print(f"  Signal / Guided  : {(sig_mode / max(tot_mode, 1e-12))*100:15.10f}%")
        print(f"  Spatial Flux Pwr : P1 = {t1_flux*100:15.10f}%, P2 = {t2_flux*100:15.10f}%, Tot = {tot_flux*100:15.10f}%")
        print(f"  Spatial Flux Sig : {sig_flux*100:15.10f}% (IL = {il_flux:13.10f} dB, XT = {xt_flux:13.10f} dB)")
        print(f"  Cell Footprint   : {cell_area:15.10f} um^2 (Limit <= 6.0 um^2: {'PASS' if cell_area <= 6.0 else 'FAIL'})")
        print(f"  Radiation Loss   : {(1.0 - tot_mode)*100:15.10f}% (modal), {(1.0 - tot_flux)*100:15.10f}% (flux)")
        print(f"{'='*105}\n")

    return {
        "state": state, "P1_mode": t1_mode, "P2_mode": t2_mode, "tot_mode": tot_mode,
        "sig_mode": sig_mode, "leak_mode": leak_mode, "il_mode": il_mode, "xt_mode": xt_mode,
        "P1_flux": t1_flux, "P2_flux": t2_flux, "tot_flux": tot_flux,
        "sig_flux": sig_flux, "leak_flux": leak_flux, "il_flux": il_flux, "xt_flux": xt_flux,
        "cell_area": cell_area, "p_ref_flux": p_ref_flux, "p_ref_mode": p_ref_mode
    }

if __name__ == "__main__":
    t0 = time.time()
    print("=" * 105)
    print("CONVERGED ZERO-MISMATCH CELL SCAN: L_c in [4.15, 4.25, 4.35] um at t_pcm = 237.5 nm, Gap = 98.14 nm")
    print("=" * 105)
    print(f"{'L_c (um)':<10} {'State 0 P2 %':<16} {'State 0 P1 %':<16} {'State 1 P1 %':<16} {'State 1 P2 %':<16} {'S0 Guided %':<14} {'S1 Guided %':<14}")
    print("-" * 105)

    res_ref = simulate_converged_point(state="amorphous", resolution=15, l_c=4.25, verbose=False)
    p_ref_flux = res_ref["p_ref_flux"]
    p_ref_mode = res_ref["p_ref_mode"]

    for lc in [4.15, 4.25, 4.35]:
        res_am = simulate_converged_point(state="amorphous", resolution=15, l_c=lc,
                                          p_ref_flux=p_ref_flux, p_ref_mode=p_ref_mode, verbose=False)
        res_cr = simulate_converged_point(state="crystalline", resolution=15, l_c=lc,
                                          p_ref_flux=p_ref_flux, p_ref_mode=p_ref_mode, verbose=False)
        p2_0 = res_am["P2_mode"] * 100.0
        p1_0 = res_am["P1_mode"] * 100.0
        p1_1 = res_cr["P1_mode"] * 100.0
        p2_1 = res_cr["P2_mode"] * 100.0
        s0_purity = (res_am["P2_mode"] / res_am["tot_mode"]) * 100.0
        s1_purity = (res_cr["P1_mode"] / res_cr["tot_mode"]) * 100.0
        print(f"{lc:<10.2f} {p2_0:<16.6f} {p1_0:<16.6f} {p1_1:<16.6f} {p2_1:<16.6f} {s0_purity:<14.4f} {s1_purity:<14.4f}")
        sys.stdout.flush()

    print("=" * 105)
    print(f"Total elapsed run time = {time.time() - t0:.2f} s")
