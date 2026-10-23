"""
SYMMETRIC ZERO-PHASE-ERROR 1x2 NON-VOLATILE OPTICAL SWITCH CELL
================================================================
Full 3D Maxwell FDTD Simulation
Architecture:
  - Input: Bare Si on WG1; WG2 enters via S-bend (gap > 300 nm -> 105.51 nm)
  - Coupler: Parallel section of length L_c at Gap = 105.51 nm
  - PCM patches: Strictly contained within [x_c_start, x_c_end] on BOTH waveguides
    - WG1: Input taper (20->240nm) -> uniform body -> output taper (240->20nm)
    - WG2: Input taper (20->240nm) -> uniform body -> output taper (240->20nm)
  - Output: Bare Si on both arms; WG2 fans out via S-bend to y = -y_pitch
  - Result: ZERO differential phase error in S-bends or access leads
"""

import sys, os, math, time
import numpy as np
import meep as mp

def make_ccw_sbend(x_start, x_end, y_start, y_end, w1, w2, n_pts=30):
    dx = x_end - x_start
    dy = y_end - y_start
    u = np.linspace(0, 1, n_pts)
    x = x_start + u * dx
    y = y_start + dy * (3 * u**2 - 2 * u**3)
    w = w1 + (w2 - w1) * (3 * u**2 - 2 * u**3)
    bot = [mp.Vector3(xi, yi - wi / 2.0, 0) for xi, yi, wi in zip(x, y, w)]
    top = [mp.Vector3(xi, yi + wi / 2.0, 0) for xi, yi, wi in reversed(list(zip(x, y, w)))]
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

def simulate_symmetric_switch(state="amorphous", l_c=5.0000000000, gap=0.1055100000,
                              t_pcm=0.1600000000, l_tip=0.8000000000, l_bend=0.8500000000,
                              y_pitch=0.4500000000, resolution=15, p_ref_mode=None, p_ref_flux=None,
                              verbose=True):
    lambda_0 = 1.0640000000
    fcen = 1.0 / lambda_0
    df = 0.10 * fcen

    w_wg = 0.2400000000
    h_wg = 0.2200000000
    core_index = 3.4800000000

    clad_mat = mp.Medium(index=1.4490000000)
    core_mat = mp.Medium(index=core_index)
    pcm_am = mp.Medium(index=2.7000000000)
    pcm_cr = mp.Medium(index=3.3000000000)

    l_in = 0.3000000000
    l_out = 0.3000000000
    iso_buf = 0.0800000000

    y1 = 0.0000000000
    y2_c = -(w_wg + gap)
    y_out2 = -y_pitch

    # Cell structure: l_in + l_bend (input) + l_c (coupler) + l_bend (output) + l_out
    cell_len = l_in + l_bend + l_c + l_bend + l_out
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
    x_c_start = x_in_end + l_bend
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

    # Reference Waveguide
    if p_ref_mode is None or p_ref_flux is None:
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

    # Integrated Switch Device
    geom_dev = [
        mp.Block(center=mp.Vector3(0, 0, -sz / 4.0), size=mp.Vector3(mp.inf, mp.inf, sz / 2.0), material=clad_mat),
        # Waveguide 1 straight silicon core through entire cell
        mp.Block(center=mp.Vector3(0, y1, h_wg / 2.0), size=mp.Vector3(mp.inf, w_wg, h_wg), material=core_mat),
    ]

    # Waveguide 2 silicon core:
    # 1. Entrance S-bend: curves from y = y_out2 in to y = y2_c
    poly_wg2_in_bend = make_ccw_sbend(x_in_end, x_c_start, y_out2, y2_c, w_wg, w_wg)
    geom_dev.append(create_prism_exact(poly_wg2_in_bend, h_wg, 0.0, core_mat))

    # 2. Parallel coupled core
    geom_dev.append(mp.Block(
        center=mp.Vector3((x_c_start + x_c_end) / 2.0, y2_c, h_wg / 2.0),
        size=mp.Vector3(l_c, w_wg, h_wg),
        material=core_mat
    ))

    # 3. Exit S-bend: curves from y = y2_c out to y = y_out2
    poly_wg2_out_bend = make_ccw_sbend(x_c_end, x_bend_end, y2_c, y_out2, w_wg, w_wg)
    geom_dev.append(create_prism_exact(poly_wg2_out_bend, h_wg, 0.0, core_mat))

    # 4. Waveguide 2 output lead to PML
    geom_dev.append(mp.Block(
        center=mp.Vector3((x_bend_end + x_right) / 2.0, y_out2, h_wg / 2.0),
        size=mp.Vector3(x_right - x_bend_end, w_wg, h_wg),
        material=core_mat
    ))

    # Active Phase-Change Material Patches: strictly on [x_c_start, x_c_end]
    mat_pcm1 = pcm_cr if state.lower() == "crystalline" else pcm_am
    mat_pcm2 = pcm_am

    # WG1 PCM patch
    poly_pcm1_in = make_ccw_linear_taper(x_c_start, x_c_start + l_tip, y1, 0.0200000000, w_wg)
    geom_dev.append(create_prism_exact(poly_pcm1_in, t_pcm, h_wg, mat_pcm1))

    l_body = l_c - 2 * l_tip
    if l_body > 0:
        geom_dev.append(mp.Block(
            center=mp.Vector3((x_c_start + l_tip + x_c_end - l_tip) / 2.0, y1, h_wg + t_pcm / 2.0),
            size=mp.Vector3(l_body, w_wg, t_pcm),
            material=mat_pcm1
        ))

    poly_pcm1_out = make_ccw_linear_taper(x_c_end - l_tip, x_c_end, y1, w_wg, 0.0200000000)
    geom_dev.append(create_prism_exact(poly_pcm1_out, t_pcm, h_wg, mat_pcm1))

    # WG2 PCM patch (strictly identical geometry to WG1)
    poly_pcm2_in = make_ccw_linear_taper(x_c_start, x_c_start + l_tip, y2_c, 0.0200000000, w_wg)
    geom_dev.append(create_prism_exact(poly_pcm2_in, t_pcm, h_wg, mat_pcm2))

    if l_body > 0:
        geom_dev.append(mp.Block(
            center=mp.Vector3((x_c_start + l_tip + x_c_end - l_tip) / 2.0, y2_c, h_wg + t_pcm / 2.0),
            size=mp.Vector3(l_body, w_wg, t_pcm),
            material=mat_pcm2
        ))

    poly_pcm2_out = make_ccw_linear_taper(x_c_end - l_tip, x_c_end, y2_c, w_wg, 0.0200000000)
    geom_dev.append(create_prism_exact(poly_pcm2_out, t_pcm, h_wg, mat_pcm2))

    # Simulation Execution
    sim_dev = mp.Simulation(cell_size=cell, boundary_layers=pml_layers, geometry=geom_dev,
                            sources=sources, resolution=resolution, default_material=clad_mat)
    f_p1 = sim_dev.add_flux(fcen, 0, 1, mp.FluxRegion(center=mp.Vector3(x_mon, y1, h_wg / 2.0), size=mp.Vector3(0, mon_w, mon_h)))
    f_p2 = sim_dev.add_flux(fcen, 0, 1, mp.FluxRegion(center=mp.Vector3(x_mon, y_out2, h_wg / 2.0), size=mp.Vector3(0, mon_w, mon_h)))

    y_target = y_out2 if state.lower() == "amorphous" else y1
    sim_dev.run(until_after_sources=mp.stop_when_fields_decayed(20, mp.Hz, mp.Vector3(x_mon, y_target, h_wg / 2.0), 1e-4))

    raw_flux_p1 = max(mp.get_fluxes(f_p1)[0], 0.0)
    raw_flux_p2 = max(mp.get_fluxes(f_p2)[0], 0.0)

    res_p1 = sim_dev.get_eigenmode_coefficients(f_p1, [1], eig_parity=mp.NO_PARITY)
    res_p2 = sim_dev.get_eigenmode_coefficients(f_p2, [1], eig_parity=mp.NO_PARITY)

    c_p1 = abs(res_p1.alpha[0, 0, 0])**2
    c_p2 = abs(res_p2.alpha[0, 0, 0])**2

    t1_mode = c_p1 / p_ref_mode
    t2_mode = c_p2 / p_ref_mode
    tot_mode = t1_mode + t2_mode

    if state.lower() == "amorphous":
        sig_mode = t2_mode
        leak_mode = t1_mode
    else:
        sig_mode = t1_mode
        leak_mode = t2_mode

    il_mode = -10.0 * math.log10(max(sig_mode, 1e-12))
    xt_mode = 10.0 * math.log10(max(leak_mode, 1e-12) / max(sig_mode, 1e-12))

    if verbose:
        print(f"\n[{state.upper()}] L_c = {l_c:.4f} um, Gap = {gap*1e3:.2f} nm, t_pcm = {t_pcm*1e3:.1f} nm:")
        print(f"  P1 (Bar)  : {t1_mode*100:15.10f}%")
        print(f"  P2 (Cross): {t2_mode*100:15.10f}%")
        print(f"  Guided Tot: {tot_mode*100:15.10f}%")
        print(f"  Signal IL : {il_mode:15.10f} dB")
        print(f"  Modal XT  : {xt_mode:15.10f} dB")
        print(f"  Footprint : {cell_area:15.10f} um^2 (<= 6.0 um^2: {'PASS' if cell_area <= 6.0 else 'FAIL'})")
        sys.stdout.flush()

    return {
        "state": state, "l_c": l_c, "P1": t1_mode, "P2": t2_mode, "tot": tot_mode,
        "sig": sig_mode, "leak": leak_mode, "il": il_mode, "xt": xt_mode,
        "area": cell_area, "p_ref_mode": p_ref_mode, "p_ref_flux": p_ref_flux
    }

if __name__ == "__main__":
    t0 = time.time()
    print("=" * 105)
    print("SYMMETRIC ZERO-PHASE-ERROR SWITCH: 3D FDTD SWEEP")
    print("=" * 105)

    ref_run = simulate_symmetric_switch(state="amorphous", l_c=4.80, verbose=False)
    p_ref_mode = ref_run["p_ref_mode"]
    p_ref_flux = ref_run["p_ref_flux"]

    print(f"{'L_c (um)':<10} {'State 0 P2 %':<16} {'State 0 P1 %':<16} {'State 1 P1 %':<16} {'State 1 P2 %':<16} {'S0 Guided %':<14} {'S1 Guided %':<14}")
    print("-" * 105)

    for lc in [4.40, 4.70, 5.00, 5.30]:
        r0 = simulate_symmetric_switch(state="amorphous", l_c=lc, p_ref_mode=p_ref_mode, p_ref_flux=p_ref_flux, verbose=False)
        r1 = simulate_symmetric_switch(state="crystalline", l_c=lc, p_ref_mode=p_ref_mode, p_ref_flux=p_ref_flux, verbose=False)
        p2_0 = r0["P2"] * 100.0
        p1_0 = r0["P1"] * 100.0
        p1_1 = r1["P1"] * 100.0
        p2_1 = r1["P2"] * 100.0
        tot_0 = r0["tot"] * 100.0
        tot_1 = r1["tot"] * 100.0
        print(f"{lc:<10.2f} {p2_0:<16.6f} {p1_0:<16.6f} {p1_1:<16.6f} {p2_1:<16.6f} {tot_0:<14.6f} {tot_1:<14.6f}")
        sys.stdout.flush()

    print("=" * 105)
    print(f"Total scan execution time = {time.time() - t0:.2f} s")
