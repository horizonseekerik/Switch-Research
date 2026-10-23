"""
PROJECT JANUS: STAGGERED-PCM DIRECTIONAL COUPLER WITH DECOUPLED BUFFER ZONE
============================================================================
KEY INNOVATION (informed by 11 prior failures):
  After the PCM-loaded coupling zone, WG2 continues as bare Si parallel to
  PCM-loaded WG1 for l_buffer before entering the S-bend. In this buffer:
    - WG1: Si + PCM (n_eff ~ 2.15 am / 2.36 cr)
    - WG2: bare Si  (n_eff ~ 1.80)
    - Delta_beta >> 2*kappa => coupling is ZERO
  This eliminates Phi_ends in BOTH states, locking the two peaks to the
  exact same physical length L_c.

ARCHITECTURE GUARDRAILS (from master failure review):
  ✓ Straight through-bus for WG1 (Lesson: 0.18% radiation loss)
  ✓ Matching PCM on BOTH WGs inside coupling zone (Lesson 3)
  ✓ Adiabatic linear tapers on all PCM terminations (Lesson 5)
  ✓ NO PCM on any curved section (Lesson 4)
  ✓ NO abrupt cross-section steps (Lesson 5)
  ✓ Eigenmode decomposition for power extraction (Lesson 10)
  ✓ 3D FDTD verification only (Lesson 6)
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

def simulate_staggered_pcm(state="amorphous", resolution=15, l_c=4.42, gap=0.098,
                           t_pcm=0.200, l_tip=0.250, l_buffer=0.350,
                           l_hold=0.350, y_pitch=0.450, l_bend=1.000,
                           p_ref_flux=None, p_ref_mode=None, verbose=True):
    """
    Staggered-PCM Directional Coupler Switch.
    
    Layout (x-axis):
      [l_in] [--- L_c (both WGs PCM) ---] [l_buffer (WG2 bare)] [l_bend (S-bend)] [l_out]
    
    WG1: Straight through-bus at y=0. PCM covers [x_c_start, x_c_end + l_buffer + l_hold].
    WG2: Parallel at y=-pitch. PCM covers [x_c_start, x_c_end] only.
         Bare Si in buffer zone. Then S-bend fanout at x_bend_start.
    """
    lambda_0 = 1.064
    fcen = 1.0 / lambda_0
    df = 0.10 * fcen

    w_wg = 0.240
    h_wg = 0.220
    core_index = 3.480

    clad_mat = mp.Medium(index=1.449)
    core_mat = mp.Medium(index=core_index)
    pcm_am = mp.Medium(index=2.70, D_conductivity=2 * math.pi * fcen * 1.0e-5)
    pcm_cr = mp.Medium(index=3.30, D_conductivity=2 * math.pi * fcen * 8.0e-5)

    l_in = 0.350
    l_out = 0.350
    iso_buf = 0.080

    y1 = 0.0                      # WG1 straight bus at y = 0
    y2_c = -(w_wg + gap)          # WG2 parallel coupled section
    y_out2 = -y_pitch             # WG2 output port

    # x-coordinates
    # [l_in][--- L_c ---][l_buffer][--- l_bend ---][l_out]
    cell_len = l_in + l_c + l_buffer + l_bend + l_out
    cell_w = y_pitch + w_wg + 2 * iso_buf
    cell_area = cell_len * cell_w

    dpml = 0.50
    pad = 0.35
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
    x_c_end = x_c_start + l_c           # End of PCM coupling zone
    x_buf_end = x_c_end + l_buffer       # End of buffer zone = start of S-bend
    x_bend_end = x_buf_end + l_bend      # End of S-bend
    x_mon = x_bend_end + 0.15            # Monitor position

    mon_w = w_wg * 2.5
    mon_h = (h_wg + t_pcm) * 2.5

    sources = [
        mp.EigenModeSource(
            src=mp.GaussianSource(fcen, fwidth=df),
            center=mp.Vector3(x_left + dpml + 0.10, y1, h_wg / 2.0),
            size=mp.Vector3(0, w_wg * 3.0, (h_wg + t_pcm) * 3.0),
            eig_band=1,
            direction=mp.X
        )
    ]

    # ── Reference straight waveguide ──────────────────────────────────────────
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
            print(f"[{state.upper()}] Reference Flux = {p_ref_flux:.6e}, Mode = {p_ref_mode:.6e}")

    # ── Full Switch Device Geometry ───────────────────────────────────────────
    # 1. Cladding + WG1 straight continuous bus
    geom_dev = [
        mp.Block(center=mp.Vector3(0, 0, -sz / 4.0), size=mp.Vector3(mp.inf, mp.inf, sz / 2.0), material=clad_mat),
        mp.Block(center=mp.Vector3(0, y1, h_wg / 2.0), size=mp.Vector3(mp.inf, w_wg, h_wg), material=core_mat),
    ]

    # 2. WG2 entrance adiabatic inverse taper (bare Si, w = 20 -> 240 nm)
    poly_wg2_in = make_ccw_linear_taper(x_c_start - l_in, x_c_start, y2_c, 0.020, w_wg)
    geom_dev.append(create_prism_exact(poly_wg2_in, h_wg, 0.0, core_mat))

    # 3. WG2 parallel coupled section (bare Si core, w = 240 nm)
    #    Extends from x_c_start all the way to x_buf_end (coupling zone + buffer zone)
    geom_dev.append(mp.Block(
        center=mp.Vector3((x_c_start + x_buf_end) / 2.0, y2_c, h_wg / 2.0),
        size=mp.Vector3(l_c + l_buffer, w_wg, h_wg),
        material=core_mat
    ))

    # 4. WG2 output S-bend fanout (bare Si, starts at x_buf_end)
    poly_wg2_bend = make_ccw_sbend(x_buf_end, x_bend_end, y2_c, y_out2, w_wg, w_wg)
    geom_dev.append(create_prism_exact(poly_wg2_bend, h_wg, 0.0, core_mat))

    # 5. WG2 output straight lead at y_out2
    geom_dev.append(mp.Block(
        center=mp.Vector3((x_bend_end + x_right) / 2.0, y_out2, h_wg / 2.0),
        size=mp.Vector3(x_right - x_bend_end, w_wg, h_wg),
        material=core_mat
    ))

    # ── PCM Placement ────────────────────────────────────────────────────────
    mat_pcm1 = pcm_cr if state.lower() == "crystalline" else pcm_am
    mat_pcm2 = pcm_am  # WG2 PCM is always amorphous

    # ── WG1 PCM: Covers coupling zone + buffer + hold ────────────────────────
    # Full extent: [x_c_start, x_buf_end + l_hold]
    # With adiabatic tapers at both ends
    x_pcm1_start = x_c_start
    x_pcm1_end = x_buf_end + l_hold  # extends through buffer + hold

    # Entrance taper
    poly_pcm1_in = make_ccw_linear_taper(x_pcm1_start, x_pcm1_start + l_tip, y1, 0.030, w_wg)
    geom_dev.append(create_prism_exact(poly_pcm1_in, t_pcm, h_wg, mat_pcm1))

    # Central body
    x_pcm1_body_start = x_pcm1_start + l_tip
    x_pcm1_body_end = x_pcm1_end - l_tip
    if x_pcm1_body_end > x_pcm1_body_start:
        geom_dev.append(mp.Block(
            center=mp.Vector3((x_pcm1_body_start + x_pcm1_body_end) / 2.0, y1, h_wg + t_pcm / 2.0),
            size=mp.Vector3(x_pcm1_body_end - x_pcm1_body_start, w_wg, t_pcm),
            material=mat_pcm1
        ))

    # Exit taper
    poly_pcm1_out = make_ccw_linear_taper(x_pcm1_end - l_tip, x_pcm1_end, y1, w_wg, 0.030)
    geom_dev.append(create_prism_exact(poly_pcm1_out, t_pcm, h_wg, mat_pcm1))

    # ── WG2 PCM: Covers coupling zone ONLY [x_c_start, x_c_end] ─────────────
    # Entrance taper
    poly_pcm2_in = make_ccw_linear_taper(x_c_start, x_c_start + l_tip, y2_c, 0.030, w_wg)
    geom_dev.append(create_prism_exact(poly_pcm2_in, t_pcm, h_wg, mat_pcm2))

    # Central body
    x_pcm2_body_start = x_c_start + l_tip
    x_pcm2_body_end = x_c_end - l_tip
    if x_pcm2_body_end > x_pcm2_body_start:
        geom_dev.append(mp.Block(
            center=mp.Vector3((x_pcm2_body_start + x_pcm2_body_end) / 2.0, y2_c, h_wg + t_pcm / 2.0),
            size=mp.Vector3(x_pcm2_body_end - x_pcm2_body_start, w_wg, t_pcm),
            material=mat_pcm2
        ))

    # Exit taper (CRITICAL: ends at x_c_end, NOT at x_buf_end)
    poly_pcm2_out = make_ccw_linear_taper(x_c_end - l_tip, x_c_end, y2_c, w_wg, 0.030)
    geom_dev.append(create_prism_exact(poly_pcm2_out, t_pcm, h_wg, mat_pcm2))

    # ── Run Device Simulation ────────────────────────────────────────────────
    sim_dev = mp.Simulation(cell_size=cell, boundary_layers=pml_layers, geometry=geom_dev,
                            sources=sources, resolution=resolution, default_material=clad_mat)
    f_p1 = sim_dev.add_flux(fcen, 0, 1, mp.FluxRegion(center=mp.Vector3(x_mon, y1, h_wg / 2.0), size=mp.Vector3(0, mon_w, mon_h)))
    f_p2 = sim_dev.add_flux(fcen, 0, 1, mp.FluxRegion(center=mp.Vector3(x_mon, y_out2, h_wg / 2.0), size=mp.Vector3(0, mon_w, mon_h)))

    y_dest = y_out2 if state.lower() == "amorphous" else y1
    sim_dev.run(until_after_sources=mp.stop_when_fields_decayed(20, mp.Hz, mp.Vector3(x_mon, y_dest, h_wg / 2.0), 1e-3))

    # ── Power Extraction: Spatial Flux ────────────────────────────────────────
    p1_flux = max(mp.get_fluxes(f_p1)[0] / p_ref_flux, 1e-12)
    p2_flux = max(mp.get_fluxes(f_p2)[0] / p_ref_flux, 1e-12)
    ptot_flux = p1_flux + p2_flux

    # ── Power Extraction: Pure TE0 Eigenmode Decomposition ────────────────────
    mode1_data = sim_dev.get_eigenmode_coefficients(f_p1, [1], eig_parity=mp.NO_PARITY)
    mode2_data = sim_dev.get_eigenmode_coefficients(f_p2, [1], eig_parity=mp.NO_PARITY)
    p1_mode = max(abs(mode1_data.alpha[0, 0, 0])**2 / p_ref_mode, 1e-12)
    p2_mode = max(abs(mode2_data.alpha[0, 0, 0])**2 / p_ref_mode, 1e-12)
    ptot_mode = p1_mode + p2_mode

    sig_mode = p2_mode if state.lower() == "amorphous" else p1_mode
    leak_mode = p1_mode if state.lower() == "amorphous" else p2_mode
    il_db = -10.0 * math.log10(max(sig_mode, 1e-12))
    xt_db = 10.0 * math.log10(max(leak_mode, 1e-12) / max(sig_mode, 1e-12))
    ratio_guided = (sig_mode / ptot_mode) * 100.0 if ptot_mode > 0 else 0.0

    if verbose:
        print(f"\n{'='*80}")
        print(f"[{state.upper()}] STAGGERED-PCM DC (L_c={l_c:.3f}um, Gap={gap*1000:.0f}nm, Buf={l_buffer*1000:.0f}nm):")
        print(f"  Modal Power (TE0): P1 = {p1_mode*100:.6f}%, P2 = {p2_mode*100:.6f}%, Total = {ptot_mode*100:.6f}%")
        print(f"  Signal Port Pwr  : {sig_mode*100:.6f}% (IL = {il_db:.4f} dB)")
        print(f"  Leakage Port Pwr : {leak_mode*100:.6f}% (XT = {xt_db:.2f} dB)")
        print(f"  Signal / Guided  : {ratio_guided:.6f}%")
        print(f"  Spatial Flux Pwr : P1 = {p1_flux*100:.4f}%, P2 = {p2_flux*100:.4f}%, Tot = {ptot_flux*100:.4f}%")
        print(f"  Cell Footprint   : {cell_area:.4f} um^2 (L={cell_len:.3f}um, W={cell_w:.3f}um)")
        print(f"  Radiation Loss   : {(1.0-ptot_mode)*100:.6f}% (modal), {(1.0-ptot_flux)*100:.4f}% (flux)")
        print(f"{'='*80}\n")

    return {
        "state": state, "l_c": l_c, "gap": gap, "l_buffer": l_buffer,
        "p1_mode": p1_mode, "p2_mode": p2_mode, "ptot_mode": ptot_mode,
        "sig_mode": sig_mode, "leak_mode": leak_mode, "il_db": il_db, "xt_db": xt_db,
        "ratio_guided": ratio_guided, "p1_flux": p1_flux, "p2_flux": p2_flux,
        "ptot_flux": ptot_flux, "cell_area": cell_area,
        "p_ref_flux": p_ref_flux, "p_ref_mode": p_ref_mode
    }

if __name__ == '__main__':
    # Quick single-point test
    res = simulate_staggered_pcm("amorphous", resolution=15, l_c=4.42, gap=0.098, verbose=True)
    simulate_staggered_pcm("crystalline", resolution=15, l_c=4.42, gap=0.098,
                           p_ref_flux=res["p_ref_flux"], p_ref_mode=res["p_ref_mode"], verbose=True)
