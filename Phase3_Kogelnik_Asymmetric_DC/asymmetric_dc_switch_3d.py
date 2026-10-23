"""
PROJECT JANUS — ASYMMETRIC DIRECTIONAL COUPLER SWITCH (NEXT CANDIDATE)
=======================================================================
Strategy: SINGLE-WG PCM with ASYMMETRIC WG WIDTHS

Core Physics Change vs. Previous Design
----------------------------------------
Previous design: BOTH waveguides coated with PCM (WG1 switchable, WG2 always amorphous).
  → Amorphous state: Δβ = 0 → full cross-coupling OK.
  → Crystalline state: Δβ = β_cr − β_am, must satisfy Δβ/2κ = √3 (Kogelnik).
  → FAILURE: Δβ in 3D is larger than MPB 2D predicts (top-only PCM asymmetry) →
             ratio Δβ/2κ > √3 consistently → crystalline bar transmission only 70-83%.

This design: ONLY WG1 gets PCM. WG2 has no PCM, but uses a slightly wider silicon
              core (w2 > w1) so that in the AMORPHOUS state the two WGs are phase-matched.

  Tuning lever:  w2 is chosen so n_eff(WG2, bare Si, w2) = n_eff(WG1, Si+PCM_am, w1).
                 Then Δβ_am = 0 → resonant coupling → Port 2 (cross) in amorphous.
                 In crystalline: WG1 index rises sharply (n_pcm_cr=3.30 vs n_pcm_am=2.70)
                 → n_eff(WG1_cr) >> n_eff(WG2) → large Δβ_cr → coupling quenched → Port 1 (bar).

Advantages
----------
1. No ambiguity about WG2 PCM state — it has none.
2. w2 can be precisely tuned (lithography-friendly, DFB-style index trim).
3. Δβ_cr is MAXIMIZED because crystalline pushes WG1 index much higher than WG2.
4. Amorphous resonance (Δβ=0) → L_c = π/(2κ) simple condition.

Design Point (starting values, to be swept)
--------------------------------------------
  w1 = 0.240 um  (WG1 core, with PCM on top)
  w2 = 0.270 um  (WG2 core, wider to match WG1_am effective index — NO PCM)
  gap = 0.085 um  (85 nm evanescent gap; must give κ = Δβ_cr / (2√3))
  L_c = 4.20 um  (straight coupling length; L_eff = L_c + 2×0.35 + 0.30 ≈ 5.20 um)
  t_pcm = 0.150 um (150 nm Sb2S3)

  λ₀ = 1064 nm, Si h=220 nm, SiO₂ cladding
  n_Si = 3.480, n_SiO2 = 1.449
  n_am  = 2.70,  k_am = 1e-5
  n_cr  = 3.30,  k_cr = 8e-5

Sweep Grid (run on Azure VM)
-----------------------------
  w2 ∈ {0.260, 0.265, 0.270, 0.275, 0.280} um
  gap ∈ {0.075, 0.080, 0.085, 0.090} um
  L_c ∈ [3.80, 4.60] um step 0.20 um

Usage
-----
  python asymmetric_dc_switch_3d.py --state both --w2 0.270 --gap 0.085 --lc 4.20 --resolution 20
  python asymmetric_dc_switch_3d.py --sweep --resolution 20
"""

import sys, os, math, time, json, csv, argparse
import numpy as np

try:
    import meep as mp
    HAS_MEEP = True
    _MEEP_ERR = None
except Exception as e:
    HAS_MEEP = False
    _MEEP_ERR = str(e)


# ─────────────────────────────────────────────────────────────────────────────
# Geometry helpers
# ─────────────────────────────────────────────────────────────────────────────

def cosine_sbend_prism(x0, x1, y0, y1, w_in, w_out, n_pts=30):
    """Cosine S-bend with width taper; returns vertex list for mp.Prism."""
    u = np.linspace(0, 1, n_pts)
    x = x0 + (x1 - x0) * u
    y = y0 + (y1 - y0) * 0.5 * (1.0 - np.cos(np.pi * u))
    w = w_in + (w_out - w_in) * 0.5 * (1.0 - np.cos(np.pi * u))
    top = [mp.Vector3(xi, yi + wi / 2.0, 0) for xi, yi, wi in zip(x, y, w)]
    bot = [mp.Vector3(xi, yi - wi / 2.0, 0)
           for xi, yi, wi in reversed(list(zip(x, y, w)))]
    return top + bot


def parabolic_pcm_taper(x0, yc, l_body, l_tip, w_pcm, w_tip=0.020, n_pts=20):
    """Parabolic PCM patch: tip→body→tip strictly inside straight section."""
    u_l = np.linspace(0, 1, n_pts)[:-1]
    xl = x0 + u_l * l_tip
    wl = w_tip + (w_pcm - w_tip) * (2 * u_l - u_l ** 2)
    xs  = x0 + l_tip
    xe  = xs + l_body
    u_r = np.linspace(0, 1, n_pts)[1:]
    xr = xe + u_r * l_tip
    wr = w_tip + (w_pcm - w_tip) * (1.0 - u_r ** 2)
    top = [mp.Vector3(xi, yc + wi / 2.0, 0) for xi, wi in zip(xl, wl)]
    top.append(mp.Vector3(xs, yc + w_pcm / 2.0, 0))
    top.append(mp.Vector3(xe, yc + w_pcm / 2.0, 0))
    top.extend([mp.Vector3(xi, yc + wi / 2.0, 0) for xi, wi in zip(xr, wr)])
    bot = [mp.Vector3(xi, yc - wi / 2.0, 0)
           for xi, wi in reversed(list(zip(xr, wr)))]
    bot.append(mp.Vector3(xe, yc - w_pcm / 2.0, 0))
    bot.append(mp.Vector3(xs, yc - w_pcm / 2.0, 0))
    bot.extend([mp.Vector3(xi, yc - wi / 2.0, 0)
                for xi, wi in reversed(list(zip(xl, wl)))])
    return top + bot


def adiabatic_input_taper(x0, x1, yc, w_full, w_tip=0.020, n_pts=20):
    """Linear adiabatic input taper for WG2 (avoids blunt-facet reflection)."""
    u = np.linspace(0, 1, n_pts)
    x = x0 + (x1 - x0) * u
    w = w_tip + (w_full - w_tip) * u
    top = [mp.Vector3(xi, yc + wi / 2.0, 0) for xi, wi in zip(x, w)]
    bot = [mp.Vector3(xi, yc - wi / 2.0, 0)
           for xi, wi in reversed(list(zip(x, w)))]
    return top + bot


# ─────────────────────────────────────────────────────────────────────────────
# Main FDTD engine
# ─────────────────────────────────────────────────────────────────────────────

def run_asym_dc_switch(
        state      = "amorphous",
        l_c        = 4.20,   # Straight coupling body length (um)
        gap        = 0.085,  # Evanescent gap (um)
        w1         = 0.240,  # WG1 width — carries PCM (um)
        w2         = 0.270,  # WG2 width — NO PCM, wider to phase-match WG1_am (um)
        w_bend     = 0.320,  # Widened S-bend width for low radiation loss (um)
        h_wg       = 0.220,  # Silicon core height (um)
        t_pcm      = 0.150,  # Sb2S3 thickness on WG1 only (um)
        l_tip      = 0.350,  # PCM taper length per side (um)
        w_tip_pcm  = 0.020,  # PCM taper lithographic tip (um)
        l_in       = 0.350,  # Input lead length (um)
        l_bend     = 1.400,  # S-bend fanout length (um) — slightly longer for wider pitch
        y_pitch    = 0.600,  # Output port pitch (um) — wider for w2=270nm
        l_out      = 0.350,  # Output lead length (um)
        iso_buf    = 0.100,  # Lateral cell isolation buffer (um)
        resolution = 20,
        verbose    = True,
):
    """
    Asymmetric Directional Coupler: WG1 has PCM (switchable), WG2 is bare Si (wider).
    Amorphous  → Δβ≈0 → resonant cross coupling → Port 2 ≥ 98%.
    Crystalline → large Δβ → coupling quenched → Port 1 ≥ 98%.
    """
    if not HAS_MEEP:
        err_msg = f"MEEP not found or failed to load (Error: {_MEEP_ERR})."
        err_msg += "\nPlease ensure your conda environment with PyMEEP is activated (e.g., 'conda activate base' or 'conda activate meep', or run via 'python' instead of 'python3')."
        raise ImportError(err_msg)

    t0 = time.time()
    lambda_0 = 1.064
    fcen = 1.0 / lambda_0
    df   = 0.10 * fcen

    # PCM materials (only WG1 gets PCM)
    clad_mat = mp.Medium(index=1.449)
    core_mat = mp.Medium(index=3.480)
    pcm_am   = mp.Medium(index=2.70, D_conductivity=2*math.pi*fcen*1.0e-5)
    pcm_cr   = mp.Medium(index=3.30, D_conductivity=2*math.pi*fcen*8.0e-5)
    mat_wg1  = pcm_cr if state.lower() == "crystalline" else pcm_am

    # Centroids (WG1 = top, source WG; WG2 = bottom, wider bare Si)
    y_sep  = (w1 + w2) / 2.0 + gap   # edge-to-edge gap preserved between centroids
    y1     =  y_sep / 2.0             # WG1 centroid
    y2     = -y_sep / 2.0             # WG2 centroid
    y_out1 =  y_pitch / 2.0
    y_out2 = -y_pitch / 2.0

    # Section lengths
    l_straight = l_c + 2 * l_tip     # full PCM section (body + 2 tapers)
    cell_len   = l_in + l_straight + l_bend + l_out
    cell_w     = y_pitch + 2 * iso_buf
    cell_area  = (l_c + 2 * l_tip) * (w1 + gap + w2)  # active footprint only

    # Simulation box
    dpml = 0.50
    pad  = 0.35
    sx   = cell_len + 2 * (dpml + pad)
    sy   = cell_w   + 2 * (dpml + pad)
    sz   = h_wg + t_pcm + 2 * (dpml + pad)

    cell      = mp.Vector3(sx, sy, sz)
    pml_layers = [mp.PML(dpml)]

    x_left  = -sx / 2.0
    x_right =  sx / 2.0
    x0_cell = -cell_len / 2.0

    x_in_end    = x0_cell  + l_in
    x_str_start = x_in_end
    x_str_end   = x_str_start + l_straight
    x_bend_end  = x_str_end  + l_bend

    # Monitor: aperture centered on output WG, 330×330 nm (clears coupling crosstalk)
    mon_w = max(w1, w2) + 0.090
    mon_h = h_wg + t_pcm + 0.090
    x_mon = x_bend_end + 0.20

    if verbose:
        print(f"\n{'='*72}")
        print(f"ASYMMETRIC DC SWITCH — [{state.upper()}]")
        print(f"  WG1 (PCM): w={w1*1e3:.0f}nm | WG2 (bare Si): w={w2*1e3:.0f}nm")
        print(f"  Gap={gap*1e3:.0f}nm  L_c={l_c:.2f}um  L_straight={l_straight:.2f}um")
        print(f"  PCM on WG1 only: {'CRYSTALLINE' if state.lower()=='crystalline' else 'AMORPHOUS'}")
        print(f"  y1={y1*1e3:.1f}nm  y2={y2*1e3:.1f}nm  (centre-to-centre sep={y_sep*1e3:.1f}nm)")
        print(f"  Active footprint: {cell_area:.4f} um^2")
        print(f"{'='*72}\n")

    # ── EigenMode Source on WG1 ───────────────────────────────────────────────
    sources = [
        mp.EigenModeSource(
            src    = mp.GaussianSource(fcen, fwidth=df),
            center = mp.Vector3(x_left + dpml + 0.10, y1, h_wg / 2.0),
            size   = mp.Vector3(0, max(w1, w2) * 4.0, (h_wg + t_pcm) * 3.0),
            eig_band  = 1,
            direction = mp.X,
        )
    ]

    # ── REFERENCE: straight WG1 (no coupling, no PCM) ─────────────────────────
    if verbose:
        print(f"[{state.upper()}] Step 1: Reference straight waveguide ...")
    geom_ref = [
        mp.Block(center=mp.Vector3(0, 0, -sz/4.0),
                 size=mp.Vector3(mp.inf, mp.inf, sz/2.0), material=clad_mat),
        mp.Block(center=mp.Vector3(0, y1, h_wg/2.0),
                 size=mp.Vector3(mp.inf, w1, h_wg), material=core_mat),
    ]
    sim_ref = mp.Simulation(cell_size=cell, boundary_layers=pml_layers,
                            geometry=geom_ref, sources=sources,
                            resolution=resolution, default_material=clad_mat)
    f_ref = sim_ref.add_flux(fcen, 0, 1,
        mp.FluxRegion(center=mp.Vector3(x_mon, y1, h_wg/2.0),
                      size=mp.Vector3(0, mon_w, mon_h)))
    sim_ref.run(until_after_sources=mp.stop_when_fields_decayed(
        20, mp.Hz, mp.Vector3(x_mon, y1, h_wg/2.0), 1e-3))
    p_ref = max(mp.get_fluxes(f_ref)[0], 1e-12)
    if verbose:
        print(f"[{state.upper()}] Reference flux = {p_ref:.6e}")

    # ── DEVICE geometry ───────────────────────────────────────────────────────
    # PCM taper on WG1 only
    pcm_pts = parabolic_pcm_taper(x_str_start, y1, l_c, l_tip, w1, w_tip_pcm)
    cx_pcm  = float(np.mean([p.x for p in pcm_pts]))
    cy_pcm  = float(np.mean([p.y for p in pcm_pts]))

    # Adiabatic WG2 input taper start (so WG2 grows from a tip into full width)
    l_si_taper   = 0.350
    x_wg2_tip    = x_str_start - l_si_taper

    geom_dev = [
        # Substrate / cladding
        mp.Block(center=mp.Vector3(0, 0, -sz/4.0),
                 size=mp.Vector3(mp.inf, mp.inf, sz/2.0), material=clad_mat),

        # WG1 input lead (full width w1 from left wall to coupling section)
        mp.Block(
            center=mp.Vector3((x_left + x_str_start) / 2.0, y1, h_wg/2.0),
            size  =mp.Vector3(x_str_start - x_left, w1, h_wg),
            material=core_mat),

        # WG2 adiabatic silicon input taper (tip 20 nm → full w2 at coupling start)
        mp.Prism(
            vertices=adiabatic_input_taper(x_wg2_tip, x_str_start, y2, w2, 0.020),
            height=h_wg, axis=mp.Vector3(0, 0, 1),
            center=mp.Vector3((x_wg2_tip + x_str_start) / 2.0, y2, h_wg/2.0),
            material=core_mat),

        # WG1 straight coupling section (core only; PCM added below)
        mp.Block(
            center=mp.Vector3((x_str_start + x_str_end) / 2.0, y1, h_wg/2.0),
            size  =mp.Vector3(l_straight, w1, h_wg),
            material=core_mat),

        # WG2 straight coupling section (wider bare Si — NO PCM)
        mp.Block(
            center=mp.Vector3((x_str_start + x_str_end) / 2.0, y2, h_wg/2.0),
            size  =mp.Vector3(l_straight, w2, h_wg),
            material=core_mat),

        # PCM patch on WG1 only (parabolic taper)
        mp.Prism(
            vertices=pcm_pts, height=t_pcm, axis=mp.Vector3(0, 0, 1),
            center  =mp.Vector3(cx_pcm, cy_pcm, h_wg + t_pcm/2.0),
            material=mat_wg1),

        # WG1 output S-bend (w1 → w_bend, fans out to Port 1 at y_out1)
        mp.Prism(
            vertices=cosine_sbend_prism(x_str_end, x_bend_end, y1, y_out1, w1, w_bend),
            height=h_wg, axis=mp.Vector3(0, 0, 1), material=core_mat),

        # WG2 output S-bend (w2 → w_bend, fans out to Port 2 at y_out2)
        mp.Prism(
            vertices=cosine_sbend_prism(x_str_end, x_bend_end, y2, y_out2, w2, w_bend),
            height=h_wg, axis=mp.Vector3(0, 0, 1), material=core_mat),

        # Output lead Port 1 (WG1 side)
        mp.Block(
            center=mp.Vector3((x_bend_end + x_right) / 2.0, y_out1, h_wg/2.0),
            size  =mp.Vector3(x_right - x_bend_end, w_bend, h_wg),
            material=core_mat),

        # Output lead Port 2 (WG2 side)
        mp.Block(
            center=mp.Vector3((x_bend_end + x_right) / 2.0, y_out2, h_wg/2.0),
            size  =mp.Vector3(x_right - x_bend_end, w_bend, h_wg),
            material=core_mat),
    ]

    # ── Device simulation ─────────────────────────────────────────────────────
    if verbose:
        print(f"[{state.upper()}] Step 2: Device simulation ...")
    sim_dev = mp.Simulation(cell_size=cell, boundary_layers=pml_layers,
                            geometry=geom_dev, sources=sources,
                            resolution=resolution, default_material=clad_mat)
    f1 = sim_dev.add_flux(fcen, 0, 1,
        mp.FluxRegion(center=mp.Vector3(x_mon, y_out1, h_wg/2.0),
                      size=mp.Vector3(0, mon_w, mon_h)))
    f2 = sim_dev.add_flux(fcen, 0, 1,
        mp.FluxRegion(center=mp.Vector3(x_mon, y_out2, h_wg/2.0),
                      size=mp.Vector3(0, mon_w, mon_h)))

    # Decay monitor at the SIGNAL port
    y_dest = y_out2 if state.lower() == "amorphous" else y_out1
    sim_dev.run(until_after_sources=mp.stop_when_fields_decayed(
        20, mp.Hz, mp.Vector3(x_mon, y_dest, h_wg/2.0), 1e-3))

    p1  = max(mp.get_fluxes(f1)[0] / p_ref, 0.0)   # Port 1 (Bar — WG1 side)
    p2  = max(mp.get_fluxes(f2)[0] / p_ref, 0.0)   # Port 2 (Cross — WG2 side)
    tot = p1 + p2

    # Single-port metrics (strict definition: unselected port = loss)
    sig  = p2 if state.lower() == "amorphous" else p1
    leak = p1 if state.lower() == "amorphous" else p2
    il_db = -10.0 * math.log10(max(sig, 1e-9))
    xt_db =  10.0 * math.log10(max(leak / max(sig, 1e-9), 1e-9))
    er_db = -xt_db

    elapsed    = time.time() - t0
    fabric_mm2 = 4096 * cell_area / 1e6

    result = {
        "design"          : "asymmetric_dc",
        "state"           : state,
        "l_c_um"          : round(l_c,  3),
        "gap_nm"          : round(gap * 1e3, 1),
        "w1_nm"           : round(w1  * 1e3, 1),
        "w2_nm"           : round(w2  * 1e3, 1),
        "l_straight_um"   : round(l_straight, 3),
        "active_area_um2" : round(cell_area, 4),
        "fabric_4096_mm2" : round(fabric_mm2, 4),
        "p1_bar_pct"      : round(p1   * 100, 3),
        "p2_cross_pct"    : round(p2   * 100, 3),
        "tot_trans_pct"   : round(tot  * 100, 3),
        "signal_pct"      : round(sig  * 100, 3),
        "leakage_pct"     : round(leak * 100, 3),
        "il_db"           : round(il_db, 4),
        "xt_db"           : round(xt_db, 4),
        "er_db"           : round(er_db, 4),
        "elapsed_s"       : round(elapsed, 2),
        "pass_il"         : il_db  <= 0.088,
        "pass_xt"         : xt_db  <= -20.0,
        "pass_area"       : cell_area < 5.0,
        "all_passed"      : (il_db <= 0.088) and (xt_db <= -20.0) and (cell_area < 5.0),
    }

    if verbose:
        flag_il   = "✓ PASS" if result["pass_il"]   else "✗ FAIL"
        flag_xt   = "✓ PASS" if result["pass_xt"]   else "✗ FAIL"
        flag_area = "✓ PASS" if result["pass_area"] else "✗ FAIL"
        print(f"\n{'='*72}")
        print(f"ASYMMETRIC DC RESULT [{state.upper()}]  L_c={l_c:.2f}um, gap={gap*1e3:.0f}nm,"
              f" w1={w1*1e3:.0f}nm, w2={w2*1e3:.0f}nm")
        print(f"  Port 1 (Bar  / WG1)  : {p1*100:.2f}%")
        print(f"  Port 2 (Cross/ WG2)  : {p2*100:.2f}%")
        print(f"  Total guided trans   : {tot*100:.2f}%")
        print(f"  Signal port          : {sig*100:.2f}%  →  IL = {il_db:.4f} dB  [{flag_il} ≤0.088dB]")
        print(f"  Leakage (unselected) : {leak*100:.2f}%  →  XT = {xt_db:.2f} dB  [{flag_xt} ≤-20dB]")
        print(f"  Active footprint     : {cell_area:.4f} um^2  [{flag_area} <5.0um^2]")
        print(f"  Fabric (4096)        : {fabric_mm2:.4f} mm^2")
        print(f"  VERDICT              : {'*** ALL TARGETS MET! ***' if result['all_passed'] else 'NOT YET — continue sweep'}")
        print(f"  Runtime              : {elapsed:.1f} s")
        print(f"{'='*72}\n")

    return result


# ─────────────────────────────────────────────────────────────────────────────
# Dual-state runner
# ─────────────────────────────────────────────────────────────────────────────

def run_both_states(l_c=4.20, gap=0.085, w2=0.270, resolution=20, **kw):
    print(f"\n{'#'*72}")
    print(f"DUAL-STATE: L_c={l_c:.2f}um, gap={gap*1e3:.0f}nm, w1=240nm, w2={w2*1e3:.0f}nm, res={resolution}")
    print(f"{'#'*72}\n")

    r_am = run_asym_dc_switch("amorphous",   l_c=l_c, gap=gap, w2=w2, resolution=resolution, **kw)
    r_cr = run_asym_dc_switch("crystalline", l_c=l_c, gap=gap, w2=w2, resolution=resolution, **kw)

    print(f"\n{'#'*72}")
    print(f"DUAL-STATE SUMMARY  L_c={l_c:.2f}um gap={gap*1e3:.0f}nm w1=240nm w2={w2*1e3:.0f}nm")
    print(f"{'#'*72}")
    for label, va, vc in [
        ("Signal Port Trans",      f"{r_am['signal_pct']:.2f}%",         f"{r_cr['signal_pct']:.2f}%"),
        ("Leakage (unselected)",   f"{r_am['leakage_pct']:.2f}%",        f"{r_cr['leakage_pct']:.2f}%"),
        ("Total Guided Trans",     f"{r_am['tot_trans_pct']:.2f}%",      f"{r_cr['tot_trans_pct']:.2f}%"),
        ("Insertion Loss",         f"{r_am['il_db']:.4f} dB",            f"{r_cr['il_db']:.4f} dB"),
        ("Crosstalk",              f"{r_am['xt_db']:.2f} dB",            f"{r_cr['xt_db']:.2f} dB"),
        ("Active Footprint",       f"{r_am['active_area_um2']:.4f} um^2",f"{r_cr['active_area_um2']:.4f} um^2"),
        ("IL Pass (≤0.088 dB)",    str(r_am['pass_il']),                 str(r_cr['pass_il'])),
        ("XT Pass (≤−20 dB)",      str(r_am['pass_xt']),                 str(r_cr['pass_xt'])),
        ("ALL TARGETS MET",        str(r_am['all_passed']),              str(r_cr['all_passed'])),
    ]:
        print(f"  {label:<30} {'Amorphous':^20} {'Crystalline':^20}")
        print(f"  {'':30} {va:^20} {vc:^20}")

    out = {"amorphous": r_am, "crystalline": r_cr}
    with open("top_asym_dc_candidates.json", "w") as f:
        json.dump(out, f, indent=2)
    print("\nSaved: top_asym_dc_candidates.json")

    csv_path = "asym_dc_results.csv"
    with open(csv_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(r_am.keys()))
        writer.writeheader()
        writer.writerow(r_am)
        writer.writerow(r_cr)
    print(f"Saved: {csv_path}")
    return r_am, r_cr


# ─────────────────────────────────────────────────────────────────────────────
# Parameter sweep
# ─────────────────────────────────────────────────────────────────────────────

def run_sweep(resolution=20):
    """
    Sweep grid chosen to find both:
      1. Phase-matching condition (Δβ_am≈0): w2 vs gap
      2. Kogelnik null (Δβ_cr/2κ = √3): L_c vs gap

    Grid:
      w2  ∈ {0.255, 0.260, 0.265, 0.270, 0.275, 0.280} um  (index-matching lever)
      gap ∈ {0.075, 0.080, 0.085, 0.090} um               (κ tuning)
      L_c ∈ {3.80, 4.00, 4.20, 4.40, 4.60} um             (L_eff = L_c + 0.70 + 0.30)
    """
    w2_vals  = [0.255, 0.260, 0.265, 0.270, 0.275, 0.280]
    gap_vals = [0.075, 0.080, 0.085, 0.090]
    lc_vals  = [3.80, 4.00, 4.20, 4.40, 4.60]

    results  = []
    total    = len(w2_vals) * len(gap_vals) * len(lc_vals)
    done     = 0

    print(f"\nSWEEP: {total} (w2, gap, L_c) combos × 2 states = {total*2} FDTD runs")
    print(f"Estimated time @ ~90s/run: {total*2*90/3600:.1f} h\n")

    csv_path  = "asym_dc_sweep_results.csv"
    json_path = "top_asym_dc_candidates.json"
    best      = []
    fieldnames = None

    with open(csv_path, "w", newline="") as csvfile:
        writer = None

        for w2 in w2_vals:
            for gap in gap_vals:
                for lc in lc_vals:
                    done += 1
                    tag = f"[{done}/{total}] w2={w2*1e3:.0f}nm gap={gap*1e3:.0f}nm L_c={lc:.2f}um"
                    print(f"\n{tag}")
                    try:
                        r_am = run_asym_dc_switch("amorphous",   l_c=lc, gap=gap, w2=w2,
                                                  resolution=resolution, verbose=False)
                        r_cr = run_asym_dc_switch("crystalline", l_c=lc, gap=gap, w2=w2,
                                                  resolution=resolution, verbose=False)
                    except Exception as e:
                        print(f"  ERROR: {e}")
                        continue

                    row = {
                        "w2_nm"       : int(w2 * 1e3),
                        "gap_nm"      : int(gap * 1e3),
                        "l_c_um"      : lc,
                        "am_signal"   : r_am["signal_pct"],
                        "cr_signal"   : r_cr["signal_pct"],
                        "am_il"       : r_am["il_db"],
                        "cr_il"       : r_cr["il_db"],
                        "am_xt"       : r_am["xt_db"],
                        "cr_xt"       : r_cr["xt_db"],
                        "am_tot"      : r_am["tot_trans_pct"],
                        "cr_tot"      : r_cr["tot_trans_pct"],
                        "active_area" : r_am["active_area_um2"],
                        "am_pass"     : r_am["all_passed"],
                        "cr_pass"     : r_cr["all_passed"],
                        "both_pass"   : r_am["all_passed"] and r_cr["all_passed"],
                    }
                    results.append(row)

                    if writer is None:
                        fieldnames = list(row.keys())
                        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
                        writer.writeheader()
                    writer.writerow(row)
                    csvfile.flush()

                    status = "*** BOTH PASS ***" if row["both_pass"] else (
                        "am PASS only" if row["am_pass"] else (
                        "cr PASS only" if row["cr_pass"] else "--"))
                    print(f"  am={r_am['signal_pct']:.2f}%(IL={r_am['il_db']:.3f}dB XT={r_am['xt_db']:.1f}dB) | "
                          f"cr={r_cr['signal_pct']:.2f}%(IL={r_cr['il_db']:.3f}dB XT={r_cr['xt_db']:.1f}dB) | {status}")

                    if row["both_pass"]:
                        best.append({"amorphous": r_am, "crystalline": r_cr})

    # Save best candidates
    with open(json_path, "w") as f:
        json.dump(best if best else results, f, indent=2)
    print(f"\nSweep complete. Results: {csv_path}  Best candidates: {json_path}")
    return results


# ─────────────────────────────────────────────────────────────────────────────
# CLI
# ─────────────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    ap = argparse.ArgumentParser(
        description="PROJECT JANUS: Asymmetric 1x2 DC Switch (WG1 PCM only, WG2 wider Si)")
    ap.add_argument("--state",      default="both",
                    choices=["amorphous", "crystalline", "both"])
    ap.add_argument("--lc",         type=float, default=4.20,  help="L_c body length (um)")
    ap.add_argument("--gap",        type=float, default=0.085, help="Evanescent gap (um)")
    ap.add_argument("--w1",         type=float, default=0.240, help="WG1 width (um)")
    ap.add_argument("--w2",         type=float, default=0.270, help="WG2 width (um)")
    ap.add_argument("--resolution", type=int,   default=20,    help="FDTD px/um")
    ap.add_argument("--sweep",      action="store_true",       help="Run full parameter sweep")
    args = ap.parse_args()

    if args.sweep:
        run_sweep(resolution=args.resolution)
    elif args.state == "both":
        run_both_states(l_c=args.lc, gap=args.gap, w2=args.w2, resolution=args.resolution,
                        w1=args.w1)
    else:
        run_asym_dc_switch(state=args.state, l_c=args.lc, gap=args.gap,
                           w1=args.w1, w2=args.w2, resolution=args.resolution)
