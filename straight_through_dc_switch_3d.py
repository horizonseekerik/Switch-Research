"""
PROJECT JANUS: 3D FULL-WAVE VECTOR MAXWELL MEEP FDTD DIRECTIONAL COUPLER SWITCH
================================================================================
ARCHITECTURE: STRAIGHT-THROUGH PORT 1 (y = 0) + S-BENT DROP PORT 2 (y = -pitch)

Core Physical Advantages:
  1. Port 1 is completely STRAIGHT (y = 0) from input to output:
     - Zero S-bend radiation loss.
     - Zero mode distortion.
     - In the Bar state (Crystalline), light travels straight into Port 1 with
       virtually 0.00 dB bend loss!
  2. Port 2 is the DROP arm:
     - Runs parallel in the interaction zone (gap = 95-105 nm).
     - Only Port 2 takes an S-bend down to y = -y_pitch.
  3. Vertical Mode Alignment:
     - Both WG1 and WG2 have matching Sb2S3 layer stacks (150 nm on top).
     - In Amorphous state: Delta_beta = 0 is EXACT, 3D vertical modes overlap
       perfectly -> 100% resonant cross-coupling into WG2 (Port 2).
     - In Crystalline state: WG1 index jumps (n=3.30 vs n=2.70) -> Delta_beta = sqrt(3)*2*kappa
       cancels evanescent transfer -> light continues straight to Port 1.

Target Specs:
  - Single-Port Transmission >= 98.0% (IL <= 0.088 dB) for BOTH states.
  - Crosstalk (unselected port) <= -20.0 dB (< 1.0%).
  - Active Footprint < 5.0 um^2.

Usage:
  python straight_through_dc_switch_3d.py --state both --lc 3.65 --gap 0.100 --resolution 20
  python straight_through_dc_switch_3d.py --sweep --resolution 20
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
# Geometry Helpers
# ─────────────────────────────────────────────────────────────────────────────

def make_cosine_sbend(x_start, x_end, y_start, y_end, w1, w2=None, n_pts=30):
    """Cosine S-bend with smooth width tapering."""
    if w2 is None:
        w2 = w1
    dx = x_end - x_start
    dy = y_end - y_start
    u = np.linspace(0, 1, n_pts)
    x = x_start + u * dx
    y = y_start + dy * (0.5 * (1.0 - np.cos(np.pi * u)))
    w = w1 + (w2 - w1) * (0.5 * (1.0 - np.cos(np.pi * u)))
    top = [mp.Vector3(xi, yi + wi / 2.0, 0) for xi, yi, wi in zip(x, y, w)]
    bot = [mp.Vector3(xi, yi - wi / 2.0, 0) for xi, yi, wi in reversed(list(zip(x, y, w)))]
    return top + bot


def make_tapered_pcm(x_start, y_c, l_body, l_t, w_p, w_t=0.020, n_pts=20):
    """Adiabatic raised-cosine PCM patch: tip -> body -> tip, strictly inside straight section."""
    u_l = np.linspace(0, 1, n_pts)[:-1]
    xl = x_start + u_l * l_t
    wl = w_t + (w_p - w_t) * 0.5 * (1.0 - np.cos(np.pi * u_l))
    x_body_start = x_start + l_t
    x_body_end = x_body_start + l_body
    u_r = np.linspace(0, 1, n_pts)[1:]
    xr = x_body_end + u_r * l_t
    wr = w_t + (w_p - w_t) * 0.5 * (1.0 + np.cos(np.pi * u_r))
    top = [mp.Vector3(xi, y_c + wi/2.0, 0) for xi, wi in zip(xl, wl)]
    top.append(mp.Vector3(x_body_start, y_c + w_p/2.0, 0))
    top.append(mp.Vector3(x_body_end, y_c + w_p/2.0, 0))
    top.extend([mp.Vector3(xi, y_c + wi/2.0, 0) for xi, wi in zip(xr, wr)])
    bot = [mp.Vector3(xi, y_c - wi/2.0, 0) for xi, wi in reversed(list(zip(xr, wr)))]
    bot.append(mp.Vector3(x_body_end, y_c - w_p/2.0, 0))
    bot.append(mp.Vector3(x_body_start, y_c - w_p/2.0, 0))
    bot.extend([mp.Vector3(xi, y_c - wi/2.0, 0) for xi, wi in reversed(list(zip(xl, wl)))])
    return top + bot


def make_wg_input_taper(x_start, x_end, y_c, w_full, w_tip=0.020, n_pts=20):
    """Adiabatic raised-cosine silicon inverse taper for WG2 entrance."""
    u = np.linspace(0, 1, n_pts)
    x = x_start + u * (x_end - x_start)
    w = w_tip + (w_full - w_tip) * 0.5 * (1.0 - np.cos(np.pi * u))
    top = [mp.Vector3(xi, y_c + wi/2.0, 0) for xi, wi in zip(x, w)]
    bot = [mp.Vector3(xi, y_c - wi/2.0, 0) for xi, wi in reversed(list(zip(x, w)))]
    return top + bot


# ─────────────────────────────────────────────────────────────────────────────
# Main Simulation Engine
# ─────────────────────────────────────────────────────────────────────────────

def run_straight_through_dc(
        state="amorphous",
        l_c=4.60,          # Straight coupling length (um) - centered on optimal crossing
        gap=0.100,         # Evanescent gap (um) - verified low-loss regime
        w_wg=0.240,        # Waveguide core width (um)
        w_bend=0.300,      # Widened low-loss S-bend width for WG2 (um)
        h_wg=0.220,        # Silicon core height (um)
        t_pcm=0.150,       # Sb2S3 PCM thickness (um)
        l_tip=0.350,       # Raised-cosine adiabatic PCM taper length per side (um)
        w_tip=0.020,       # Lithographic tip width (um)
        l_in=0.350,        # Input lead length (um)
        l_bend=1.300,      # Cosine S-bend fanout length for WG2 (um)
        y_pitch=0.560,     # Output port pitch (um)
        l_out=0.350,       # Output lead length (um)
        iso_buf=0.100,     # Lateral cell isolation buffer (um)
        resolution=20,     # FDTD grid resolution (px/um)
        verbose=True
):
    """
    Simulates the Straight-Through Directional Coupler Switch:
      - Port 1 (Through) is centered at y = 0.0 (completely straight).
      - Port 2 (Cross) is the drop arm: starts at y = -(w_wg + gap), bends to -y_pitch.
    """
    if not HAS_MEEP:
        err_msg = f"MEEP not found or failed to load (Error: {_MEEP_ERR})."
        err_msg += "\nPlease ensure your conda environment with PyMEEP is active."
        raise ImportError(err_msg)

    t0 = time.time()
    lambda_0 = 1.064
    fcen = 1.0 / lambda_0
    df = 0.10 * fcen

    # Material properties
    clad_mat = mp.Medium(index=1.449)
    core_mat = mp.Medium(index=3.480)
    pcm_am = mp.Medium(index=2.70, D_conductivity=2 * math.pi * fcen * 1.0e-5)
    pcm_cr = mp.Medium(index=3.30, D_conductivity=2 * math.pi * fcen * 8.0e-5)

    # WG1 (straight bus) has switchable PCM; WG2 (drop arm) has amorphous PCM
    mat_wg1_pcm = pcm_cr if state.lower() == "crystalline" else pcm_am
    mat_wg2_pcm = pcm_am

    # Geometry coordinates
    y1 = 0.0                      # WG1 is completely straight at y = 0
    y2 = -(w_wg + gap)            # WG2 coupling position
    y_out1 = 0.0                  # Port 1 (Bar) output is straight at y = 0
    y_out2 = -y_pitch             # Port 2 (Cross) output drops to -y_pitch

    l_straight = l_c + 2 * l_tip  # Coupled straight section
    cell_len = l_in + l_straight + l_bend + l_out
    cell_w = y_pitch + w_bend + 2 * iso_buf
    cell_area = l_straight * (w_wg * 2 + gap)  # Active footprint

    # Simulation domain size
    dpml = 0.50
    pad = 0.35
    sx = cell_len + 2 * (dpml + pad)
    sy = cell_w + 2 * (dpml + pad)
    sz = h_wg + t_pcm + 2 * (dpml + pad)

    cell = mp.Vector3(sx, sy, sz)
    pml_layers = [mp.PML(dpml)]

    x_left = -sx / 2.0
    x_right = sx / 2.0
    x0_cell = -cell_len / 2.0

    x_in_end = x0_cell + l_in
    x_str_start = x_in_end
    x_str_end = x_str_start + l_straight
    x_bend_end = x_str_end + l_bend

    # Monitor placement
    mon_w = 0.320
    mon_h = 0.320
    x_mon = x_bend_end + 0.15

    # Source on straight WG1 at y = 0
    sources = [
        mp.EigenModeSource(
            src=mp.GaussianSource(fcen, fwidth=df),
            center=mp.Vector3(x_left + dpml + 0.10, y1, h_wg / 2.0),
            size=mp.Vector3(0, w_wg * 3.0, (h_wg + t_pcm) * 3.0),
            eig_band=1,
            direction=mp.X
        )
    ]

    # ── 1. Reference: Straight Waveguide at y = 0 ─────────────────────────────
    if verbose:
        print(f"[{state.upper()}] 1. Straight Reference Waveguide (y=0)...")
    geom_ref = [
        mp.Block(center=mp.Vector3(0, 0, -sz/4.0), size=mp.Vector3(mp.inf, mp.inf, sz/2.0), material=clad_mat),
        mp.Block(center=mp.Vector3(0, y1, h_wg/2.0), size=mp.Vector3(mp.inf, w_wg, h_wg), material=core_mat)
    ]
    sim_ref = mp.Simulation(cell_size=cell, boundary_layers=pml_layers, geometry=geom_ref,
                            sources=sources, resolution=resolution, default_material=clad_mat)
    f_ref = sim_ref.add_flux(fcen, 0, 1, mp.FluxRegion(
        center=mp.Vector3(x_mon, y1, h_wg/2.0), size=mp.Vector3(0, mon_w, mon_h)))
    sim_ref.run(until_after_sources=mp.stop_when_fields_decayed(
        20, mp.Hz, mp.Vector3(x_mon, y1, h_wg/2.0), 1e-3))
    p_ref = max(mp.get_fluxes(f_ref)[0], 1e-12)
    if verbose:
        print(f"[{state.upper()}] Reference Flux = {p_ref:.6e}")

    # ── 2. Device Geometry ────────────────────────────────────────────────────
    pcm_pts1 = make_tapered_pcm(x_str_start, y1, l_c, l_tip, w_wg, w_tip)
    cx1 = float(np.mean([p.x for p in pcm_pts1]))
    cy1 = float(np.mean([p.y for p in pcm_pts1]))

    pcm_pts2 = make_tapered_pcm(x_str_start, y2, l_c, l_tip, w_wg, w_tip)
    cx2 = float(np.mean([p.x for p in pcm_pts2]))
    cy2 = float(np.mean([p.y for p in pcm_pts2]))

    l_si_taper = 0.350
    x_wg2_tip = x_str_start - l_si_taper

    geom_dev = [
        # Substrate
        mp.Block(center=mp.Vector3(0, 0, -sz/4.0), size=mp.Vector3(mp.inf, mp.inf, sz/2.0), material=clad_mat),

        # ── WG1 (STRAIGHT BUS — NO BENDS): Perfectly continuous along y = 0
        mp.Block(center=mp.Vector3(0, y1, h_wg/2.0),
                 size=mp.Vector3(mp.inf, w_wg, h_wg), material=core_mat),

        # ── WG2 (DROP ARM):
        # Adiabatic inverse taper at entrance
        mp.Prism(vertices=make_wg_input_taper(x_wg2_tip, x_str_start, y2, w_wg, w_tip),
                 height=h_wg, axis=mp.Vector3(0, 0, 1),
                 center=mp.Vector3((x_wg2_tip + x_str_start)/2.0, y2, h_wg/2.0), material=core_mat),
        # Straight coupled section
        mp.Block(center=mp.Vector3((x_str_start + x_str_end)/2.0, y2, h_wg/2.0),
                 size=mp.Vector3(l_straight, w_wg, h_wg), material=core_mat),
        # S-bend diverging away to y_out2
        mp.Prism(vertices=make_cosine_sbend(x_str_end, x_bend_end, y2, y_out2, w_wg, w_bend),
                 height=h_wg, axis=mp.Vector3(0, 0, 1), material=core_mat),
        # Output lead at y_out2
        mp.Block(center=mp.Vector3((x_bend_end + x_right)/2.0, y_out2, h_wg/2.0),
                 size=mp.Vector3(x_right - x_bend_end, w_bend, h_wg), material=core_mat),

        # ── PCM Patches (top of WG1 and WG2)
        mp.Prism(vertices=pcm_pts1, height=t_pcm, axis=mp.Vector3(0, 0, 1),
                 center=mp.Vector3(cx1, cy1, h_wg + t_pcm/2.0), material=mat_wg1_pcm),
        mp.Prism(vertices=pcm_pts2, height=t_pcm, axis=mp.Vector3(0, 0, 1),
                 center=mp.Vector3(cx2, cy2, h_wg + t_pcm/2.0), material=mat_wg2_pcm),
    ]

    # ── 3. Device Simulation ──────────────────────────────────────────────────
    if verbose:
        print(f"[{state.upper()}] 2. Device Simulation (L_c={l_c:.2f}um, gap={gap*1000:.0f}nm)...")
    sim_dev = mp.Simulation(cell_size=cell, boundary_layers=pml_layers, geometry=geom_dev,
                            sources=sources, resolution=resolution, default_material=clad_mat)
    f1 = sim_dev.add_flux(fcen, 0, 1, mp.FluxRegion(
        center=mp.Vector3(x_mon, y_out1, h_wg/2.0), size=mp.Vector3(0, mon_w, mon_h)))
    f2 = sim_dev.add_flux(fcen, 0, 1, mp.FluxRegion(
        center=mp.Vector3(x_mon, y_out2, h_wg/2.0), size=mp.Vector3(0, mon_w, mon_h)))

    y_dest = y_out2 if state.lower() == "amorphous" else y_out1
    sim_dev.run(until_after_sources=mp.stop_when_fields_decayed(
        20, mp.Hz, mp.Vector3(x_mon, y_dest, h_wg/2.0), 1e-3))

    p1 = max(mp.get_fluxes(f1)[0] / p_ref, 0.0)
    p2 = max(mp.get_fluxes(f2)[0] / p_ref, 0.0)
    tot = p1 + p2

    # Single-port metrics:
    # Amorphous state -> Target is Port 2 (Cross)
    # Crystalline state -> Target is Port 1 (Straight Through)
    sig = p2 if state.lower() == "amorphous" else p1
    leak = p1 if state.lower() == "amorphous" else p2
    il_db = -10.0 * math.log10(max(sig, 1e-9))
    xt_db = 10.0 * math.log10(max(leak / max(sig, 1e-9), 1e-9))
    er_db = -xt_db

    elapsed = time.time() - t0
    fabric_mm2 = 4096 * cell_area / 1.0e6

    result = {
        "design": "straight_through_dc",
        "state": state,
        "l_c_um": round(l_c, 3),
        "gap_nm": round(gap * 1000, 1),
        "l_straight_um": round(l_straight, 3),
        "active_area_um2": round(cell_area, 4),
        "fabric_4096_mm2": round(fabric_mm2, 4),
        "p1_bar_pct": round(p1 * 100, 3),
        "p2_cross_pct": round(p2 * 100, 3),
        "tot_trans_pct": round(tot * 100, 3),
        "signal_port_pct": round(sig * 100, 3),
        "leakage_port_pct": round(leak * 100, 3),
        "il_db": round(il_db, 4),
        "xt_db": round(xt_db, 4),
        "er_db": round(er_db, 4),
        "elapsed_s": round(elapsed, 2),
        "pass_il": il_db <= 0.088,
        "pass_xt": xt_db <= -20.0,
        "pass_area": cell_area < 5.0,
        "all_passed": (il_db <= 0.088) and (xt_db <= -20.0) and (cell_area < 5.0),
    }

    if verbose:
        flag_il = "PASS" if result["pass_il"] else "FAIL"
        flag_xt = "PASS" if result["pass_xt"] else "FAIL"
        flag_area = "PASS" if result["pass_area"] else "FAIL"
        print(f"\n{'='*72}")
        print(f"STRAIGHT-THROUGH DC [{state.upper()}]  L_c={l_c:.2f}um, gap={gap*1000:.0f}nm")
        print(f"  Port 1 (Through / y=0)   : {p1*100:.2f}%")
        print(f"  Port 2 (Drop / y=-pitch) : {p2*100:.2f}%")
        print(f"  Total Guided Trans       : {tot*100:.2f}%")
        print(f"  Signal Port              : {sig*100:.2f}% -> IL = {il_db:.4f} dB  [{flag_il} <=0.088dB]")
        print(f"  Leakage (Unselected)     : {leak*100:.2f}% -> XT = {xt_db:.2f} dB  [{flag_xt} <=-20dB]")
        print(f"  Active Footprint         : {cell_area:.3f} um^2  [{flag_area} <5.0um^2]")
        print(f"  4096-Switch Fabric       : {fabric_mm2:.4f} mm^2")
        print(f"  RESULT                   : {'*** ALL TARGETS MET! ***' if result['all_passed'] else 'Targets not yet met'}")
        print(f"  Runtime                  : {elapsed:.1f} s")
        print(f"{'='*72}\n")

    return result


# ─────────────────────────────────────────────────────────────────────────────
# Both States Runner
# ─────────────────────────────────────────────────────────────────────────────

def run_both_states(l_c=4.60, gap=0.100, resolution=20, **kwargs):
    print(f"\n{'#'*72}")
    print(f"STRAIGHT-THROUGH DUAL-STATE VERIFICATION: L_c={l_c:.2f}um, gap={gap*1000:.0f}nm, res={resolution}")
    print(f"{'#'*72}\n")

    r_am = run_straight_through_dc("amorphous", l_c=l_c, gap=gap, resolution=resolution, **kwargs)
    r_cr = run_straight_through_dc("crystalline", l_c=l_c, gap=gap, resolution=resolution, **kwargs)

    print(f"\n{'#'*72}")
    print(f"DUAL-STATE SUMMARY  (L_c={l_c:.2f}um, gap={gap*1000:.0f}nm)")
    print(f"{'#'*72}")
    hdr = f"{'Metric':<28} {'Amorphous (State 0)':<23} {'Crystalline (State 1)'}"
    print(hdr)
    print("-" * 75)
    rows = [
        ("Signal Port Trans", f"{r_am['signal_port_pct']:.2f}%", f"{r_cr['signal_port_pct']:.2f}%"),
        ("Leakage (Unselected)", f"{r_am['leakage_port_pct']:.2f}%", f"{r_cr['leakage_port_pct']:.2f}%"),
        ("Total Guided Trans", f"{r_am['tot_trans_pct']:.2f}%", f"{r_cr['tot_trans_pct']:.2f}%"),
        ("Insertion Loss", f"{r_am['il_db']:.4f} dB", f"{r_cr['il_db']:.4f} dB"),
        ("Crosstalk", f"{r_am['xt_db']:.2f} dB", f"{r_cr['xt_db']:.2f} dB"),
        ("Active Footprint", f"{r_am['active_area_um2']:.3f} um^2", f"{r_cr['active_area_um2']:.3f} um^2"),
        ("4096-Fabric Area", f"{r_am['fabric_4096_mm2']:.4f} mm^2", f"{r_cr['fabric_4096_mm2']:.4f} mm^2"),
        ("IL Pass (<= 0.088dB)", str(r_am['pass_il']), str(r_cr['pass_il'])),
        ("XT Pass (<= -20dB)", str(r_am['pass_xt']), str(r_cr['pass_xt'])),
        ("ALL TARGETS MET", str(r_am['all_passed']), str(r_cr['all_passed'])),
    ]
    for label, va, vc in rows:
        print(f"  {label:<26} {va:<23} {vc}")
    print(f"{'#'*72}\n")

    # Save JSON
    out = {"amorphous": r_am, "crystalline": r_cr}
    with open("top_straight_through_dc_candidates.json", "w") as f:
        json.dump(out, f, indent=2)
    print("Saved: top_straight_through_dc_candidates.json")

    # Save CSV
    csv_path = "straight_through_dc_results.csv"
    with open(csv_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(r_am.keys()))
        writer.writeheader()
        writer.writerow(r_am)
        writer.writerow(r_cr)
    print(f"Saved: {csv_path}")
    return r_am, r_cr


# ─────────────────────────────────────────────────────────────────────────────
# Parameter Sweep
# ─────────────────────────────────────────────────────────────────────────────

def run_sweep(lc_min=4.40, lc_max=4.80, lc_step=0.10,
              gap_list=None, resolution=20):
    if gap_list is None:
        gap_list = [0.100]

    results = []
    print(f"\nSWEEP: L_c in [{lc_min:.2f}, {lc_max:.2f}] step {lc_step:.2f} | gaps={[int(g*1000) for g in gap_list]}nm")
    lc_vals = np.arange(lc_min, lc_max + lc_step / 2, lc_step)

    csv_path = "straight_through_dc_sweep_results.csv"
    best = []

    with open(csv_path, "w", newline="") as csvfile:
        writer = None
        for gap in gap_list:
            for lc in lc_vals:
                r_am = run_straight_through_dc("amorphous", l_c=lc, gap=gap, resolution=resolution, verbose=False)
                r_cr = run_straight_through_dc("crystalline", l_c=lc, gap=gap, resolution=resolution, verbose=False)
                row = {
                    "l_c_um": round(lc, 2),
                    "gap_nm": int(gap * 1000),
                    "am_signal": r_am["signal_port_pct"],
                    "cr_signal": r_cr["signal_port_pct"],
                    "am_il": r_am["il_db"],
                    "cr_il": r_cr["il_db"],
                    "am_xt": r_am["xt_db"],
                    "cr_xt": r_cr["xt_db"],
                    "active_area": r_am["active_area_um2"],
                    "am_pass": r_am["all_passed"],
                    "cr_pass": r_cr["all_passed"],
                    "both_pass": r_am["all_passed"] and r_cr["all_passed"],
                }
                results.append(row)
                if writer is None:
                    writer = csv.DictWriter(csvfile, fieldnames=list(row.keys()))
                    writer.writeheader()
                writer.writerow(row)
                csvfile.flush()

                status = "*** BOTH PASS ***" if row["both_pass"] else (
                    "am PASS" if row["am_pass"] else ("cr PASS" if row["cr_pass"] else "--"))
                print(f"  L_c={lc:.2f}um, gap={int(gap*1000)}nm | "
                      f"am={r_am['signal_port_pct']:.2f}%(IL={r_am['il_db']:.3f}dB,XT={r_am['xt_db']:.1f}dB) | "
                      f"cr={r_cr['signal_port_pct']:.2f}%(IL={r_cr['il_db']:.3f}dB,XT={r_cr['xt_db']:.1f}dB) | {status}")

                if row["both_pass"]:
                    best.append({"amorphous": r_am, "crystalline": r_cr})

    with open("top_straight_through_dc_candidates.json", "w") as f:
        json.dump(best if best else results, f, indent=2)

    print(f"\nSweep complete! Saved: {csv_path}")
    return results


# ─────────────────────────────────────────────────────────────────────────────
# CLI Entry Point
# ─────────────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="PROJECT JANUS: Straight-Through 1x2 DC Switch 3D MEEP")
    parser.add_argument("--state", default="both", choices=["amorphous", "crystalline", "both"])
    parser.add_argument("--lc", type=float, default=4.60, help="Straight coupling length (um)")
    parser.add_argument("--gap", type=float, default=0.100, help="Evanescent gap (um)")
    parser.add_argument("--resolution", type=int, default=20, help="FDTD grid resolution (px/um)")
    parser.add_argument("--sweep", action="store_true", help="Run parameter sweep")
    parser.add_argument("--lc_min", type=float, default=4.40, help="Sweep: min L_c (um)")
    parser.add_argument("--lc_max", type=float, default=4.80, help="Sweep: max L_c (um)")
    parser.add_argument("--lc_step", type=float, default=0.10, help="Sweep: L_c step (um)")
    parser.add_argument("--gaps", type=float, nargs="+", default=[0.100], help="Sweep: list of gaps (um)")
    args = parser.parse_args()

    if args.sweep:
        run_sweep(lc_min=args.lc_min, lc_max=args.lc_max, lc_step=args.lc_step,
                  gap_list=args.gaps, resolution=args.resolution)
    elif args.state == "both":
        run_both_states(l_c=args.lc, gap=args.gap, resolution=args.resolution)
    else:
        run_straight_through_dc(state=args.state, l_c=args.lc, gap=args.gap, resolution=args.resolution)
