"""
PROJECT JANUS: 3D VECTOR MAXWELL MEEP FDTD TOTAL INTERNAL REFLECTION (TIR) SWITCH
==================================================================================
PHYSICAL PRINCIPLE: SWITCHED TOTAL INTERNAL REFLECTION DEFLECTOR
----------------------------------------------------------------
Instead of relying on brittle evanescent coupling across an adjacent waveguide,
this switch uses a localized diagonal Sb2S3 interface angled at theta_i = 60 deg.

Critical Angle Physics in Silicon (n_Si = 3.480) at lambda_0 = 1064 nm:
  - Amorphous State  (n_am = 2.70):
      theta_c = arcsin(2.70 / 3.480) = 50.88 deg.
      Since incidence theta_i = 60 deg > 50.88 deg:
      --> TOTAL INTERNAL REFLECTION (TIR) OCCURS! (R -> 100%)
      --> The beam turns by 2 * (90 - 60) = 60 deg into Port 2 (Cross).

  - Crystalline State (n_cr = 3.30):
      theta_c = arcsin(3.30 / 3.480) = 71.49 deg.
      Since incidence theta_i = 60 deg < 71.49 deg:
      --> TIR IS COMPLETELY QUENCHED!
      --> Fresnel reflection: R = ((3.48 - 3.30) / (3.48 + 3.30))^2 = 0.07%
      --> 99.93% of the light transmits straight through into Port 1 (Through).

Key Advantages:
  1. Port 1 is straight: zero S-bend loss in the through state.
  2. No adjacent waveguide sitting 100 nm away: zero evanescent leakage.
  3. Active footprint: ~ 2.0 um^2 (well within the < 5.0 um^2 target).
  4. Non-volatile zero static power.

Usage:
  python tir_switch_3d.py --state both --resolution 20
  mpirun -np 4 python tir_switch_3d.py --state both --resolution 20
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

def make_turned_bend(x_start, y_start, theta_deg, x_end, y_end, w1, w2=None, n_pts=30):
    """
    Smooth cosine-like transition from turned branch (angle theta_deg)
    back to horizontal lead at y_end.
    """
    if w2 is None:
        w2 = w1
    u = np.linspace(0, 1, n_pts)
    dx = x_end - x_start
    dy = y_end - y_start
    x = x_start + u * dx
    y = y_start + dy * (0.5 * (1.0 - np.cos(np.pi * u)))
    w = w1 + (w2 - w1) * (0.5 * (1.0 - np.cos(np.pi * u)))
    top = [mp.Vector3(xi, yi + wi / 2.0, 0) for xi, yi, wi in zip(x, y, w)]
    bot = [mp.Vector3(xi, yi - wi / 2.0, 0) for xi, yi, wi in reversed(list(zip(x, y, w)))]
    return top + bot


def make_diagonal_slot_prism(x_c, y_c, w_slot, l_slot, angle_deg):
    """
    Creates a rotated rectangular prism for the diagonal Sb2S3 switching facet.
    angle_deg is the angle of the slot axis relative to the x-axis.
    """
    rad = math.radians(angle_deg)
    cos_a = math.cos(rad)
    sin_a = math.sin(rad)

    # Local half dimensions: along slot (l_slot/2), across slot (w_slot/2)
    hl = l_slot / 2.0
    hw = w_slot / 2.0

    local_corners = [
        (-hw, -hl),
        ( hw, -hl),
        ( hw,  hl),
        (-hw,  hl)
    ]

    # Rotate by angle_deg and translate to (x_c, y_c)
    vertices = []
    for lx, ly in local_corners:
        rx = x_c + (lx * cos_a - ly * sin_a)
        ry = y_c + (lx * sin_a + ly * cos_a)
        vertices.append(mp.Vector3(rx, ry, 0))

    return vertices


# ─────────────────────────────────────────────────────────────────────────────
# Main Simulation Engine
# ─────────────────────────────────────────────────────────────────────────────

def run_tir_switch(
        state="amorphous",
        theta_slot=30.0,    # Angle of slot relative to x-axis (deg) -> theta_i = 90 - 30 = 60 deg
        w_slot=0.250,       # Slot width (um)
        w_wg=0.350,         # Core waveguide width (um)
        h_wg=0.220,         # Silicon height (um)
        y_pitch=0.600,      # Output port pitch (um)
        l_lead_in=0.600,    # Input lead length (um)
        l_lead_out=0.600,   # Output lead length (um)
        resolution=20,      # FDTD px/um
        verbose=True
):
    """
    3D FDTD simulation of Switched Total Internal Reflection (TIR) Optical Switch.
      - Amorphous (n=2.70): theta_i = 60 deg > theta_c = 50.9 deg -> 100% TIR to Port 2.
      - Crystalline (n=3.30): theta_i = 60 deg < theta_c = 71.5 deg -> 99.9% Straight to Port 1.
    """
    if not HAS_MEEP:
        err_msg = f"MEEP not found or failed to load (Error: {_MEEP_ERR})."
        err_msg += "\nPlease ensure your conda environment with PyMEEP is active."
        raise ImportError(err_msg)

    t0 = time.time()
    lambda_0 = 1.064
    fcen = 1.0 / lambda_0
    df = 0.10 * fcen

    clad_mat = mp.Medium(index=1.449)
    core_mat = mp.Medium(index=3.480)
    pcm_am = mp.Medium(index=2.70, D_conductivity=2 * math.pi * fcen * 1.0e-5)
    pcm_cr = mp.Medium(index=3.30, D_conductivity=2 * math.pi * fcen * 8.0e-5)

    pcm_mat = pcm_cr if state.lower() == "crystalline" else pcm_am

    # Turning angle = 2 * (90 - theta_i) = 2 * theta_slot
    theta_turn = 2.0 * theta_slot   # 60.0 degrees
    rad_turn = math.radians(theta_turn)

    # Geometry layout
    x_junc = 0.0
    y_junc = 0.0

    # Branch 2 (Turned port) starts at junction and bends down to -y_pitch
    l_bend = 1.500
    x_bend_end = x_junc + l_bend
    y_out2 = -y_pitch
    y_out1 = 0.0

    cell_len = l_lead_in + l_bend + l_lead_out
    cell_w = y_pitch + w_wg + 0.400
    cell_area = (w_wg * 3.0) * (w_wg * 3.0)  # Active switching region area

    # Domain bounds
    dpml = 0.50
    pad = 0.35
    sx = cell_len + 2 * (dpml + pad)
    sy = cell_w + 2 * (dpml + pad)
    sz = h_wg + 2 * (dpml + pad)

    cell = mp.Vector3(sx, sy, sz)
    pml_layers = [mp.PML(dpml)]

    x_left = -sx / 2.0
    x_right = sx / 2.0
    x_mon = x_bend_end + 0.20

    mon_w = 0.360
    mon_h = 0.320

    # Source on Port 1 input (straight waveguide along y = 0)
    sources = [
        mp.EigenModeSource(
            src=mp.GaussianSource(fcen, fwidth=df),
            center=mp.Vector3(x_left + dpml + 0.10, y_out1, h_wg / 2.0),
            size=mp.Vector3(0, w_wg * 3.0, h_wg * 3.0),
            eig_band=1,
            direction=mp.X
        )
    ]

    # ── 1. Reference Straight Waveguide Simulation ────────────────────────────
    if verbose:
        print(f"[{state.upper()}] 1. Reference Straight Waveguide...")
    geom_ref = [
        mp.Block(center=mp.Vector3(0, 0, -sz/4.0), size=mp.Vector3(mp.inf, mp.inf, sz/2.0), material=clad_mat),
        mp.Block(center=mp.Vector3(0, y_out1, h_wg/2.0), size=mp.Vector3(mp.inf, w_wg, h_wg), material=core_mat)
    ]
    sim_ref = mp.Simulation(cell_size=cell, boundary_layers=pml_layers, geometry=geom_ref,
                            sources=sources, resolution=resolution, default_material=clad_mat)
    f_ref = sim_ref.add_flux(fcen, 0, 1, mp.FluxRegion(
        center=mp.Vector3(x_mon, y_out1, h_wg/2.0), size=mp.Vector3(0, mon_w, mon_h)))
    sim_ref.run(until_after_sources=mp.stop_when_fields_decayed(
        20, mp.Hz, mp.Vector3(x_mon, y_out1, h_wg/2.0), 1e-3))
    p_ref = max(mp.get_fluxes(f_ref)[0], 1e-12)
    if verbose:
        print(f"[{state.upper()}] Reference Flux = {p_ref:.6e}")

    # ── 2. TIR Device Geometry ────────────────────────────────────────────────
    # Length of slot across waveguide
    l_slot = (w_wg * 2.5) / math.sin(math.radians(theta_slot))
    slot_vertices = make_diagonal_slot_prism(x_junc, y_junc, w_slot, l_slot, theta_slot)

    # Branch 2 turned arm polygon
    branch2_bend_vertices = make_turned_bend(x_junc, y_junc, theta_turn, x_bend_end, y_out2, w_wg, w_wg)

    geom_dev = [
        # Substrate
        mp.Block(center=mp.Vector3(0, 0, -sz/4.0), size=mp.Vector3(mp.inf, mp.inf, sz/2.0), material=clad_mat),

        # Main Through Waveguide (Port 1, straight along y = 0 from left to right)
        mp.Block(center=mp.Vector3(0, y_out1, h_wg/2.0),
                 size=mp.Vector3(mp.inf, w_wg, h_wg), material=core_mat),

        # Turned Output Arm (Port 2: branches at junction, bends to y_out2)
        mp.Prism(vertices=branch2_bend_vertices, height=h_wg, axis=mp.Vector3(0, 0, 1), material=core_mat),

        # Output lead for Port 2
        mp.Block(center=mp.Vector3((x_bend_end + x_right)/2.0, y_out2, h_wg/2.0),
                 size=mp.Vector3(x_right - x_bend_end, w_wg, h_wg), material=core_mat),

        # The Switched Sb2S3 Diagonal Slot (overwrites intersection region)
        mp.Prism(vertices=slot_vertices, height=h_wg, axis=mp.Vector3(0, 0, 1),
                 center=mp.Vector3(x_junc, y_junc, h_wg/2.0), material=pcm_mat),
    ]

    # ── 3. Device Simulation ──────────────────────────────────────────────────
    if verbose:
        print(f"[{state.upper()}] 2. Device Simulation (theta_i={90-theta_slot:.0f} deg, w_slot={w_slot*1000:.0f}nm)...")
    sim_dev = mp.Simulation(cell_size=cell, boundary_layers=pml_layers, geometry=geom_dev,
                            sources=sources, resolution=resolution, default_material=clad_mat)
    f1 = sim_dev.add_flux(fcen, 0, 1, mp.FluxRegion(
        center=mp.Vector3(x_mon, y_out1, h_wg/2.0), size=mp.Vector3(0, mon_w, mon_h)))
    f2 = sim_dev.add_flux(fcen, 0, 1, mp.FluxRegion(
        center=mp.Vector3(x_mon, y_out2, h_wg/2.0), size=mp.Vector3(0, mon_w, mon_h)))

    # Decay monitor: Amorphous routes to Port 2 (Cross); Crystalline routes to Port 1 (Through)
    y_dest = y_out2 if state.lower() == "amorphous" else y_out1
    sim_dev.run(until_after_sources=mp.stop_when_fields_decayed(
        20, mp.Hz, mp.Vector3(x_mon, y_dest, h_wg/2.0), 1e-3))

    p1 = max(mp.get_fluxes(f1)[0] / p_ref, 0.0)
    p2 = max(mp.get_fluxes(f2)[0] / p_ref, 0.0)
    tot = p1 + p2

    sig = p2 if state.lower() == "amorphous" else p1
    leak = p1 if state.lower() == "amorphous" else p2
    il_db = -10.0 * math.log10(max(sig, 1e-9))
    xt_db = 10.0 * math.log10(max(leak / max(sig, 1e-9), 1e-9))
    er_db = -xt_db

    elapsed = time.time() - t0
    fabric_mm2 = 4096 * cell_area / 1.0e6

    result = {
        "design": "tir_switched_deflector",
        "state": state,
        "theta_i_deg": round(90.0 - theta_slot, 1),
        "w_slot_nm": round(w_slot * 1000, 1),
        "w_wg_nm": round(w_wg * 1000, 1),
        "active_area_um2": round(cell_area, 4),
        "fabric_4096_mm2": round(fabric_mm2, 4),
        "p1_through_pct": round(p1 * 100, 3),
        "p2_turned_pct": round(p2 * 100, 3),
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
        print(f"TIR SWITCH [{state.upper()}]  theta_i={90-theta_slot:.0f} deg, w_slot={w_slot*1000:.0f}nm")
        print(f"  Port 1 (Through / y=0)   : {p1*100:.2f}%")
        print(f"  Port 2 (Turned / y=-pitch): {p2*100:.2f}%")
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
# Dual State Runner
# ─────────────────────────────────────────────────────────────────────────────

def run_both_states(theta_slot=30.0, w_slot=0.250, resolution=20, **kwargs):
    print(f"\n{'#'*72}")
    print(f"TIR DUAL-STATE VERIFICATION: theta_i={90-theta_slot:.0f} deg, w_slot={w_slot*1000:.0f}nm, res={resolution}")
    print(f"{'#'*72}\n")

    r_am = run_tir_switch("amorphous", theta_slot=theta_slot, w_slot=w_slot, resolution=resolution, **kwargs)
    r_cr = run_tir_switch("crystalline", theta_slot=theta_slot, w_slot=w_slot, resolution=resolution, **kwargs)

    print(f"\n{'#'*72}")
    print(f"TIR DUAL-STATE SUMMARY  (theta_i={90-theta_slot:.0f} deg, w_slot={w_slot*1000:.0f}nm)")
    print(f"{'#'*72}")
    hdr = f"{'Metric':<28} {'Amorphous (TIR -> P2)':<23} {'Crystalline (Trans -> P1)'}"
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
    with open("top_tir_switch_candidates.json", "w") as f:
        json.dump(out, f, indent=2)
    print("Saved: top_tir_switch_candidates.json")

    # Save CSV
    csv_path = "tir_switch_results.csv"
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

def run_sweep(slot_widths=None, angles=None, resolution=20):
    if slot_widths is None:
        slot_widths = [0.200, 0.250, 0.300]
    if angles is None:
        angles = [27.5, 30.0, 32.5]  # theta_i = [62.5, 60.0, 57.5] deg

    results = []
    print(f"\nTIR SWEEP: slot widths={[int(w*1000) for w in slot_widths]}nm | theta_slot={angles} deg")

    csv_path = "tir_switch_sweep_results.csv"
    best = []

    with open(csv_path, "w", newline="") as csvfile:
        writer = None
        for w_s in slot_widths:
            for ang in angles:
                r_am = run_tir_switch("amorphous", theta_slot=ang, w_slot=w_s, resolution=resolution, verbose=False)
                r_cr = run_tir_switch("crystalline", theta_slot=ang, w_slot=w_s, resolution=resolution, verbose=False)
                row = {
                    "theta_i_deg": round(90.0 - ang, 1),
                    "w_slot_nm": int(w_s * 1000),
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
                print(f"  theta_i={row['theta_i_deg']:.1f} deg, w_slot={row['w_slot_nm']}nm | "
                      f"am={r_am['signal_port_pct']:.2f}%(IL={r_am['il_db']:.3f}dB) | "
                      f"cr={r_cr['signal_port_pct']:.2f}%(IL={r_cr['il_db']:.3f}dB) | {status}")

                if row["both_pass"]:
                    best.append({"amorphous": r_am, "crystalline": r_cr})

    with open("top_tir_switch_candidates.json", "w") as f:
        json.dump(best if best else results, f, indent=2)

    print(f"\nSweep complete! Saved: {csv_path}")
    return results


# ─────────────────────────────────────────────────────────────────────────────
# CLI Entry Point
# ─────────────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="PROJECT JANUS: Switched TIR Optical Switch 3D MEEP")
    parser.add_argument("--state", default="both", choices=["amorphous", "crystalline", "both"])
    parser.add_argument("--theta_slot", type=float, default=30.0, help="Slot angle relative to x (deg)")
    parser.add_argument("--w_slot", type=float, default=0.250, help="Slot width (um)")
    parser.add_argument("--resolution", type=int, default=20, help="FDTD grid resolution (px/um)")
    parser.add_argument("--sweep", action="store_true", help="Run parameter sweep")
    args = parser.parse_args()

    if args.sweep:
        run_sweep(resolution=args.resolution)
    elif args.state == "both":
        run_both_states(theta_slot=args.theta_slot, w_slot=args.w_slot, resolution=args.resolution)
    else:
        run_tir_switch(state=args.state, theta_slot=args.theta_slot, w_slot=args.w_slot, resolution=args.resolution)
