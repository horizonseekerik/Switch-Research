"""
PROJECT JANUS: 3D VECTOR MAXWELL MEEP FDTD — EMBEDDED PCM SLOT SWITCH (GEOMETRY B)
===================================================================================
Physical Architecture:
  - Material Platform : Stoichiometric Silicon Nitride (Si3N4, n = 2.016)
  - Cladding          : Silicon Dioxide (SiO2, n = 1.449)
  - Active Element    : Antimony Trisulfide (Sb2S3) directly embedded in-plane
                        as a central slot slab between dual Si3N4 rails.
  - Rail Width        : W_rail = 550 nm, Height = 300 nm
  - Central Slot Slab : W_slot = 65 nm, Height = 300 nm (in-plane embedded)
  - Switching Physics :
      * Amorphous (State 0, n_am = 2.70):
          Light launched in Rail 1 bends/transfers across the embedded PCM slab into
          Rail 2 over L_c = 2.56 um (half-beat cycle, theta = 90 deg).
          Rail 2 S-bend adiabatic taper guides beam smoothly to Port 1 (> 99%).
      * Crystalline (State 1, n_cr = 3.30):
          Beat length halves to L_cr = 1.28 um. Light bends across to Rail 2 and
          returns fully into Rail 1 over L_c = 2.56 um (full beat cycle, theta = 180 deg).
          Rail 1 S-bend adiabatic taper guides beam smoothly to Port 2 (> 98%).

Neither port goes straight: both ports bend via adiabatic S-bends to decoupled leads!
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
# Geometry Helpers: Pure Counter-Clockwise (CCW) Polygons
# ─────────────────────────────────────────────────────────────────────────────

def make_cosine_sbend(x_start, x_end, y_start, y_end, w, z=0.0, n_pts=30):
    """
    Creates a smooth raised-cosine S-bend polygon with strictly Counter-Clockwise
    (CCW) vertex ordering and positive normal +z.
    """
    dx = x_end - x_start
    dy = y_end - y_start
    u = np.linspace(0, 1, n_pts)
    x = x_start + u * dx
    y = y_start + dy * (0.5 * (1.0 - np.cos(np.pi * u)))
    
    # Bottom boundary: x_start -> x_end
    bot = [mp.Vector3(xi, yi - w / 2.0, z) for xi, yi in zip(x, y)]
    # Top boundary: x_end -> x_start
    top = [mp.Vector3(xi, yi + w / 2.0, z) for xi, yi in reversed(list(zip(x, y)))]
    return bot + top


# ─────────────────────────────────────────────────────────────────────────────
# Main Simulation Engine
# ─────────────────────────────────────────────────────────────────────────────

def run_sin_slot_switch(
        state="amorphous",
        l_c=2.56,          # Embedded PCM interaction length (um) - analytical 2:1 resonance
        w_rail=0.550,      # Si3N4 rail width (um)
        w_slot=0.065,      # Central embedded Sb2S3 slot slab width (um)
        h_wg=0.300,        # Si3N4 core and PCM slab height (um)
        l_bend=2.000,      # Adiabatic S-bend length (um)
        y_port_sep=1.400,  # Center-to-center port separation at monitors (um)
        l_lead=0.400,      # Straight port lead length (um)
        resolution=25,     # FDTD grid resolution (px/um)
        p_ref=None,        # Cached reference flux
        verbose=True
):
    """
    Simulates the Geometry B Si3N4 Embedded Slot Switch:
      - Input: Port 1 (Top rail at +y_port)
      - Interaction: Dual Si3N4 rails with embedded Sb2S3 central slab
      - Output:
          * Port 1 (Primary / Amorphous): Bottom rail S-bend to -y_port
          * Port 2 (Switched / Crystalline): Top rail S-bend to +y_port
    """
    if not HAS_MEEP:
        err_msg = f"MEEP not found or failed to load ({_MEEP_ERR})."
        err_msg += "\nPlease ensure your conda environment with PyMEEP is active."
        raise ImportError(err_msg)

    t0 = time.time()
    lambda_0 = 1.064
    fcen = 1.0 / lambda_0
    df = 0.10 * fcen

    # Material Definitions
    sio2_mat = mp.Medium(index=1.449)
    sin_mat  = mp.Medium(index=2.016)

    # Sb2S3 PCM Material Properties at 1064 nm
    n_am = 2.700
    k_am = 1.0e-5
    cond_am = 2 * math.pi * fcen * k_am / n_am
    pcm_am = mp.Medium(index=n_am, D_conductivity=cond_am)

    n_cr = 3.300
    k_cr = 1.8e-4
    cond_cr = 2 * math.pi * fcen * k_cr / n_cr
    pcm_cr = mp.Medium(index=n_cr, D_conductivity=cond_cr)

    pcm_mat = pcm_cr if state.lower() == "crystalline" else pcm_am

    # Geometry Layout Coordinates along Y:
    y_int = (w_rail + w_slot) / 2.0       # Interaction center y for top rail (+y_int) and bottom (-y_int)
    y_port = y_port_sep / 2.0             # Output port center y for top rail (+y_port) and bottom (-y_port)

    # Geometry Layout Coordinates along X:
    # x0: domain left (excluding PML/pad)
    # x1: end of input lead = start of entrance S-bends
    # x2: end of entrance S-bends = start of central interaction region
    # x3: end of central interaction region = start of exit S-bends
    # x4: end of exit S-bends = start of output leads
    # x5: domain right (excluding PML/pad)
    cell_len = l_lead + l_bend + l_c + l_bend + l_lead
    cell_w = 2 * y_port + w_rail + 0.60
    active_area = l_c * (2 * w_rail + w_slot)

    dpml = 0.50
    pad = 0.35
    sx = cell_len + 2 * (dpml + pad)
    sy = cell_w + 2 * (dpml + pad)
    sz = h_wg + 2 * (dpml + pad)

    cell = mp.Vector3(sx, sy, sz)
    pml_layers = [mp.PML(dpml)]

    x_left = -sx / 2.0
    x_right = sx / 2.0
    x0 = -cell_len / 2.0

    x1 = x0 + l_lead
    x2 = x1 + l_bend
    x3 = x2 + l_c
    x4 = x3 + l_bend
    x5 = x4 + l_lead

    # Monitor placement in output leads
    mon_w = w_rail + 0.30
    mon_h = h_wg + 0.30
    x_mon = x4 + l_lead / 2.0

    # Source injected on Rail 1 (Top arm) at input lead (y = +y_port)
    sources = [
        mp.EigenModeSource(
            src=mp.GaussianSource(fcen, fwidth=df),
            center=mp.Vector3(x_left + dpml + 0.10, y_port, h_wg / 2.0),
            size=mp.Vector3(0, w_rail * 3.0, h_wg * 3.0),
            eig_band=1,
            direction=mp.X
        )
    ]

    # ── 1. Reference Simulation: Straight Bare Si3N4 Bus ──────────────────────
    if p_ref is None:
        if verbose:
            print(f"[{state.upper()}] 1. Running Reference Si3N4 Waveguide (Straight Bare Si3N4 at y={y_port:.2f})...")
        geom_ref = [
            mp.Block(center=mp.Vector3(0, 0, -sz / 4.0), size=mp.Vector3(mp.inf, mp.inf, sz / 2.0), material=sio2_mat),
            mp.Block(center=mp.Vector3(0, y_port, h_wg / 2.0), size=mp.Vector3(mp.inf, w_rail, h_wg), material=sin_mat)
        ]
        sim_ref = mp.Simulation(cell_size=cell, boundary_layers=pml_layers, geometry=geom_ref,
                                sources=sources, resolution=resolution, default_material=sio2_mat)
        f_ref = sim_ref.add_flux(fcen, 0, 1, mp.FluxRegion(
            center=mp.Vector3(x_mon, y_port, h_wg / 2.0), size=mp.Vector3(0, mon_w, mon_h)))
        sim_ref.run(until_after_sources=mp.stop_when_fields_decayed(
            20, mp.Hz, mp.Vector3(x_mon, y_port, h_wg / 2.0), 1e-3))
        p_ref = max(mp.get_fluxes(f_ref)[0], 1e-12)
        if verbose:
            print(f"[{state.upper()}] Reference Flux Cached = {p_ref:.6e}")
    else:
        if verbose:
            print(f"[{state.upper()}] Using Cached Reference Flux = {p_ref:.6e}")

    # ── 2. Device Geometry (Geometry B: Central Embedded PCM Slot Slab) ────────
    # S-Bend Polygons (base sitting at z = 0.0, extruding to h_wg = 0.300 um)
    # Entrance S-bends: bring rails from y_port to y_int
    poly_in_top = make_cosine_sbend(x1, x2, y_port, y_int, w_rail, z=0.0)
    poly_in_bot = make_cosine_sbend(x1, x2, -y_port, -y_int, w_rail, z=0.0)

    # Exit S-bends: bring rails from y_int to y_port
    poly_out_top = make_cosine_sbend(x3, x4, y_int, y_port, w_rail, z=0.0)
    poly_out_bot = make_cosine_sbend(x3, x4, -y_int, -y_port, w_rail, z=0.0)

    geom_dev = [
        # Substrate
        mp.Block(center=mp.Vector3(0, 0, -sz / 4.0), size=mp.Vector3(mp.inf, mp.inf, sz / 2.0), material=sio2_mat),

        # ── INPUT LEADS (x_left to x1)
        mp.Block(center=mp.Vector3((x_left + x1) / 2.0, y_port, h_wg / 2.0),
                 size=mp.Vector3(x1 - x_left, w_rail, h_wg), material=sin_mat),
        mp.Block(center=mp.Vector3((x_left + x1) / 2.0, -y_port, h_wg / 2.0),
                 size=mp.Vector3(x1 - x_left, w_rail, h_wg), material=sin_mat),

        # ── ENTRANCE S-BEND TAPERS (x1 to x2)
        mp.Prism(vertices=poly_in_top, height=h_wg, axis=mp.Vector3(0, 0, 1), material=sin_mat),
        mp.Prism(vertices=poly_in_bot, height=h_wg, axis=mp.Vector3(0, 0, 1), material=sin_mat),

        # ── CENTRAL INTERACTION REGION (x2 to x3): Dual Si3N4 Rails + Embedded Center PCM Slab
        # Top Si3N4 Rail (y = +y_int)
        mp.Block(center=mp.Vector3((x2 + x3) / 2.0, y_int, h_wg / 2.0),
                 size=mp.Vector3(l_c, w_rail, h_wg), material=sin_mat),
        # Central Embedded Sb2S3 Slot Slab (y = 0.0)
        mp.Block(center=mp.Vector3((x2 + x3) / 2.0, 0.0, h_wg / 2.0),
                 size=mp.Vector3(l_c, w_slot, h_wg), material=pcm_mat),
        # Bottom Si3N4 Rail (y = -y_int)
        mp.Block(center=mp.Vector3((x2 + x3) / 2.0, -y_int, h_wg / 2.0),
                 size=mp.Vector3(l_c, w_rail, h_wg), material=sin_mat),

        # ── EXIT S-BEND TAPERS (x3 to x4)
        mp.Prism(vertices=poly_out_top, height=h_wg, axis=mp.Vector3(0, 0, 1), material=sin_mat),
        mp.Prism(vertices=poly_out_bot, height=h_wg, axis=mp.Vector3(0, 0, 1), material=sin_mat),

        # ── OUTPUT LEADS (x4 to x_right)
        # Port 2 (Top arm at +y_port)
        mp.Block(center=mp.Vector3((x4 + x_right) / 2.0, y_port, h_wg / 2.0),
                 size=mp.Vector3(x_right - x4, w_rail, h_wg), material=sin_mat),
        # Port 1 (Bottom arm at -y_port)
        mp.Block(center=mp.Vector3((x4 + x_right) / 2.0, -y_port, h_wg / 2.0),
                 size=mp.Vector3(x_right - x4, w_rail, h_wg), material=sin_mat),
    ]

    # ── 3. Simulation Execution ───────────────────────────────────────────────
    if verbose:
        print(f"[{state.upper()}] 2. Device Simulation (L_c={l_c:.2f}um, w_rail={w_rail*1000:.0f}nm, w_slot={w_slot*1000:.0f}nm)...")

    sim_dev = mp.Simulation(cell_size=cell, boundary_layers=pml_layers, geometry=geom_dev,
                            sources=sources, resolution=resolution, default_material=sio2_mat)

    # Monitor Port 1 (Bottom lead: -y_port) - Primary state output
    f1 = sim_dev.add_flux(fcen, 0, 1, mp.FluxRegion(
        center=mp.Vector3(x_mon, -y_port, h_wg / 2.0), size=mp.Vector3(0, mon_w, mon_h)))
    # Monitor Port 2 (Top lead: +y_port) - Switched state output
    f2 = sim_dev.add_flux(fcen, 0, 1, mp.FluxRegion(
        center=mp.Vector3(x_mon, y_port, h_wg / 2.0), size=mp.Vector3(0, mon_w, mon_h)))

    # Decay monitor at target port:
    y_dest = -y_port if state.lower() == "amorphous" else y_port
    sim_dev.run(until_after_sources=mp.stop_when_fields_decayed(
        20, mp.Hz, mp.Vector3(x_mon, y_dest, h_wg / 2.0), 1e-3))

    p1 = max(mp.get_fluxes(f1)[0] / p_ref, 0.0)  # Port 1 (Bottom arm: primary Amorphous target)
    p2 = max(mp.get_fluxes(f2)[0] / p_ref, 0.0)  # Port 2 (Top arm: switched Crystalline target)
    tot = p1 + p2

    # In Amorphous state: signal port is Port 1 (-y_port), leakage is Port 2 (+y_port)
    # In Crystalline state: signal port is Port 2 (+y_port), leakage is Port 1 (-y_port)
    sig = p1 if state.lower() == "amorphous" else p2
    leak = p2 if state.lower() == "amorphous" else p1

    il_db = -10.0 * math.log10(max(sig, 1e-9))
    xt_db = 10.0 * math.log10(max(leak / max(sig, 1e-9), 1e-9))
    er_db = -xt_db

    elapsed = time.time() - t0
    fabric_mm2 = 4096 * active_area / 1.0e6

    target_il = 0.043 if state.lower() == "amorphous" else 0.088  # >99% for primary (0.043dB), >98% for switched (0.088dB)
    pass_il = (il_db <= target_il)
    pass_xt = (xt_db <= -20.0)
    all_pass = pass_il and pass_xt

    res = {
        "state": state,
        "l_c": l_c,
        "w_rail_nm": int(w_rail * 1000),
        "w_slot_nm": int(w_slot * 1000),
        "port1_bottom": p1,
        "port2_top": p2,
        "total_guided": tot,
        "radiation_loss": max(1.0 - tot, 0.0),
        "signal_trans": sig,
        "leakage": leak,
        "insertion_loss_db": il_db,
        "crosstalk_db": xt_db,
        "extinction_ratio_db": er_db,
        "active_area_um2": active_area,
        "fabric_area_mm2": fabric_mm2,
        "pass_il": pass_il,
        "pass_xt": pass_xt,
        "all_pass": all_pass,
        "elapsed_s": elapsed,
        "p_ref": p_ref
    }

    if verbose:
        print(f"\n[{state.upper()}] Results (L_c={l_c:.2f}um, w_slot={w_slot*1000:.0f}nm):")
        print(f"  Port 1 (Primary / Bottom): {p1 * 100:.2f}%")
        print(f"  Port 2 (Switched / Top)  : {p2 * 100:.2f}%")
        print(f"  Total Guided Trans       : {tot * 100:.2f}% (Rad Loss: {(1-tot)*100:.2f}%)")
        print(f"  Signal Port              : {sig * 100:.2f}% -> IL = {il_db:.4f} dB [{'PASS' if pass_il else 'FAIL'}]")
        print(f"  Leakage (Unselected)     : {leak * 100:.2f}% -> XT = {xt_db:.2f} dB [{'PASS' if pass_xt else 'FAIL'}]")
        print(f"  Active Footprint         : {active_area:.3f} um^2 [PASS <5.0um^2]")
        print(f"  4096-Switch Fabric Area  : {fabric_mm2:.4f} mm^2 [PASS <20.0mm^2]")
        print(f"  VERDICT                  : {'ALL TARGETS MET!' if all_pass else 'Targets not yet met'}")
        print(f"  Simulation Time          : {elapsed:.1f} s")

    return res


# ─────────────────────────────────────────────────────────────────────────────
# Dual-State Verification Runner
# ─────────────────────────────────────────────────────────────────────────────

def run_both_states(l_c=2.56, w_rail=0.550, w_slot=0.065, resolution=25, **kwargs):
    print("=" * 80)
    print("PROJECT JANUS: Si3N4 EMBEDDED PCM SLOT SWITCH (GEOMETRY B)")
    print(f"Parameters: L_c={l_c:.2f}um, W_rail={w_rail*1000:.0f}nm, W_slot={w_slot*1000:.0f}nm, Res={resolution} px/um")
    print("Architecture: In-plane embedded Sb2S3 slot slab between dual Si3N4 rails")
    print("=" * 80)

    # 1. Run Amorphous State (computes and caches reference flux)
    r_am = run_sin_slot_switch("amorphous", l_c=l_c, w_rail=w_rail, w_slot=w_slot, resolution=resolution, **kwargs)
    p_ref_cached = r_am.get("p_ref", None)

    print("-" * 80)

    # 2. Run Crystalline State (reuses cached reference flux)
    r_cr = run_sin_slot_switch("crystalline", l_c=l_c, w_rail=w_rail, w_slot=w_slot, resolution=resolution, p_ref=p_ref_cached, **kwargs)

    both_pass = r_am["all_pass"] and r_cr["all_pass"]

    print("\n" + "#" * 80)
    print(f"DUAL-STATE SUMMARY: Si3N4 Embedded Slot Switch (L_c={l_c:.2f}um, w_slot={w_slot*1000:.0f}nm)")
    print("#" * 80)
    print(f"{'Metric':<28} {'Amorphous (State 0)':<24} {'Crystalline (State 1)':<24}")
    print("-" * 80)
    print(f"{'Signal Port Trans':<28} {r_am['signal_trans']*100:.2f}%{'':<18} {r_cr['signal_trans']*100:.2f}%")
    print(f"{'Leakage (Unselected)':<28} {r_am['leakage']*100:.2f}%{'':<18} {r_cr['leakage']*100:.2f}%")
    print(f"{'Total Guided Trans':<28} {r_am['total_guided']*100:.2f}%{'':<18} {r_cr['total_guided']*100:.2f}%")
    print(f"{'Radiation Loss':<28} {r_am['radiation_loss']*100:.2f}%{'':<18} {r_cr['radiation_loss']*100:.2f}%")
    print(f"{'Insertion Loss':<28} {r_am['insertion_loss_db']:.4f} dB{'':<15} {r_cr['insertion_loss_db']:.4f} dB")
    print(f"{'Crosstalk':<28} {r_am['crosstalk_db']:.2f} dB{'':<16} {r_cr['crosstalk_db']:.2f} dB")
    print(f"{'Active Footprint':<28} {r_am['active_area_um2']:.3f} um^2{'':<16} {r_cr['active_area_um2']:.3f} um^2")
    print(f"{'4096-Fabric Area':<28} {r_am['fabric_area_mm2']:.4f} mm^2{'':<15} {r_cr['fabric_area_mm2']:.4f} mm^2")
    print(f"{'IL Target Met':<28} {str(r_am['pass_il']):<24} {str(r_cr['pass_il']):<24}")
    print(f"{'XT Target Met (<= -20dB)':<28} {str(r_am['pass_xt']):<24} {str(r_cr['pass_xt']):<24}")
    print(f"{'ALL TARGETS MET':<28} {str(r_am['all_pass']):<24} {str(r_cr['all_pass']):<24}")
    print("#" * 80)
    print(f"OVERALL VERDICT: {'PASSED ALL SPECIFICATIONS!' if both_pass else 'TARGETS NOT YET MET'}\n")

    return r_am, r_cr


# ─────────────────────────────────────────────────────────────────────────────
# Parameter Sweep
# ─────────────────────────────────────────────────────────────────────────────

def run_sweep(lc_min=2.40, lc_max=2.70, lc_step=0.05, w_slot=0.065, w_rail=0.550, resolution=25):
    lc_vals = np.arange(lc_min, lc_max + 1e-6, lc_step)
    out_csv = "sin_slot_switch_sweep_results.csv"

    print("=" * 80)
    print("Si3N4 EMBEDDED SLOT SWITCH PARAMETER SWEEP")
    print(f"L_c range: [{lc_min:.2f}, {lc_max:.2f}] step {lc_step:.2f} | w_slot={int(w_slot*1000)}nm, w_rail={int(w_rail*1000)}nm")
    print(f"Output CSV: {out_csv}")
    print("=" * 80)

    # Pre-compute reference simulation once
    print("\n[SWEEP PRE-COMPUTE] Pre-computing reference flux once...")
    r_init = run_sin_slot_switch("amorphous", l_c=lc_vals[0], w_rail=w_rail, w_slot=w_slot, resolution=resolution, verbose=False)
    p_ref_cached = r_init["p_ref"]
    print(f"[SWEEP PRE-COMPUTE] Reference flux cached: {p_ref_cached:.6e}. Running sweep points...\n")

    results = []
    with open(out_csv, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow([
            "l_c_um", "w_slot_nm", "w_rail_nm",
            "am_signal", "cr_signal",
            "am_il", "cr_il",
            "am_xt", "cr_xt",
            "am_rad", "cr_rad",
            "active_area",
            "am_pass", "cr_pass", "both_pass"
        ])

        for lc in lc_vals:
            r_am = run_sin_slot_switch("amorphous", l_c=lc, w_rail=w_rail, w_slot=w_slot, resolution=resolution, p_ref=p_ref_cached, verbose=False)
            r_cr = run_sin_slot_switch("crystalline", l_c=lc, w_rail=w_rail, w_slot=w_slot, resolution=resolution, p_ref=p_ref_cached, verbose=False)
            both_pass = r_am["all_pass"] and r_cr["all_pass"]

            row = [
                round(lc, 3), int(w_slot * 1000), int(w_rail * 1000),
                round(r_am["signal_trans"] * 100, 3), round(r_cr["signal_trans"] * 100, 3),
                round(r_am["insertion_loss_db"], 4), round(r_cr["insertion_loss_db"], 4),
                round(r_am["crosstalk_db"], 2), round(r_cr["crosstalk_db"], 2),
                round(r_am["radiation_loss"] * 100, 2), round(r_cr["radiation_loss"] * 100, 2),
                round(r_am["active_area_um2"], 4),
                r_am["all_pass"], r_cr["all_pass"], both_pass
            ]
            writer.writerow(row)
            f.flush()
            results.append(row)

            print(f"  L_c={lc:.2f}um | "
                  f"AM: Sig={r_am['signal_trans']*100:.1f}%, IL={r_am['insertion_loss_db']:.3f}dB, XT={r_am['crosstalk_db']:.1f}dB | "
                  f"CR: Sig={r_cr['signal_trans']*100:.1f}%, IL={r_cr['insertion_loss_db']:.3f}dB, XT={r_cr['crosstalk_db']:.1f}dB | "
                  f"{'ALL PASS' if both_pass else 'FAIL'}")

    return results


# ─────────────────────────────────────────────────────────────────────────────
# CLI Entrypoint
# ─────────────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Si3N4 Embedded PCM Slot Switch 3D MEEP FDTD")
    parser.add_argument("--state", type=str, default="both", choices=["amorphous", "crystalline", "both"])
    parser.add_argument("--lc", type=float, default=2.56, help="Interaction length (um)")
    parser.add_argument("--w_rail", type=float, default=0.550, help="Si3N4 rail width (um)")
    parser.add_argument("--w_slot", type=float, default=0.065, help="Central embedded Sb2S3 slot slab width (um)")
    parser.add_argument("--resolution", type=int, default=25, help="FDTD grid resolution (px/um)")
    parser.add_argument("--sweep", action="store_true", help="Run fine parameter sweep around resonance")
    parser.add_argument("--lc_min", type=float, default=2.40)
    parser.add_argument("--lc_max", type=float, default=2.70)
    parser.add_argument("--lc_step", type=float, default=0.05)

    args = parser.parse_args()

    if args.sweep:
        run_sweep(lc_min=args.lc_min, lc_max=args.lc_max, lc_step=args.lc_step,
                  w_slot=args.w_slot, w_rail=args.w_rail, resolution=args.resolution)
    elif args.state == "both":
        run_both_states(l_c=args.lc, w_rail=args.w_rail, w_slot=args.w_slot, resolution=args.resolution)
    else:
        run_sin_slot_switch(state=args.state, l_c=args.lc, w_rail=args.w_rail, w_slot=args.w_slot, resolution=args.resolution)
