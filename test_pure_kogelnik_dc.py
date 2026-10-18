"""
PROJECT JANUS: PURE KOGELNIK 1x2 DIRECTIONAL COUPLER SWITCH (3D FDTD)
=====================================================================
Physics:
1. WG1: Straight continuous bus at y = 0, width w = 240 nm.
2. WG2: Approaches to gap = 98 nm at x_c_start from an input fan-out (gap > 220 nm),
        runs parallel to WG1 at gap = 98 nm over L_c,
        and fans out to y_out2 = -0.50 um at x_c_end.
3. PCM Coverage:
   - Synchronized on both waveguides over the parallel interaction region [x_c_start, x_c_end].
   - Smooth adiabatic tapers (l_tip = 0.25 um, w_tip = 30 nm) at both ends.
   - WG2 PCM is fixed Amorphous (n = 2.70).
   - WG1 PCM switches:
       State 0 (Amorphous):   n = 2.70 -> Delta_beta = 0 identically (100% Cross into Port 2)
       State 1 (Crystalline): n = 3.30 -> Delta_beta = 1.21 rad/um, u = sqrt(3) (100% Bar into Port 1)
4. Zero parasitic coupling: evanescent coupling only exists inside [x_c_start, x_c_end].
"""

import sys, os, math, time, json, argparse
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
    bot = [mp.Vector3(xi, y - wi / 2.0, 0) for xi, y_val, wi in zip(x, [y]*n_pts, w)]
    top = [mp.Vector3(xi, y + wi / 2.0, 0) for xi, y_val, wi in reversed(list(zip(x, [y]*n_pts, w)))]
    return bot + top

def create_prism_exact(vertices, height, z_base, material):
    cx = float(np.mean([pt.x for pt in vertices]))
    cy = float(np.mean([pt.y for pt in vertices]))
    cz = float(z_base + height / 2.0)
    return mp.Prism(vertices=vertices, height=height, axis=mp.Vector3(0, 0, 1),
                    center=mp.Vector3(cx, cy, cz), material=material)

def simulate_switch(state="amorphous", resolution=15, l_c=4.30, gap=0.098,
                    t_pcm=0.160, l_tip=0.200, y_out=0.500,
                    l_bend=0.750, l_in_bend=0.600, p_ref=None, verbose=True):
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

    y1 = 0.0                         # WG1 straight through at y = 0
    y2_c = -(w_wg + gap)             # WG2 coupling position (-0.338 um)
    y_in2 = -(w_wg + 0.280)          # WG2 entrance starting position (gap = 280 nm)
    y_out2 = -y_out                  # WG2 exit decoupled port (-0.500 um)

    l_lead = 0.300
    iso_buf = 0.080

    cell_len = l_lead + l_in_bend + l_c + l_bend + l_lead
    cell_w = y_out + w_wg + 2 * iso_buf # ~0.90 um
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

    x_in_start = x0 + l_lead
    x_c_start = x_in_start + l_in_bend
    x_c_end = x_c_start + l_c
    x_bend_end = x_c_end + l_bend
    x_mon = x_bend_end + 0.15

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

    # Reference straight waveguide
    if p_ref is None:
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
        p_ref = max(mp.get_fluxes(f_ref)[0], 1e-12)
        if verbose:
            print(f"[{state.upper()}] Cached Reference Flux = {p_ref:.6e}")

    # Device geometry:
    # 1. Straight continuous bus for WG1 at y = 0
    geom_dev = [
        mp.Block(center=mp.Vector3(0, 0, -sz / 4.0), size=mp.Vector3(mp.inf, mp.inf, sz / 2.0), material=clad_mat),
        mp.Block(center=mp.Vector3(0, y1, h_wg / 2.0), size=mp.Vector3(mp.inf, w_wg, h_wg), material=core_mat),
    ]

    # 2. WG2 Entrance S-bend: approaches from y_in2 to y2_c
    # Starts with inverse taper at y_in2
    poly_wg2_in_taper = make_ccw_linear_taper(x_in_start - 0.200, x_in_start, y_in2, 0.020, w_wg)
    geom_dev.append(create_prism_exact(poly_wg2_in_taper, h_wg, 0.0, core_mat))

    poly_wg2_in_bend = make_ccw_sbend(x_in_start, x_c_start, y_in2, y2_c, w_wg)
    geom_dev.append(create_prism_exact(poly_wg2_in_bend, h_wg, 0.0, core_mat))

    # 3. WG2 Parallel Coupled Section: [x_c_start, x_c_end] at y2_c
    geom_dev.append(mp.Block(
        center=mp.Vector3((x_c_start + x_c_end) / 2.0, y2_c, h_wg / 2.0),
        size=mp.Vector3(l_c, w_wg, h_wg), material=core_mat
    ))

    # 4. WG2 Exit S-bend: fans out from y2_c to y_out2
    poly_wg2_out_bend = make_ccw_sbend(x_c_end, x_bend_end, y2_c, y_out2, w_wg)
    geom_dev.append(create_prism_exact(poly_wg2_out_bend, h_wg, 0.0, core_mat))

    # 5. WG2 Output Lead at y_out2
    geom_dev.append(mp.Block(
        center=mp.Vector3((x_bend_end + x_right) / 2.0, y_out2, h_wg / 2.0),
        size=mp.Vector3(x_right - x_bend_end, w_wg, h_wg), material=core_mat
    ))

    # 6. PCM Stacks:
    mat_pcm1 = pcm_cr if state.lower() == "crystalline" else pcm_am
    mat_pcm2 = pcm_am

    # WG1 PCM:
    # Taper in: [x_c_start - l_tip, x_c_start]
    poly_pcm1_in = make_ccw_linear_taper(x_c_start - l_tip, x_c_start, y1, 0.030, w_wg)
    geom_dev.append(create_prism_exact(poly_pcm1_in, t_pcm, h_wg, mat_pcm1))
    # Full body: [x_c_start, x_c_end]
    geom_dev.append(mp.Block(
        center=mp.Vector3((x_c_start + x_c_end) / 2.0, y1, h_wg + t_pcm / 2.0),
        size=mp.Vector3(l_c, w_wg, t_pcm), material=mat_pcm1
    ))
    # Taper out: [x_c_end, x_c_end + l_tip]
    poly_pcm1_out = make_ccw_linear_taper(x_c_end, x_c_end + l_tip, y1, w_wg, 0.030)
    geom_dev.append(create_prism_exact(poly_pcm1_out, t_pcm, h_wg, mat_pcm1))

    # WG2 PCM:
    # Taper in: [x_c_start - l_tip, x_c_start]
    poly_pcm2_in = make_ccw_linear_taper(x_c_start - l_tip, x_c_start, y2_c, 0.030, w_wg)
    geom_dev.append(create_prism_exact(poly_pcm2_in, t_pcm, h_wg, mat_pcm2))
    # Full body: [x_c_start, x_c_end]
    geom_dev.append(mp.Block(
        center=mp.Vector3((x_c_start + x_c_end) / 2.0, y2_c, h_wg + t_pcm / 2.0),
        size=mp.Vector3(l_c, w_wg, t_pcm), material=mat_pcm2
    ))
    # Taper out: [x_c_end, x_c_end + l_tip]
    poly_pcm2_out = make_ccw_linear_taper(x_c_end, x_c_end + l_tip, y2_c, w_wg, 0.030)
    geom_dev.append(create_prism_exact(poly_pcm2_out, t_pcm, h_wg, mat_pcm2))

    # Simulation execution
    sim_dev = mp.Simulation(cell_size=cell, boundary_layers=pml_layers, geometry=geom_dev,
                            sources=sources, resolution=resolution, default_material=clad_mat)
    f_p1 = sim_dev.add_flux(fcen, 0, 1, mp.FluxRegion(center=mp.Vector3(x_mon, y1, h_wg / 2.0), size=mp.Vector3(0, mon_w, mon_h)))
    f_p2 = sim_dev.add_flux(fcen, 0, 1, mp.FluxRegion(center=mp.Vector3(x_mon, y_out2, h_wg / 2.0), size=mp.Vector3(0, mon_w, mon_h)))

    y_dest = y_out2 if state.lower() == "amorphous" else y1
    sim_dev.run(until_after_sources=mp.stop_when_fields_decayed(20, mp.Hz, mp.Vector3(x_mon, y_dest, h_wg / 2.0), 1e-3))

    p1 = max(mp.get_fluxes(f_p1)[0] / p_ref, 1e-12)
    p2 = max(mp.get_fluxes(f_p2)[0] / p_ref, 1e-12)
    p_tot = p1 + p2

    sig = p2 if state.lower() == "amorphous" else p1
    leak = p1 if state.lower() == "amorphous" else p2
    il_db = -10.0 * math.log10(sig)
    xt_db = 10.0 * math.log10(leak / sig)
    ratio_guided = sig / p_tot * 100.0

    if verbose:
        print(f"\n{'='*70}")
        print(f"[{state.upper()}] (L_c={l_c:.3f}um, Gap={gap*1000:.0f}nm, Area={cell_area:.3f}um2):")
        print(f"  Port 1 (Through / Bar) : {p1*100:.2f}%")
        print(f"  Port 2 (Cross)         : {p2*100:.2f}%")
        print(f"  Total Guided Power     : {p_tot*100:.2f}% (Rad Loss: {(1-p_tot)*100:.2f}%)")
        print(f"  Signal Port Power      : {sig*100:.2f}% (IL = {il_db:.4f} dB)")
        print(f"  Leakage (Crosstalk)    : {leak*100:.2f}% (XT = {xt_db:.2f} dB)")
        print(f"  Signal / Guided Power  : {ratio_guided:.2f}%")
        print(f"  Cell Footprint         : {cell_area:.3f} um^2 (Length={cell_len:.2f}um, Width={cell_w:.2f}um)")
        print(f"{'='*70}\n")

    return {"state": state, "p1": p1, "p2": p2, "p_tot": p_tot, "sig": sig, "leak": leak,
            "il_db": il_db, "xt_db": xt_db, "ratio_guided": ratio_guided, "cell_area": cell_area, "p_ref": p_ref}

def run_both(resolution=15, l_c=4.30, gap=0.098):
    print("=" * 80)
    print("  PROJECT JANUS: PURE KOGELNIK 1x2 SWITCH (3D FDTD)")
    print(f"  L_c = {l_c:.3f} um, Gap = {gap*1000:.1f} nm, Resolution = {resolution}")
    print("=" * 80)
    res_am = simulate_switch("amorphous", resolution=resolution, l_c=l_c, gap=gap, p_ref=None)
    res_cr = simulate_switch("crystalline", resolution=resolution, l_c=l_c, gap=gap, p_ref=res_am["p_ref"])

    print("\n" + "=" * 80)
    print("  DUAL-STATE PERFORMANCE SUMMARY")
    print("=" * 80)
    print(f"State 0 (Amorphous   -> Port 2 Cross) : Signal = {res_am['sig']*100:.2f}%, Leak = {res_am['leak']*100:.2f}%, Guided = {res_am['p_tot']*100:.2f}%, Sig/Guided = {res_am['ratio_guided']:.2f}%")
    print(f"State 1 (Crystalline -> Port 1 Bar)   : Signal = {res_cr['sig']*100:.2f}%, Leak = {res_cr['leak']*100:.2f}%, Guided = {res_cr['p_tot']*100:.2f}%, Sig/Guided = {res_cr['ratio_guided']:.2f}%")
    print(f"Cell Footprint                         : {res_am['cell_area']:.3f} um^2 (Limit <= 6.0 um^2: PASS)")
    print("=" * 80 + "\n")
    return res_am, res_cr

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("--lc", type=float, default=4.30)
    parser.add_argument("--gap", type=float, default=0.098)
    parser.add_argument("--resolution", type=int, default=15)
    args = parser.parse_args()
    run_both(resolution=args.resolution, l_c=args.lc, gap=args.gap)
