"""
PROJECT JANUS: PRODUCTION SUB-6 UM^2 1x2 OPTICAL SWITCH CELL (3D FDTD)
======================================================================
Operating Wavelength: lambda_0 = 1064 nm
Waveguide Core      : Silicon (n = 3.480), w = 240 nm, h = 220 nm
Active PCM Stack    : Sb2S3 (n_am = 2.70, n_cr = 3.30), t_pcm = 160 nm
Cladding            : Silicon Dioxide (n = 1.449)

Physical Switching Mechanism:
  - Exact Kogelnik Phase-Mismatch Switching Invariant:
      u = Delta_beta / (2 * kappa) = sqrt(3) ~ 1.73205
  - Symmetric Amorphous Stack: Matching Sb2S3 patch on WG2 enforces Delta_beta = 0 in State 0.
  - State 0 (Amorphous)  : Resonant Cross-Coupling -> Port 2 (Cross)
  - State 1 (Crystalline): Exact Zero-Beat Quenching -> Port 1 (Through / Bar)
  - Footprint            : <= 6.0 um^2 Total Cell Footprint
"""

import sys, os, math, time, json, argparse
import numpy as np
import meep as mp

def make_ccw_sbend(x_start, x_end, y_start, y_end, w1, w2=None, n_pts=30):
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

def make_ccw_pcm_patch(x_center, y_center, l_body, l_tip, w, w_tip=0.030, n_pts=15):
    x_l_tip = x_center - l_body / 2.0 - l_tip
    x_l_body = x_center - l_body / 2.0
    x_r_body = x_center + l_body / 2.0
    x_r_tip = x_center + l_body / 2.0 + l_tip

    u_l = np.linspace(0, 1, n_pts)[:-1]
    xl = x_l_tip + u_l * l_tip
    wl = w_tip + (w - w_tip) * (2 * u_l - u_l**2)

    u_r = np.linspace(0, 1, n_pts)[1:]
    xr = x_r_body + u_r * l_tip
    wr = w_tip + (w - w_tip) * (1.0 - u_r**2)

    x_all = list(xl) + [x_l_body, x_r_body] + list(xr)
    w_all = list(wl) + [w, w] + list(wr)

    bot = [mp.Vector3(xi, y_center - wi / 2.0, 0) for xi, wi in zip(x_all, w_all)]
    top = [mp.Vector3(xi, y_center + wi / 2.0, 0) for xi, wi in reversed(list(zip(x_all, w_all)))]
    return bot + top

def create_prism_exact(vertices, height, z_base, material):
    cx = float(np.mean([pt.x for pt in vertices]))
    cy = float(np.mean([pt.y for pt in vertices]))
    cz = float(z_base + height / 2.0)
    return mp.Prism(vertices=vertices, height=height, axis=mp.Vector3(0, 0, 1),
                    center=mp.Vector3(cx, cy, cz), material=material)

def simulate_switch_state(state="amorphous", resolution=20, l_c=4.30, gap=0.098,
                          t_pcm=0.160, y_pitch=0.450, l_bend=1.150, w_bend=0.280,
                          p_ref=None, verbose=True):
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

    cell_len = l_in + l_c + l_bend + l_out
    cell_w = y_pitch + w_bend/2.0 + w_wg/2.0 + 2 * iso_buf
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
    x_c_end = x_c_start + l_c
    x_bend_end = x_c_end + l_bend

    x_mon = x_bend_end + 0.15
    mon_w = max(w_wg, w_bend) * 2.0
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
        geom_ref = [
            mp.Block(center=mp.Vector3(0, 0, -sz / 4.0), size=mp.Vector3(mp.inf, mp.inf, sz / 2.0), material=clad_mat),
            mp.Block(center=mp.Vector3(0, y1, h_wg / 2.0), size=mp.Vector3(mp.inf, w_wg, h_wg), material=core_mat)
        ]
        sim_ref = mp.Simulation(cell_size=cell, boundary_layers=pml_layers, geometry=geom_ref,
                                sources=sources, resolution=resolution, default_material=clad_mat)
        f_ref = sim_ref.add_flux(fcen, 0, 1, mp.FluxRegion(center=mp.Vector3(x_mon, y1, h_wg / 2.0), size=mp.Vector3(0, mon_w, mon_h)))
        sim_ref.run(until_after_sources=mp.stop_when_fields_decayed(20, mp.Hz, mp.Vector3(x_mon, y1, h_wg / 2.0), 1e-3))
        p_ref = max(mp.get_fluxes(f_ref)[0], 1e-12)

    # Device geometry:
    geom_dev = [
        mp.Block(center=mp.Vector3(0, 0, -sz / 4.0), size=mp.Vector3(mp.inf, mp.inf, sz / 2.0), material=clad_mat),
        mp.Block(center=mp.Vector3(0, y1, h_wg / 2.0), size=mp.Vector3(mp.inf, w_wg, h_wg), material=core_mat),
    ]

    # WG2 entrance adiabatic inverse taper
    taper_len = min(0.350, l_in)
    poly_wg2_in = make_ccw_linear_taper(x_c_start - taper_len, x_c_start, y2_c, 0.020, w_wg)
    geom_dev.append(create_prism_exact(poly_wg2_in, h_wg, 0.0, core_mat))

    # WG2 straight coupled section
    geom_dev.append(mp.Block(
        center=mp.Vector3((x_c_start + x_c_end) / 2.0, y2_c, h_wg / 2.0),
        size=mp.Vector3(l_c, w_wg, h_wg),
        material=core_mat
    ))

    # WG2 output S-bend fanout with smooth widening to w_bend
    poly_wg2_bend = make_ccw_sbend(x_c_end, x_bend_end, y2_c, y_out2, w_wg, w_bend)
    geom_dev.append(create_prism_exact(poly_wg2_bend, h_wg, 0.0, core_mat))

    # WG2 output lead at y_out2
    geom_dev.append(mp.Block(
        center=mp.Vector3((x_bend_end + x_right) / 2.0, y_out2, h_wg / 2.0),
        size=mp.Vector3(x_right - x_bend_end, w_bend, h_wg),
        material=core_mat
    ))

    # PCM patches strictly inside coupling region
    mat_pcm1 = pcm_cr if state.lower() == "crystalline" else pcm_am
    l_tip = 0.300
    l_body = max(l_c - 2.0 * l_tip, 0.100)
    x_pcm_center = (x_c_start + x_c_end) / 2.0
    poly_pcm1 = make_ccw_pcm_patch(x_pcm_center, y1, l_body, l_tip, w_wg)
    geom_dev.append(create_prism_exact(poly_pcm1, t_pcm, h_wg, mat_pcm1))

    # Matching Amorphous PCM on WG2 for exact phase matching in State 0
    poly_pcm2 = make_ccw_pcm_patch(x_pcm_center, y2_c, l_body, l_tip, w_wg)
    geom_dev.append(create_prism_exact(poly_pcm2, t_pcm, h_wg, pcm_am))

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

    if verbose:
        print(f"\n{'='*70}")
        print(f"[{state.upper()}] (L_c={l_c:.3f}um, Gap={gap*1000:.0f}nm, t_pcm={t_pcm*1000:.0f}nm, Area={cell_area:.3f}um2):")
        print(f"  Port 1 (Through / Bar) : {p1*100:.2f}%")
        print(f"  Port 2 (Cross)         : {p2*100:.2f}%")
        print(f"  Total Guided Power     : {p_tot*100:.2f}% (Rad Loss: {(1-p_tot)*100:.2f}%)")
        print(f"  Signal Port Power      : {sig*100:.2f}% (IL = {il_db:.4f} dB)")
        print(f"  Leakage (Crosstalk)    : {leak*100:.2f}% (XT = {xt_db:.2f} dB)")
        print(f"  Cell Footprint         : {cell_area:.3f} um^2 (Length={cell_len:.2f}um, Width={cell_w:.2f}um)")
        print(f"{'='*70}\n")

    return {
        "state": state, "p1": p1, "p2": p2, "p_tot": p_tot, "sig": sig, "leak": leak,
        "il_db": il_db, "xt_db": xt_db, "cell_area": cell_area, "cell_len": cell_len,
        "cell_w": cell_w, "p_ref": p_ref
    }

def main():
    parser = argparse.ArgumentParser(description="3D Vector FDTD Simulation of Sub-6 um^2 Optical Switch")
    parser.add_argument("--state", type=str, default="both", choices=["amorphous", "crystalline", "both"])
    parser.add_argument("--resolution", type=int, default=20)
    parser.add_argument("--lc", type=float, default=4.30)
    parser.add_argument("--gap", type=float, default=0.098)
    parser.add_argument("--tpcm", type=float, default=0.160)
    parser.add_argument("--ypitch", type=float, default=0.450)
    parser.add_argument("--save", type=str, default="verified_sub6um_switch_3d.json")
    args = parser.parse_args()

    print("=" * 80)
    print("  PROJECT JANUS: 3D VECTOR MAXWELL FDTD PRODUCTION VERIFIER")
    print(f"  Resolution: {args.resolution} px/um | L_c: {args.lc:.3f} um | Gap: {args.gap*1000:.0f} nm | t_pcm: {args.tpcm*1000:.0f} nm")
    print("=" * 80)

    if args.state == "both":
        res_am = simulate_switch_state("amorphous", resolution=args.resolution, l_c=args.lc,
                                       gap=args.gap, t_pcm=args.tpcm, y_pitch=args.ypitch)
        res_cr = simulate_switch_state("crystalline", resolution=args.resolution, l_c=args.lc,
                                       gap=args.gap, t_pcm=args.tpcm, y_pitch=args.ypitch,
                                       p_ref=res_am["p_ref"])

        metrics = {
            "device": "Janus Sub-6 um^2 1x2 Non-Volatile Optical Switch Cell",
            "wavelength_nm": 1064,
            "architecture": "Phase-Mismatched Directional Coupler with Symmetric Amorphous Stack",
            "parameters": {
                "l_c_um": args.lc,
                "gap_nm": round(args.gap * 1000, 1),
                "t_pcm_nm": round(args.tpcm * 1000, 1),
                "w_wg_nm": 240,
                "h_wg_nm": 220,
                "w_bend_nm": 280,
                "y_pitch_nm": round(args.ypitch * 1000, 1),
                "resolution": args.resolution
            },
            "metrics": {
                "state_0_amorphous": {
                    "signal_port": "Port 2 (Cross)",
                    "transmission_pct": round(res_am["sig"] * 100, 3),
                    "insertion_loss_db": round(res_am["il_db"], 4),
                    "crosstalk_leakage_pct": round(res_am["leak"] * 100, 3),
                    "crosstalk_db": round(res_am["xt_db"], 2),
                    "radiation_loss_pct": round((1 - res_am["p_tot"]) * 100, 3)
                },
                "state_1_crystalline": {
                    "signal_port": "Port 1 (Through / Bar)",
                    "transmission_pct": round(res_cr["sig"] * 100, 3),
                    "insertion_loss_db": round(res_cr["il_db"], 4),
                    "crosstalk_leakage_pct": round(res_cr["leak"] * 100, 3),
                    "crosstalk_db": round(res_cr["xt_db"], 2),
                    "radiation_loss_pct": round((1 - res_cr["p_tot"]) * 100, 3)
                },
                "footprint": {
                    "length_um": round(res_am["cell_len"], 3),
                    "width_um": round(res_am["cell_w"], 3),
                    "cell_area_um2": round(res_am["cell_area"], 3),
                    "area_ceiling_um2": 6.0,
                    "area_compliant": bool(res_am["cell_area"] <= 6.0)
                }
            }
        }

        with open(args.save, "w") as f:
            json.dump(metrics, f, indent=2)

        print("\n" + "=" * 80)
        print("  FINAL 3D FDTD COMPLIANCE VERIFICATION SUMMARY")
        print("=" * 80)
        print(f"State 0 (Amorphous   -> Port 2 Cross) : Signal = {res_am['sig']*100:.2f}% | IL = {res_am['il_db']:.4f} dB | XT = {res_am['xt_db']:.2f} dB")
        print(f"State 1 (Crystalline -> Port 1 Bar)   : Signal = {res_cr['sig']*100:.2f}% | IL = {res_cr['il_db']:.4f} dB | XT = {res_cr['xt_db']:.2f} dB")
        print(f"Cell Area                             : {res_am['cell_area']:.3f} um^2 (Limit <= 6.0 um^2: {'PASS' if res_am['cell_area'] <= 6.0 else 'FAIL'})")
        print(f"Verified JSON Export                  : {args.save}")
        print("=" * 80)

    else:
        res = simulate_switch_state(args.state, resolution=args.resolution, l_c=args.lc,
                                    gap=args.gap, t_pcm=args.tpcm, y_pitch=args.ypitch)

if __name__ == '__main__':
    main()
