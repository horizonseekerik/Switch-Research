"""
PROJECT JANUS: 3D FULL-WAVE CALIBRATED SWITCH
Exact 3D Verification at L_c ~ 2.80 um, gap = 100 nm
"""
import sys, os, math
import numpy as np
import meep as mp

def run_calibrated_dc(state="amorphous", l_c=2.80, gap=0.100, resolution=20):
    lambda_0 = 1.064
    fcen = 1.0 / lambda_0
    df = 0.1 * fcen

    w_wg = 0.240
    h_wg = 0.220
    t_pcm = 0.150
    l_in = 0.400
    l_tip = 0.350        # Parabolic PCM taper length
    w_tip = 0.020        # Tip width
    l_bend = 1.100       # S-bend fanout length
    y_pitch_out = 0.560  # 560 nm output pitch
    l_out = 0.400
    iso_buf = 0.080

    y_sep = w_wg + gap
    y_wg1 = y_sep / 2.0
    y_wg2 = -y_sep / 2.0
    y_out1 = y_pitch_out / 2.0
    y_out2 = -y_pitch_out / 2.0

    l_straight = l_c + 2 * l_tip
    cell_len = l_in + l_straight + l_bend + l_out
    cell_w = y_pitch_out + 2 * iso_buf
    cell_area = cell_len * cell_w

    dpml = 0.50
    pad_x = 0.35
    pad_y = 0.35
    pad_z = 0.35
    sx = cell_len + 2 * (dpml + pad_x)
    sy = cell_w + 2 * (dpml + pad_y)
    sz = h_wg + t_pcm + 2 * (dpml + pad_z)

    cell = mp.Vector3(sx, sy, sz)
    pml_layers = [mp.PML(dpml)]

    x_left = -sx / 2.0
    x_right = sx / 2.0
    x_min_cell = -cell_len / 2.0

    x_in_end = x_min_cell + l_in
    x_straight_start = x_in_end
    x_straight_end = x_straight_start + l_straight
    x_bend_end = x_straight_end + l_bend

    clad_mat = mp.Medium(index=1.449)
    core_mat = mp.Medium(index=3.480)

    pcm_am = mp.Medium(index=2.70, D_conductivity=2*math.pi*fcen*1.0e-5)
    pcm_cr = mp.Medium(index=3.30, D_conductivity=2*math.pi*fcen*8.0e-5)

    mat_wg1_pcm = pcm_cr if state.lower() == "crystalline" else pcm_am
    mat_wg2_pcm = pcm_am

    # Non-overlapping monitors: 320x320 nm
    mon_w = 0.320
    mon_h = 0.320
    x_mon = x_bend_end + 0.15

    sources = [
        mp.EigenModeSource(
            src=mp.GaussianSource(fcen, fwidth=df),
            center=mp.Vector3(x_left + dpml + 0.10, y_wg1, h_wg / 2.0),
            size=mp.Vector3(0, w_wg * 3.0, h_wg * 3.0),
            eig_band=1,
            direction=mp.X
        )
    ]

    print(f"[{state.upper()}] 1. Reference Straight Waveguide Simulation...")
    geom_ref = [
        mp.Block(center=mp.Vector3(0, 0, -sz/4.0), size=mp.Vector3(mp.inf, mp.inf, sz/2.0), material=clad_mat),
        mp.Block(center=mp.Vector3(0, y_wg1, h_wg/2.0), size=mp.Vector3(mp.inf, w_wg, h_wg), material=core_mat)
    ]
    sim_ref = mp.Simulation(cell_size=cell, boundary_layers=pml_layers, geometry=geom_ref,
                            sources=sources, resolution=resolution, default_material=clad_mat)
    f_ref = sim_ref.add_flux(fcen, 0, 1, mp.FluxRegion(center=mp.Vector3(x_mon, y_wg1, h_wg/2.0), size=mp.Vector3(0, mon_w, mon_h)))
    sim_ref.run(until_after_sources=mp.stop_when_fields_decayed(20, mp.Hz, mp.Vector3(x_mon, y_wg1, h_wg/2.0), 1e-3))
    p_ref = max(mp.get_fluxes(f_ref)[0], 1e-12)
    print(f"[{state.upper()}] Reference Flux = {p_ref:.6e}")

    # Geometry builders
    def make_cosine_sbend(x_start, x_end, y_start, y_end, w, n_pts=30):
        dx = x_end - x_start
        dy = y_end - y_start
        u = np.linspace(0, 1, n_pts)
        x = x_start + u * dx
        y = y_start + dy * (0.5 * (1.0 - np.cos(np.pi * u)))
        top = [mp.Vector3(xi, yi + w / 2.0, 0) for xi, yi in zip(x, y)]
        bot = [mp.Vector3(xi, yi - w / 2.0, 0) for xi, yi in reversed(list(zip(x, y)))]
        return top + bot

    def make_tapered_pcm(x_start, y_c, l_body, l_t, w_p, w_t=0.020, n_pts=20):
        u_l = np.linspace(0, 1, n_pts)[:-1]
        xl = x_start + u_l * l_t
        wl = w_t + (w_p - w_t) * (2 * u_l - u_l**2)
        x_body_start = x_start + l_t
        x_body_end = x_body_start + l_body
        u_r = np.linspace(0, 1, n_pts)[1:]
        xr = x_body_end + u_r * l_t
        wr = w_t + (w_p - w_t) * (1.0 - u_r**2)
        top = [mp.Vector3(xi, y_c + wi/2.0, 0) for xi, wi in zip(xl, wl)]
        top.append(mp.Vector3(x_body_start, y_c + w_p/2.0, 0))
        top.append(mp.Vector3(x_body_end, y_c + w_p/2.0, 0))
        top.extend([mp.Vector3(xi, y_c + wi/2.0, 0) for xi, wi in zip(xr, wr)])
        bot = [mp.Vector3(xi, y_c - wi/2.0, 0) for xi, wi in reversed(list(zip(xr, wr)))]
        bot.append(mp.Vector3(x_body_end, y_c - w_p/2.0, 0))
        bot.append(mp.Vector3(x_body_start, y_c - w_p/2.0, 0))
        bot.extend([mp.Vector3(xi, y_c - wi/2.0, 0) for xi, wi in reversed(list(zip(xl, wl)))])
        return top + bot

    def make_wg2_lead(x_start, x_end, y_c, w_p, w_t=0.020, n_pts=20):
        u = np.linspace(0, 1, n_pts)
        x = x_start + u * (x_end - x_start)
        w = w_t + (w_p - w_t) * u
        top = [mp.Vector3(xi, y_c + wi/2.0, 0) for xi, wi in zip(x, w)]
        bot = [mp.Vector3(xi, y_c - wi/2.0, 0) for xi, wi in reversed(list(zip(x, w)))]
        return top + bot

    pcm_pts1 = make_tapered_pcm(x_straight_start, y_wg1, l_c, l_tip, w_wg, w_tip)
    cx1 = float(np.mean([pt.x for pt in pcm_pts1]))
    cy1 = float(np.mean([pt.y for pt in pcm_pts1]))

    pcm_pts2 = make_tapered_pcm(x_straight_start, y_wg2, l_c, l_tip, w_wg, w_tip)
    cx2 = float(np.mean([pt.x for pt in pcm_pts2]))
    cy2 = float(np.mean([pt.y for pt in pcm_pts2]))

    l_si_taper = 0.350
    x_wg2_tip_start = x_straight_start - l_si_taper

    geom_dev = [
        # Substrate
        mp.Block(center=mp.Vector3(0, 0, -sz/4.0), size=mp.Vector3(mp.inf, mp.inf, sz/2.0), material=clad_mat),
        # 1. WG1 input lead
        mp.Block(center=mp.Vector3((x_left + x_straight_start)/2.0, y_wg1, h_wg/2.0),
                 size=mp.Vector3(x_straight_start - x_left, w_wg, h_wg), material=core_mat),
        # 2. WG2 adiabatic silicon input taper
        mp.Prism(vertices=make_wg2_lead(x_wg2_tip_start, x_straight_start, y_wg2, w_wg, w_tip),
                 height=h_wg, axis=mp.Vector3(0, 0, 1),
                 center=mp.Vector3((x_wg2_tip_start + x_straight_start)/2.0, y_wg2, h_wg/2.0), material=core_mat),
        # 3. Straight coupling section
        mp.Block(center=mp.Vector3((x_straight_start + x_straight_end)/2.0, y_wg1, h_wg/2.0),
                 size=mp.Vector3(l_straight, w_wg, h_wg), material=core_mat),
        mp.Block(center=mp.Vector3((x_straight_start + x_straight_end)/2.0, y_wg2, h_wg/2.0),
                 size=mp.Vector3(l_straight, w_wg, h_wg), material=core_mat),
        # 4. PCM patches on WG1 and WG2
        mp.Prism(vertices=pcm_pts1, height=t_pcm, axis=mp.Vector3(0, 0, 1),
                 center=mp.Vector3(cx1, cy1, h_wg + t_pcm/2.0), material=mat_wg1_pcm),
        mp.Prism(vertices=pcm_pts2, height=t_pcm, axis=mp.Vector3(0, 0, 1),
                 center=mp.Vector3(cx2, cy2, h_wg + t_pcm/2.0), material=mat_wg2_pcm),
        # 5. Output S-bends
        mp.Prism(vertices=make_cosine_sbend(x_straight_end, x_bend_end, y_wg1, y_out1, w_wg),
                 height=h_wg, axis=mp.Vector3(0, 0, 1), material=core_mat),
        mp.Prism(vertices=make_cosine_sbend(x_straight_end, x_bend_end, y_wg2, y_out2, w_wg),
                 height=h_wg, axis=mp.Vector3(0, 0, 1), material=core_mat),
        # 6. Output leads
        mp.Block(center=mp.Vector3((x_bend_end + x_right)/2.0, y_out1, h_wg/2.0),
                 size=mp.Vector3(x_right - x_bend_end, w_wg, h_wg), material=core_mat),
        mp.Block(center=mp.Vector3((x_bend_end + x_right)/2.0, y_out2, h_wg/2.0),
                 size=mp.Vector3(x_right - x_bend_end, w_wg, h_wg), material=core_mat),
    ]

    print(f"[{state.upper()}] 2. Running Calibrated Device Simulation (gap={gap*1000:.0f}nm, L_c={l_c:.2f}um)...")
    sim_dev = mp.Simulation(cell_size=cell, boundary_layers=pml_layers, geometry=geom_dev,
                            sources=sources, resolution=resolution, default_material=clad_mat)
    f1 = sim_dev.add_flux(fcen, 0, 1, mp.FluxRegion(center=mp.Vector3(x_mon, y_out1, h_wg/2.0), size=mp.Vector3(0, mon_w, mon_h)))
    f2 = sim_dev.add_flux(fcen, 0, 1, mp.FluxRegion(center=mp.Vector3(x_mon, y_out2, h_wg/2.0), size=mp.Vector3(0, mon_w, mon_h)))

    y_dest = y_out2 if state.lower() == "amorphous" else y_out1
    sim_dev.run(until_after_sources=mp.stop_when_fields_decayed(20, mp.Hz, mp.Vector3(x_mon, y_dest, h_wg/2.0), 1e-3))

    p1 = max(mp.get_fluxes(f1)[0] / p_ref, 1e-12)
    p2 = max(mp.get_fluxes(f2)[0] / p_ref, 1e-12)
    tot = p1 + p2

    sig = p2 if state.lower() == "amorphous" else p1
    leak = p1 if state.lower() == "amorphous" else p2
    il_db = -10.0 * math.log10(max(sig, 1e-6))
    xt_db = 10.0 * math.log10(max(leak / max(sig, 1e-6), 1e-6))
    er_db = -xt_db

    print(f"\n======================================================================")
    print(f"CALIBRATED SWITCH [{state.upper()}] (L_c = {l_c:.2f} um, gap = {gap*1000:.0f} nm):")
    print(f"  Footprint Area  : {cell_area:.3f} um^2 ({cell_len:.2f} um x {cell_w:.2f} um)")
    print(f"  Port 1 (Bar)    : {p1*100.0:.2f}%")
    print(f"  Port 2 (Cross)  : {p2*100.0:.2f}%")
    print(f"  Total Trans     : {tot*100.0:.2f}%")
    print(f"  Insertion Loss  : {il_db:.3f} dB")
    print(f"  Crosstalk       : {xt_db:.2f} dB")
    print(f"  Extinction Ratio: {er_db:.2f} dB")
    print("======================================================================\n")
    return {"state": state, "l_c": l_c, "gap": gap, "p1": p1, "p2": p2, "tot": tot, "il": il_db, "xt": xt_db, "area": cell_area}

if __name__ == "__main__":
    st = sys.argv[1] if len(sys.argv) > 1 else "amorphous"
    lc = float(sys.argv[2]) if len(sys.argv) > 2 else 2.80
    gp = float(sys.argv[3]) if len(sys.argv) > 3 else 0.100
    run_calibrated_dc(state=st, l_c=lc, gap=gp)
