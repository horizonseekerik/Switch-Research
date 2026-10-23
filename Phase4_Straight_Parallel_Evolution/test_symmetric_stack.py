import sys, os, math, time, json
import numpy as np
import meep as mp

def make_ccw_sbend(x_start, x_end, y_start, y_end, w, n_pts=30):
    dx = x_end - x_start
    dy = y_end - y_start
    u = np.linspace(0, 1, n_pts)
    x = x_start + u * dx
    y = y_start + dy * (3 * u**2 - 2 * u**3)
    bot = [mp.Vector3(xi, yi - w / 2.0, 0) for xi, yi in zip(x, y)]
    top = [mp.Vector3(xi, yi + w / 2.0, 0) for xi, yi in reversed(list(zip(x, y)))]
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

def simulate_switch(state="amorphous", resolution=20, l_c=1.50, gap=0.045, symmetric_pcm=True):
    lambda_0 = 1.064
    fcen = 1.0 / lambda_0
    df = 0.10 * fcen

    w_wg = 0.250
    h_wg = 0.220
    l_bend = 0.950
    y_pitch = 0.600
    t_pcm = 0.020
    core_index = 3.480

    clad_mat = mp.Medium(index=1.449)
    core_mat = mp.Medium(index=core_index)
    pcm_am = mp.Medium(index=2.70, D_conductivity=2 * math.pi * fcen * 1.0e-5)
    pcm_cr = mp.Medium(index=3.30, D_conductivity=2 * math.pi * fcen * 8.0e-5)

    l_in = 0.350
    l_out = 0.350
    iso_buf = 0.150

    y1 = (w_wg + gap) / 2.0
    y2_c = -y1
    y_out2 = -(y_pitch - y1)

    cell_len = l_in + l_c + l_bend + l_out
    cell_w = y_pitch + 2 * iso_buf
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

    # Reference waveguide
    geom_ref = [
        mp.Block(center=mp.Vector3(0, 0, -sz / 4.0), size=mp.Vector3(mp.inf, mp.inf, sz / 2.0), material=clad_mat),
        mp.Block(center=mp.Vector3(0, y1, h_wg / 2.0), size=mp.Vector3(mp.inf, w_wg, h_wg), material=core_mat)
    ]
    sim_ref = mp.Simulation(cell_size=cell, boundary_layers=pml_layers, geometry=geom_ref,
                            sources=sources, resolution=resolution, default_material=clad_mat)
    f_ref = sim_ref.add_flux(fcen, 0, 1, mp.FluxRegion(center=mp.Vector3(x_mon, y1, h_wg / 2.0), size=mp.Vector3(0, mon_w, mon_h)))
    sim_ref.run(until_after_sources=mp.stop_when_fields_decayed(20, mp.Hz, mp.Vector3(x_mon, y1, h_wg / 2.0), 1e-3))
    p_ref = max(mp.get_fluxes(f_ref)[0], 1e-12)

    # Device geometry
    geom_dev = [
        mp.Block(center=mp.Vector3(0, 0, -sz / 4.0), size=mp.Vector3(mp.inf, mp.inf, sz / 2.0), material=clad_mat),
        mp.Block(center=mp.Vector3(0, y1, h_wg / 2.0), size=mp.Vector3(mp.inf, w_wg, h_wg), material=core_mat),
    ]

    # WG2 entrance adiabatic taper
    taper_len = min(0.350, l_in)
    poly_wg2_in = make_ccw_linear_taper(x_c_start - taper_len, x_c_start, y2_c, 0.020, w_wg)
    geom_dev.append(create_prism_exact(poly_wg2_in, h_wg, 0.0, core_mat))

    # WG2 straight coupled section
    geom_dev.append(mp.Block(
        center=mp.Vector3((x_c_start + x_c_end) / 2.0, y2_c, h_wg / 2.0),
        size=mp.Vector3(l_c, w_wg, h_wg),
        material=core_mat
    ))

    # WG2 output S-bend fanout
    poly_wg2_bend = make_ccw_sbend(x_c_end, x_bend_end, y2_c, y_out2, w_wg)
    geom_dev.append(create_prism_exact(poly_wg2_bend, h_wg, 0.0, core_mat))

    # WG2 output lead
    geom_dev.append(mp.Block(
        center=mp.Vector3((x_bend_end + x_right) / 2.0, y_out2, h_wg / 2.0),
        size=mp.Vector3(x_right - x_bend_end, w_wg, h_wg),
        material=core_mat
    ))

    # Perfectly aligned adiabatic PCM patches strictly inside [x_c_start, x_c_end]
    mat_pcm1 = pcm_cr if state.lower() == "crystalline" else pcm_am
    l_tip = 0.200
    l_body = max(l_c - 2.0 * l_tip, 0.100)
    x_pcm_center = (x_c_start + x_c_end) / 2.0
    poly_pcm1 = make_ccw_pcm_patch(x_pcm_center, y1, l_body, l_tip, w_wg)
    geom_dev.append(create_prism_exact(poly_pcm1, t_pcm, h_wg, mat_pcm1))

    # Matching Amorphous PCM on WG2 for exact phase matching
    if symmetric_pcm:
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

    print(f"\n[{state.upper()}] (symmetric_pcm={symmetric_pcm}, L_c={l_c:.3f}um, Gap={gap*1000:.0f}nm):")
    print(f"  Port 1 (Through / Bar) : {p1*100:.2f}%")
    print(f"  Port 2 (Cross)         : {p2*100:.2f}%")
    print(f"  Total Guided Power     : {p_tot*100:.2f}% (Rad Loss: {(1-p_tot)*100:.2f}%)")
    print(f"  Signal Port Power      : {sig*100:.2f}% (IL = {il_db:.4f} dB)")
    print(f"  Leakage (Crosstalk)    : {leak*100:.2f}% (XT = {xt_db:.2f} dB)")
    print(f"  Cell Footprint         : {cell_area:.3f} um^2 (Length={cell_len:.2f}um, Width={cell_w:.2f}um)")

    return {"p1": p1, "p2": p2, "p_tot": p_tot, "sig": sig, "leak": leak, "il_db": il_db, "xt_db": xt_db, "cell_area": cell_area}

if __name__ == '__main__':
    st = sys.argv[1] if len(sys.argv) > 1 else "amorphous"
    res = int(sys.argv[2]) if len(sys.argv) > 2 else 15
    lc = float(sys.argv[3]) if len(sys.argv) > 3 else 1.50
    sym = sys.argv[4].lower() == "true" if len(sys.argv) > 4 else True
    simulate_switch(state=st, resolution=res, l_c=lc, symmetric_pcm=sym)
