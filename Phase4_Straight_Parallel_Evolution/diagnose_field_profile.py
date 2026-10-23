import sys, os, math, time
import numpy as np
import meep as mp

def run_diagnostic(state="crystalline", resolution=15, l_c=3.80, gap=0.100, t_pcm=0.200):
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
    l_bend = 1.000
    y_pitch = 0.450
    iso_buf = 0.080
    l_tip = 0.250
    l_hold = 0.350

    y1 = 0.0
    y2_c = -(w_wg + gap)
    y_out2 = -y_pitch

    cell_len = l_in + l_c + l_bend + l_out
    cell_w = y_pitch + w_wg + 2 * iso_buf

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
    geom_ref = [
        mp.Block(center=mp.Vector3(0, 0, -sz / 4.0), size=mp.Vector3(mp.inf, mp.inf, sz / 2.0), material=clad_mat),
        mp.Block(center=mp.Vector3(0, y1, h_wg / 2.0), size=mp.Vector3(mp.inf, w_wg, h_wg), material=core_mat)
    ]
    sim_ref = mp.Simulation(cell_size=cell, boundary_layers=pml_layers, geometry=geom_ref,
                            sources=sources, resolution=resolution, default_material=clad_mat)
    f_ref = sim_ref.add_flux(fcen, 0, 1, mp.FluxRegion(center=mp.Vector3(x_bend_end + 0.15, y1, h_wg / 2.0), size=mp.Vector3(0, mon_w, mon_h)))
    sim_ref.run(until_after_sources=mp.stop_when_fields_decayed(20, mp.Hz, mp.Vector3(x_bend_end + 0.15, y1, h_wg / 2.0), 1e-3))
    p_ref = max(mp.get_fluxes(f_ref)[0], 1e-12)

    # Helper polygons
    def make_sbend(x1, x2, y1_val, y2_val, w):
        u = np.linspace(0, 1, 35)
        x = x1 + u * (x2 - x1)
        y = y1_val + (y2_val - y1_val) * (3 * u**2 - 2 * u**3)
        bot = [mp.Vector3(xi, yi - w / 2.0, 0) for xi, yi in zip(x, y)]
        top = [mp.Vector3(xi, yi + w / 2.0, 0) for xi, yi in reversed(list(zip(x, y)))]
        return bot + top

    def make_taper(x1, x2, y, w1, w2):
        u = np.linspace(0, 1, 20)
        x = x1 + u * (x2 - x1)
        w = w1 + u * (w2 - w1)
        bot = [mp.Vector3(xi, y - wi / 2.0, 0) for xi, wi in zip(x, w)]
        top = [mp.Vector3(xi, y + wi / 2.0, 0) for xi, wi in reversed(list(zip(x, w)))]
        return bot + top

    def make_prism(vertices, height, z_base, material):
        cx = float(np.mean([pt.x for pt in vertices]))
        cy = float(np.mean([pt.y for pt in vertices]))
        cz = float(z_base + height / 2.0)
        return mp.Prism(vertices=vertices, height=height, axis=mp.Vector3(0, 0, 1),
                        center=mp.Vector3(cx, cy, cz), material=material)

    geom_dev = [
        mp.Block(center=mp.Vector3(0, 0, -sz / 4.0), size=mp.Vector3(mp.inf, mp.inf, sz / 2.0), material=clad_mat),
        mp.Block(center=mp.Vector3(0, y1, h_wg / 2.0), size=mp.Vector3(mp.inf, w_wg, h_wg), material=core_mat),
        make_prism(make_taper(x_c_start - l_in, x_c_start, y2_c, 0.020, w_wg), h_wg, 0.0, core_mat),
        mp.Block(center=mp.Vector3((x_c_start + x_c_end)/2.0, y2_c, h_wg/2.0), size=mp.Vector3(l_c, w_wg, h_wg), material=core_mat),
        make_prism(make_sbend(x_c_end, x_bend_end, y2_c, y_out2, w_wg), h_wg, 0.0, core_mat),
        mp.Block(center=mp.Vector3((x_bend_end + x_right)/2.0, y_out2, h_wg/2.0), size=mp.Vector3(x_right - x_bend_end, w_wg, h_wg), material=core_mat),
    ]

    mat_pcm1 = pcm_cr if state.lower() == "crystalline" else pcm_am
    mat_pcm2 = pcm_am

    # WG1 PCM
    geom_dev.append(make_prism(make_taper(x_c_start, x_c_start + l_tip, y1, 0.030, w_wg), t_pcm, h_wg, mat_pcm1))
    geom_dev.append(mp.Block(center=mp.Vector3((x_c_start + l_tip + x_c_end + l_hold)/2.0, y1, h_wg + t_pcm/2.0),
                             size=mp.Vector3((x_c_end + l_hold) - (x_c_start + l_tip), w_wg, t_pcm), material=mat_pcm1))
    geom_dev.append(make_prism(make_taper(x_c_end + l_hold, x_c_end + l_hold + l_tip, y1, w_wg, 0.030), t_pcm, h_wg, mat_pcm1))

    # WG2 PCM
    geom_dev.append(make_prism(make_taper(x_c_start, x_c_start + l_tip, y2_c, 0.030, w_wg), t_pcm, h_wg, mat_pcm2))
    geom_dev.append(mp.Block(center=mp.Vector3((x_c_start + l_tip + x_c_end - l_tip)/2.0, y2_c, h_wg + t_pcm/2.0),
                             size=mp.Vector3((x_c_end - l_tip) - (x_c_start + l_tip), w_wg, t_pcm), material=mat_pcm2))
    geom_dev.append(make_prism(make_taper(x_c_end - l_tip, x_c_end, y2_c, w_wg, 0.030), t_pcm, h_wg, mat_pcm2))

    sim_dev = mp.Simulation(cell_size=cell, boundary_layers=pml_layers, geometry=geom_dev,
                            sources=sources, resolution=resolution, default_material=clad_mat)

    # Add 7 longitudinal cross-section monitors along the device
    x_checkpoints = [
        ("Input Entrance", x_c_start),
        ("25% Parallel", x_c_start + 0.25 * l_c),
        ("50% Parallel", x_c_start + 0.50 * l_c),
        ("75% Parallel", x_c_start + 0.75 * l_c),
        ("Parallel Exit", x_c_end),
        ("Mid S-Bend", x_c_end + 0.50 * l_bend),
        ("Output Ports", x_bend_end + 0.15)
    ]

    monitors = []
    for name, x_pos in x_checkpoints:
        # Determine y position of WG2 at this x
        if x_pos <= x_c_end:
            y2_pos = y2_c
        elif x_pos >= x_bend_end:
            y2_pos = y_out2
        else:
            u = (x_pos - x_c_end) / l_bend
            y2_pos = y2_c + (y_out2 - y2_c) * (3 * u**2 - 2 * u**3)
        
        f1 = sim_dev.add_flux(fcen, 0, 1, mp.FluxRegion(center=mp.Vector3(x_pos, y1, h_wg / 2.0), size=mp.Vector3(0, mon_w, mon_h)))
        f2 = sim_dev.add_flux(fcen, 0, 1, mp.FluxRegion(center=mp.Vector3(x_pos, y2_pos, h_wg / 2.0), size=mp.Vector3(0, mon_w, mon_h)))
        monitors.append((name, x_pos, f1, f2))

    y_dest = y_out2 if state.lower() == "amorphous" else y1
    sim_dev.run(until_after_sources=mp.stop_when_fields_decayed(20, mp.Hz, mp.Vector3(x_bend_end + 0.15, y_dest, h_wg / 2.0), 1e-3))

    print("\n" + "=" * 90)
    print(f"  LONGITUDINAL POWER FLOW DIAGNOSTIC: [{state.upper()}] (L_c = {l_c:.2f} um, Gap = {gap*1000:.0f} nm)")
    print("=" * 90)
    print(f"{'Checkpoint':<18} {'x(um)':<8} {'P_WG1(%)':<14} {'P_WG2(%)':<14} {'P_Total(%)':<14} {'Rad Loss(%)':<14}")
    print("-" * 90)

    for name, x_pos, f1, f2 in monitors:
        p1 = max(mp.get_fluxes(f1)[0] / p_ref, 1e-12) * 100.0
        p2 = max(mp.get_fluxes(f2)[0] / p_ref, 1e-12) * 100.0
        p_tot = p1 + p2
        rad_loss = max(0.0, 100.0 - p_tot)
        print(f"{name:<18} {x_pos:<8.3f} {p1:<14.2f} {p2:<14.2f} {p_tot:<14.2f} {rad_loss:<14.2f}")

    print("=" * 90 + "\n")

if __name__ == '__main__':
    run_diagnostic("crystalline", l_c=3.80, gap=0.100, t_pcm=0.200)
    run_diagnostic("amorphous", l_c=3.80, gap=0.100, t_pcm=0.200)
