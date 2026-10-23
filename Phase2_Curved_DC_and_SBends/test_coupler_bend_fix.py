"""
Test fast-divergence angled bend to suppress fringe coupling in 1x2 switch
"""
import sys, math
import numpy as np
import meep as mp

def test_fast_bend(angle_deg=10.0, l_c=4.65, gap=0.040, state="crystalline"):
    lambda_0 = 1.064
    fcen = 1.0 / lambda_0
    df = 0.1 * fcen
    w_wg = 0.240
    h_wg = 0.220
    t_pcm = 0.150

    y_wg1 = 0.000
    y_wg2_c = -(w_wg + gap) # -0.280 um
    y_out2 = -0.550          # 550 nm separation at output

    # Linear tilt angle from coupler end
    # dy = y_out2 - y_wg2_c = -0.270 um
    # with angle theta = 12 deg: dx_tilt = abs(dy) / tan(theta) = 0.270 / 0.212 = 1.27 um
    dx_tilt = abs(y_out2 - y_wg2_c) / math.tan(math.radians(angle_deg))
    l_bend = dx_tilt
    l_in, l_out = 0.25, 0.25

    cell_len = l_in + l_c + l_bend + l_out
    cell_w = abs(y_out2) + w_wg + 0.16
    sx = cell_len + 1.6
    sy = cell_w + 1.2
    sz = h_wg + t_pcm + 1.2
    cell = mp.Vector3(sx, sy, sz)
    pml = [mp.PML(0.5)]

    clad = mp.Medium(index=1.449)
    core = mp.Medium(index=3.480)
    pcm_am = mp.Medium(index=2.70)
    pcm_cr = mp.Medium(index=3.30)
    mat_pcm1 = pcm_cr if state == "crystalline" else pcm_am
    mat_pcm2 = pcm_am

    x_left = -sx / 2.0
    x_right = sx / 2.0
    x_min_cell = -cell_len / 2.0
    x_in_end = x_min_cell + l_in
    x_c_end = x_in_end + l_c
    x_bend_end = x_c_end + l_bend
    x_mon = x_bend_end + 0.12

    # Linear bend polygon for WG2
    def make_linear_bend(x0, x1, y0, y1, w):
        p1 = mp.Vector3(x0, y0 + w/2.0, 0)
        p2 = mp.Vector3(x1, y1 + w/2.0, 0)
        p3 = mp.Vector3(x1, y1 - w/2.0, 0)
        p4 = mp.Vector3(x0, y0 - w/2.0, 0)
        return [p1, p2, p3, p4]

    # Input sharp taper for WG2
    def make_taper(x0, x1, y, w_in, w_out):
        return [
            mp.Vector3(x0, y + w_in/2.0, 0),
            mp.Vector3(x1, y + w_out/2.0, 0),
            mp.Vector3(x1, y - w_out/2.0, 0),
            mp.Vector3(x0, y - w_in/2.0, 0),
        ]

    sources = [
        mp.EigenModeSource(src=mp.GaussianSource(fcen, fwidth=df),
                           center=mp.Vector3(x_left + 0.6, y_wg1, h_wg/2.0),
                           size=mp.Vector3(0, w_wg*3, h_wg*3), eig_band=1, direction=mp.X)
    ]

    # 1. Straight ref
    geom_ref = [
        mp.Block(center=mp.Vector3(0, 0, -sz/4.0), size=mp.Vector3(mp.inf, mp.inf, sz/2.0), material=clad),
        mp.Block(center=mp.Vector3(0, y_wg1, h_wg/2.0), size=mp.Vector3(mp.inf, w_wg, h_wg), material=core)
    ]
    sim_ref = mp.Simulation(cell_size=cell, boundary_layers=pml, geometry=geom_ref, sources=sources, resolution=20, default_material=clad)
    f_ref = sim_ref.add_flux(fcen, 0, 1, mp.FluxRegion(center=mp.Vector3(x_mon, y_wg1, h_wg/2.0), size=mp.Vector3(0, 0.26, 0.30)))
    sim_ref.run(until_after_sources=mp.stop_when_fields_decayed(20, mp.Hz, mp.Vector3(x_mon, y_wg1, h_wg/2.0), 1e-3))
    p_ref = max(mp.get_fluxes(f_ref)[0], 1e-12)

    # 2. Device
    geom = [
        mp.Block(center=mp.Vector3(0, 0, -sz/4.0), size=mp.Vector3(mp.inf, mp.inf, sz/2.0), material=clad),
        # WG1 (Straight)
        mp.Block(center=mp.Vector3(0, y_wg1, h_wg/2.0), size=mp.Vector3(mp.inf, w_wg, h_wg), material=core),
        # WG2
        mp.Prism(vertices=make_taper(x_in_end - 0.35, x_in_end, y_wg2_c, 0.025, w_wg), height=h_wg, axis=mp.Vector3(0,0,1), material=core),
        mp.Block(center=mp.Vector3((x_in_end + x_c_end)/2.0, y_wg2_c, h_wg/2.0), size=mp.Vector3(l_c, w_wg, h_wg), material=core),
        mp.Prism(vertices=make_linear_bend(x_c_end, x_bend_end, y_wg2_c, y_out2, w_wg), height=h_wg, axis=mp.Vector3(0,0,1), material=core),
        mp.Block(center=mp.Vector3((x_bend_end + x_right)/2.0, y_out2, h_wg/2.0), size=mp.Vector3(x_right - x_bend_end, w_wg, h_wg), material=core),
        # PCM on WG1
        mp.Block(center=mp.Vector3((x_in_end + x_c_end)/2.0, y_wg1, h_wg + t_pcm/2.0), size=mp.Vector3(l_c, w_wg, t_pcm), material=mat_pcm1),
        # PCM on WG2
        mp.Block(center=mp.Vector3((x_in_end + x_c_end)/2.0, y_wg2_c, h_wg + t_pcm/2.0), size=mp.Vector3(l_c, w_wg, t_pcm), material=mat_pcm2),
    ]
    sim = mp.Simulation(cell_size=cell, boundary_layers=pml, geometry=geom, sources=sources, resolution=20, default_material=clad)
    f1 = sim.add_flux(fcen, 0, 1, mp.FluxRegion(center=mp.Vector3(x_mon, y_wg1, h_wg/2.0), size=mp.Vector3(0, 0.26, 0.30)))
    f2 = sim.add_flux(fcen, 0, 1, mp.FluxRegion(center=mp.Vector3(x_mon, y_out2, h_wg/2.0), size=mp.Vector3(0, 0.26, 0.30)))

    sim.run(until_after_sources=mp.stop_when_fields_decayed(20, mp.Hz, mp.Vector3(x_mon, y_wg1, h_wg/2.0), 1e-3))
    p1 = max(mp.get_fluxes(f1)[0] / p_ref, 1e-12)
    p2 = max(mp.get_fluxes(f2)[0] / p_ref, 1e-12)

    sig = p1 if state == "crystalline" else p2
    leak = p2 if state == "crystalline" else p1
    il = -10 * math.log10(max(sig, 1e-6))
    xt = 10 * math.log10(max(leak / max(sig, 1e-6), 1e-6))
    print(f"[{state.upper()}] Angle = {angle_deg} deg: P1(Bar) = {p1*100:.2f}%, P2(Cross) = {p2*100:.2f}%, Total = {(p1+p2)*100:.2f}%, IL = {il:.3f} dB, XT = {xt:.2f} dB")
    return p1, p2, il, xt

if __name__ == '__main__':
    test_fast_bend(angle_deg=12.0, state="crystalline")
