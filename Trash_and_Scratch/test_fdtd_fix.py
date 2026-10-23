import math
import time
import meep as mp
import numpy as np

def run_test():
    lambda_0 = 1.064
    fcen = 1.0 / lambda_0
    w_wg = 0.250
    h_wg = 0.220
    t_pcm = 0.020
    gap = 0.045
    l_c = 2.400
    l_ext = 0.500
    l_bend = 0.950
    y_pitch_out = 0.600
    l_in = 0.350
    l_out = 0.350
    resolution = 20

    y_sep_coupler = w_wg + gap
    y_wg1 = y_sep_coupler / 2.0
    y_wg2 = -y_sep_coupler / 2.0
    y_out1 = y_pitch_out / 2.0
    y_out2 = -y_pitch_out / 2.0

    cell_len = l_in + l_c + l_bend + l_out
    cell_w = y_pitch_out + 2 * 0.200

    dpml = 0.50
    pad_x = 0.40
    pad_y = 0.40
    pad_z = 0.40
    sx = cell_len + 2 * (dpml + pad_x)
    sy = cell_w + 2 * (dpml + pad_y)
    sz = h_wg + t_pcm + 2 * (dpml + pad_z)

    cell = mp.Vector3(sx, sy, sz)
    pml_layers = [mp.PML(dpml)]

    x_left = -sx / 2.0
    x_right = sx / 2.0
    x_min_cell = -cell_len / 2.0
    x_in_end = x_min_cell + l_in
    x_coupler_end = x_in_end + l_c
    x_bend_end = x_coupler_end + l_bend
    x_max_cell = x_bend_end + l_out

    clad_mat = mp.Medium(index=1.449)
    core_mat = mp.Medium(index=3.480)

    # 1. Reference Simulation with continuous waveguide
    mon_w = w_wg * 2.5
    mon_h = h_wg * 2.5
    x_mon = x_bend_end + 0.15

    df = 0.1 * fcen
    sources = [
        mp.EigenModeSource(
            src=mp.GaussianSource(fcen, fwidth=df),
            center=mp.Vector3(x_left + dpml + 0.10, y_wg1, h_wg / 2.0),
            size=mp.Vector3(0, w_wg * 3.0, h_wg * 3.0),
            eig_band=1,
            direction=mp.X
        )
    ]

    print("Running reference straight waveguide...")
    sim_ref = mp.Simulation(
        cell_size=cell,
        boundary_layers=pml_layers,
        geometry=[
            mp.Block(center=mp.Vector3(0, 0, -sz/4.0), size=mp.Vector3(mp.inf, mp.inf, sz/2.0), material=clad_mat),
            mp.Block(center=mp.Vector3(0, y_wg1, h_wg/2.0), size=mp.Vector3(mp.inf, w_wg, h_wg), material=core_mat)
        ],
        sources=sources,
        resolution=resolution,
        default_material=clad_mat
    )
    flux_ref = sim_ref.add_flux(fcen, 0, 1, mp.FluxRegion(center=mp.Vector3(x_mon, y_wg1, h_wg / 2.0), size=mp.Vector3(0, mon_w, mon_h)))
    sim_ref.run(until_after_sources=mp.stop_when_fields_decayed(20, mp.Hz, mp.Vector3(x_mon, y_wg1, h_wg / 2.0), 1e-3))
    p_ref = mp.get_fluxes(flux_ref)[0]
    print(f"Reference Power = {p_ref:.6e}")

    # 2. Device simulation with waveguides extending into PML
    def make_sbend(x_start, x_end, y_start, y_end, w, n_pts=20):
        dx = x_end - x_start
        dy = y_end - y_start
        u = np.linspace(0, 1, n_pts)
        x = x_start + u * dx
        y = y_start + dy * (3 * u**2 - 2 * u**3)
        top = [mp.Vector3(xi, yi + w / 2.0, 0) for xi, yi in zip(x, y)]
        bot = [mp.Vector3(xi, yi - w / 2.0, 0) for xi, yi in reversed(list(zip(x, y)))]
        return top + bot

    geom = [
        mp.Block(center=mp.Vector3(0, 0, -sz/4.0), size=mp.Vector3(mp.inf, mp.inf, sz/2.0), material=clad_mat),
        # Input lead extending all the way through left PML
        mp.Block(center=mp.Vector3((x_left + x_in_end)/2.0, y_wg1, h_wg/2.0), size=mp.Vector3(x_in_end - x_left, w_wg, h_wg), material=core_mat),
        # Coupler section
        mp.Block(center=mp.Vector3((x_in_end + x_coupler_end)/2.0, y_wg1, h_wg/2.0), size=mp.Vector3(l_c, w_wg, h_wg), material=core_mat),
        mp.Block(center=mp.Vector3((x_in_end + x_coupler_end)/2.0, y_wg2, h_wg/2.0), size=mp.Vector3(l_c, w_wg, h_wg), material=core_mat),
        # S-bends
        mp.Prism(vertices=make_sbend(x_coupler_end, x_bend_end, y_wg1, y_out1, w_wg), height=h_wg, axis=mp.Vector3(0,0,1), material=core_mat),
        mp.Prism(vertices=make_sbend(x_coupler_end, x_bend_end, y_wg2, y_out2, w_wg), height=h_wg, axis=mp.Vector3(0,0,1), material=core_mat),
        # Output leads extending all the way through right PML
        mp.Block(center=mp.Vector3((x_bend_end + x_right)/2.0, y_out1, h_wg/2.0), size=mp.Vector3(x_right - x_bend_end, w_wg, h_wg), material=core_mat),
        mp.Block(center=mp.Vector3((x_bend_end + x_right)/2.0, y_out2, h_wg/2.0), size=mp.Vector3(x_right - x_bend_end, w_wg, h_wg), material=core_mat),
    ]

    print("\nRunning device simulation with PML-extended leads...")
    sim_dev = mp.Simulation(
        cell_size=cell,
        boundary_layers=pml_layers,
        geometry=geom,
        sources=sources,
        resolution=resolution,
        default_material=clad_mat
    )
    flux_p1 = sim_dev.add_flux(fcen, 0, 1, mp.FluxRegion(center=mp.Vector3(x_mon, y_out1, h_wg / 2.0), size=mp.Vector3(0, mon_w, mon_h)))
    flux_p2 = sim_dev.add_flux(fcen, 0, 1, mp.FluxRegion(center=mp.Vector3(x_mon, y_out2, h_wg / 2.0), size=mp.Vector3(0, mon_w, mon_h)))

    sim_dev.run(until_after_sources=mp.stop_when_fields_decayed(20, mp.Hz, mp.Vector3(x_mon, y_out2, h_wg / 2.0), 1e-3))
    p1 = mp.get_fluxes(flux_p1)[0]
    p2 = mp.get_fluxes(flux_p2)[0]
    print(f"\nRESULTS:")
    print(f"P1 (Bar)   = {p1:.6e} (Normalized = {p1/p_ref*100:.2f}%)")
    print(f"P2 (Cross) = {p2:.6e} (Normalized = {p2/p_ref*100:.2f}%)")
    print(f"Total Transmission = {(p1+p2)/p_ref*100:.2f}%")

if __name__ == "__main__":
    run_test()
