"""
Test pure S-bend transmission for 240 nm silicon waveguide at 1064 nm in 3D MEEP
"""
import sys, math
import numpy as np
import meep as mp

def run_sbend_test(l_bend=1.100, dy=0.110, resolution=20):
    lambda_0 = 1.064
    fcen = 1.0 / lambda_0
    df = 0.1 * fcen
    w_wg = 0.240
    h_wg = 0.220
    l_in = 0.50
    l_out = 0.80

    cell_len = l_in + l_bend + l_out
    cell_w = 1.20
    dpml = 0.50
    pad = 0.35
    sx = cell_len + 2 * (dpml + pad)
    sy = cell_w + 2 * (dpml + pad)
    sz = h_wg + 2 * (dpml + pad)
    cell = mp.Vector3(sx, sy, sz)
    pml_layers = [mp.PML(dpml)]

    clad_mat = mp.Medium(index=1.449)
    core_mat = mp.Medium(index=3.480)

    x_left = -sx / 2.0
    x_in_end = -cell_len / 2.0 + l_in
    x_bend_end = x_in_end + l_bend
    x_mon = x_bend_end + 0.40

    mon_w = 0.40
    mon_h = 0.40

    sources = [
        mp.EigenModeSource(
            src=mp.GaussianSource(fcen, fwidth=df),
            center=mp.Vector3(x_left + dpml + 0.10, 0, h_wg / 2.0),
            size=mp.Vector3(0, w_wg * 3.0, h_wg * 3.0),
            eig_band=1,
            direction=mp.X
        )
    ]

    # Reference straight waveguide
    geom_ref = [
        mp.Block(center=mp.Vector3(0, 0, -sz/4.0), size=mp.Vector3(mp.inf, mp.inf, sz/2.0), material=clad_mat),
        mp.Block(center=mp.Vector3(0, 0, h_wg/2.0), size=mp.Vector3(mp.inf, w_wg, h_wg), material=core_mat)
    ]
    sim_ref = mp.Simulation(cell_size=cell, boundary_layers=pml_layers, geometry=geom_ref,
                            sources=sources, resolution=resolution, default_material=clad_mat)
    f_ref = sim_ref.add_flux(fcen, 0, 1, mp.FluxRegion(center=mp.Vector3(x_mon, 0, h_wg/2.0), size=mp.Vector3(0, mon_w, mon_h)))
    sim_ref.run(until_after_sources=mp.stop_when_fields_decayed(20, mp.Hz, mp.Vector3(x_mon, 0, h_wg/2.0), 1e-3))
    p_ref = max(mp.get_fluxes(f_ref)[0], 1e-12)

    # S-bend waveguide
    def make_cosine_sbend(x_start, x_end, y_start, y_end, w, n_pts=30):
        dx = x_end - x_start
        dy_val = y_end - y_start
        u = np.linspace(0, 1, n_pts)
        x = x_start + u * dx
        y = y_start + dy_val * (0.5 * (1.0 - np.cos(np.pi * u)))
        top = [mp.Vector3(xi, yi + w / 2.0, 0) for xi, yi in zip(x, y)]
        bot = [mp.Vector3(xi, yi - w / 2.0, 0) for xi, yi in reversed(list(zip(x, y)))]
        return top + bot

    geom_bend = [
        mp.Block(center=mp.Vector3(0, 0, -sz/4.0), size=mp.Vector3(mp.inf, mp.inf, sz/2.0), material=clad_mat),
        mp.Block(center=mp.Vector3((x_left + x_in_end)/2.0, 0, h_wg/2.0), size=mp.Vector3(x_in_end - x_left, w_wg, h_wg), material=core_mat),
        mp.Prism(vertices=make_cosine_sbend(x_in_end, x_bend_end, 0, dy, w_wg), height=h_wg, axis=mp.Vector3(0, 0, 1), material=core_mat),
        mp.Block(center=mp.Vector3((x_bend_end + sx/2.0)/2.0, dy, h_wg/2.0), size=mp.Vector3(sx/2.0 - x_bend_end, w_wg, h_wg), material=core_mat),
    ]
    sim_bend = mp.Simulation(cell_size=cell, boundary_layers=pml_layers, geometry=geom_bend,
                             sources=sources, resolution=resolution, default_material=clad_mat)
    f_bend = sim_bend.add_flux(fcen, 0, 1, mp.FluxRegion(center=mp.Vector3(x_mon, dy, h_wg/2.0), size=mp.Vector3(0, mon_w, mon_h)))
    sim_bend.run(until_after_sources=mp.stop_when_fields_decayed(20, mp.Hz, mp.Vector3(x_mon, dy, h_wg/2.0), 1e-3))
    p_bend = max(mp.get_fluxes(f_bend)[0], 1e-12)

    trans = p_bend / p_ref
    loss_db = -10.0 * math.log10(trans)
    print(f"S-bend (L={l_bend:.2f}um, dy={dy*1000:.0f}nm): Trans = {trans*100:.2f}%, Loss = {loss_db:.3f} dB")
    return trans

if __name__ == "__main__":
    for lb in [0.80, 1.00, 1.20, 1.50]:
        run_sbend_test(l_bend=lb)
