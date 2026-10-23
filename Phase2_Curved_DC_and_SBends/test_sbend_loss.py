import sys, math
import numpy as np
import meep as mp

def make_ccw_sbend(x_start, x_end, y_start, y_end, w1, w2, n_pts=30):
    dx = x_end - x_start
    dy = y_end - y_start
    u = np.linspace(0, 1, n_pts)
    x = x_start + u * dx
    y = y_start + dy * (3 * u**2 - 2 * u**3)
    w = w1 + (w2 - w1) * (3 * u**2 - 2 * u**3)
    bot = [mp.Vector3(xi, yi - wi / 2.0, 0) for xi, yi, wi in zip(x, y, w)]
    top = [mp.Vector3(xi, yi + wi / 2.0, 0) for xi, yi, wi in reversed(list(zip(x, y, w)))]
    return bot + top

def create_prism_exact(vertices, height, z_base, material):
    cx = float(np.mean([pt.x for pt in vertices]))
    cy = float(np.mean([pt.y for pt in vertices]))
    cz = float(z_base + height / 2.0)
    return mp.Prism(vertices=vertices, height=height, axis=mp.Vector3(0, 0, 1),
                    center=mp.Vector3(cx, cy, cz), material=material)

def test_sbend(l_bend=0.75, dy=0.2045, res=15):
    lambda_0 = 1.064
    fcen = 1.0 / lambda_0
    w_wg = 0.240
    h_wg = 0.220
    core_mat = mp.Medium(index=3.480)
    clad_mat = mp.Medium(index=1.449)

    l_lead = 0.50
    sx = l_lead + l_bend + l_lead + 1.0
    sy = abs(dy) + 1.2
    sz = h_wg + 1.0
    cell = mp.Vector3(sx, sy, sz)
    pml = [mp.PML(0.5)]

    x_start = -l_bend / 2.0
    x_end = l_bend / 2.0
    x_src = -sx/2.0 + 0.6
    x_mon = sx/2.0 - 0.6

    # Straight reference
    geom_ref = [
        mp.Block(center=mp.Vector3(0, 0, 0), size=mp.Vector3(mp.inf, w_wg, h_wg), material=core_mat)
    ]
    sim_ref = mp.Simulation(cell_size=cell, boundary_layers=pml, geometry=geom_ref, resolution=res, default_material=clad_mat)
    src = [mp.EigenModeSource(src=mp.GaussianSource(fcen, fwidth=0.1*fcen),
                              center=mp.Vector3(x_src, 0, 0), size=mp.Vector3(0, 1.0, 1.0), eig_band=1, direction=mp.X)]
    sim_ref.sources = src
    f_ref = sim_ref.add_flux(fcen, 0, 1, mp.FluxRegion(center=mp.Vector3(x_mon, 0, 0), size=mp.Vector3(0, 0.6, 0.6)))
    sim_ref.run(until_after_sources=mp.stop_when_fields_decayed(20, mp.Hz, mp.Vector3(x_mon, 0, 0), 1e-4))
    res_ref = sim_ref.get_eigenmode_coefficients(f_ref, [1])
    pref = abs(res_ref.alpha[0,0,0])**2

    # S-bend device
    poly = make_ccw_sbend(x_start, x_end, 0, dy, w_wg, w_wg)
    geom_bend = [
        mp.Block(center=mp.Vector3((x_src + x_start)/2.0, 0, 0), size=mp.Vector3(x_start - (-sx/2.0), w_wg, h_wg), material=core_mat),
        create_prism_exact(poly, h_wg, -h_wg/2.0, core_mat),
        mp.Block(center=mp.Vector3((x_end + sx/2.0)/2.0, dy, 0), size=mp.Vector3(sx/2.0 - x_end, w_wg, h_wg), material=core_mat)
    ]
    sim_bend = mp.Simulation(cell_size=cell, boundary_layers=pml, geometry=geom_bend, resolution=res, default_material=clad_mat)
    sim_bend.sources = src
    f_bend = sim_bend.add_flux(fcen, 0, 1, mp.FluxRegion(center=mp.Vector3(x_mon, dy, 0), size=mp.Vector3(0, 0.6, 0.6)))
    sim_bend.run(until_after_sources=mp.stop_when_fields_decayed(20, mp.Hz, mp.Vector3(x_mon, dy, 0), 1e-4))
    res_bend = sim_bend.get_eigenmode_coefficients(f_bend, [1])
    pbend = abs(res_bend.alpha[0,0,0])**2
    t = pbend / pref
    print(f"L_bend = {l_bend:.2f} um, dy = {dy*1e3:.1f} nm -> Transmission = {t*100:.6f}%, Loss = {(1-t)*100:.6f}%")

if __name__ == "__main__":
    for lb in [0.75, 1.00, 1.50, 2.00]:
        test_sbend(l_bend=lb, dy=0.2045)
