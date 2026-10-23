import sys, math
import numpy as np
import meep as mp

def make_ccw_parabolic_taper(x_start, x_end, y, w_start, w_end, n_pts=30, p=2.0):
    u = np.linspace(0, 1, n_pts)
    x = x_start + u * (x_end - x_start)
    w = w_start + (u**p) * (w_end - w_start)
    bot = [mp.Vector3(xi, y - wi / 2.0, 0) for xi, wi in zip(x, w)]
    top = [mp.Vector3(xi, y + wi / 2.0, 0) for xi, wi in reversed(list(zip(x, w)))]
    return bot + top

def create_prism_exact(vertices, height, z_base, material):
    cx = float(np.mean([pt.x for pt in vertices]))
    cy = float(np.mean([pt.y for pt in vertices]))
    cz = float(z_base + height / 2.0)
    return mp.Prism(vertices=vertices, height=height, axis=mp.Vector3(0, 0, 1),
                    center=mp.Vector3(cx, cy, cz), material=material)

def test_taper(profile="linear", l_tip=1.0, t_pcm=0.160, power=2.0, res=20):
    lambda_0 = 1.064
    fcen = 1.0 / lambda_0
    w_wg = 0.240
    h_wg = 0.220
    core_mat = mp.Medium(index=3.480)
    clad_mat = mp.Medium(index=1.449)
    pcm_cr = mp.Medium(index=3.300)

    l_body = 1.0
    l_lead = 0.4
    total_len = 2*l_lead + 2*l_tip + l_body
    dpml = 0.5
    sx = total_len + 2*dpml
    sy = 1.5
    sz = h_wg + t_pcm + 1.0
    cell = mp.Vector3(sx, sy, sz)
    pml = [mp.PML(dpml)]

    x_src = -sx/2.0 + dpml + 0.1
    x_mon = sx/2.0 - dpml - 0.1
    x1 = -l_body/2.0 - l_tip
    x2 = -l_body/2.0
    x3 = l_body/2.0
    x4 = l_body/2.0 + l_tip

    # Reference bare waveguide
    geom_ref = [mp.Block(center=mp.Vector3(0, 0, h_wg/2.0), size=mp.Vector3(mp.inf, w_wg, h_wg), material=core_mat)]
    src = [mp.EigenModeSource(src=mp.GaussianSource(fcen, fwidth=0.1*fcen),
                              center=mp.Vector3(x_src, 0, h_wg/2.0), size=mp.Vector3(0, 1.0, 1.0), eig_band=1, direction=mp.X)]
    sim_ref = mp.Simulation(cell_size=cell, boundary_layers=pml, geometry=geom_ref, resolution=res, default_material=clad_mat)
    sim_ref.sources = src
    f_ref = sim_ref.add_flux(fcen, 0, 1, mp.FluxRegion(center=mp.Vector3(x_mon, 0, h_wg/2.0), size=mp.Vector3(0, 0.8, 0.8)))
    sim_ref.run(until_after_sources=mp.stop_when_fields_decayed(20, mp.Hz, mp.Vector3(x_mon, 0, h_wg/2.0), 1e-4))
    res_ref = sim_ref.get_eigenmode_coefficients(f_ref, [1])
    pref = abs(res_ref.alpha[0,0,0])**2

    # Tapered device
    geom_dev = [mp.Block(center=mp.Vector3(0, 0, h_wg/2.0), size=mp.Vector3(mp.inf, w_wg, h_wg), material=core_mat)]
    # In taper
    poly_in = make_ccw_parabolic_taper(x1, x2, 0, 0.015, w_wg, p=power if profile=="parabolic" else 1.0)
    geom_dev.append(create_prism_exact(poly_in, t_pcm, h_wg, pcm_cr))
    # Body
    geom_dev.append(mp.Block(center=mp.Vector3(0, 0, h_wg + t_pcm/2.0), size=mp.Vector3(l_body, w_wg, t_pcm), material=pcm_cr))
    # Out taper
    poly_out = make_ccw_parabolic_taper(x3, x4, 0, w_wg, 0.015, p=1.0/power if profile=="parabolic" else 1.0)
    geom_dev.append(create_prism_exact(poly_out, t_pcm, h_wg, pcm_cr))

    sim_dev = mp.Simulation(cell_size=cell, boundary_layers=pml, geometry=geom_dev, resolution=res, default_material=clad_mat)
    sim_dev.sources = src
    f_dev = sim_dev.add_flux(fcen, 0, 1, mp.FluxRegion(center=mp.Vector3(x_mon, 0, h_wg/2.0), size=mp.Vector3(0, 0.8, 0.8)))
    sim_dev.run(until_after_sources=mp.stop_when_fields_decayed(20, mp.Hz, mp.Vector3(x_mon, 0, h_wg/2.0), 1e-4))
    res_dev = sim_dev.get_eigenmode_coefficients(f_dev, [1])
    pdev = abs(res_dev.alpha[0,0,0])**2
    t = pdev / pref
    print(f"Profile: {profile:<10} L_tip: {l_tip:.2f} um, t_pcm: {t_pcm*1e3:.0f} nm -> Trans = {t*100:12.8f}%, Loss = {(1-t)*100:12.8f}%")
    sys.stdout.flush()

if __name__ == "__main__":
    for tp in [0.120, 0.140, 0.160]:
        for lt in [0.80, 1.00, 1.20]:
            test_taper(profile="linear", l_tip=lt, t_pcm=tp)
            test_taper(profile="parabolic", l_tip=lt, t_pcm=tp, power=2.0)
