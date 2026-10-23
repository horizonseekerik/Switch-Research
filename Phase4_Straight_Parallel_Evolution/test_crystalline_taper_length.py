import sys, math, time
import numpy as np
import meep as mp

def make_ccw_linear_taper(x_start, x_end, y, w_start, w_end, n_pts=30):
    u = np.linspace(0, 1, n_pts)
    x = x_start + u * (x_end - x_start)
    w = w_start + u * (w_end - w_start)
    bot = [mp.Vector3(xi, y - wi / 2.0, 0) for xi, wi in zip(x, w)]
    top = [mp.Vector3(xi, y + wi / 2.0, 0) for xi, wi in reversed(list(zip(x, w)))]
    return bot + top

def create_prism_exact(vertices, height, z_base, material):
    cx = float(np.mean([pt.x for pt in vertices]))
    cy = float(np.mean([pt.y for pt in vertices]))
    cz = float(z_base + height / 2.0)
    return mp.Prism(vertices=vertices, height=height, axis=mp.Vector3(0, 0, 1),
                    center=mp.Vector3(cx, cy, cz), material=material)

def test_taper_length(l_tip=1.0000000000, t_pcm=0.1600000000, resolution=20):
    lambda_0 = 1.0640000000
    fcen = 1.0 / lambda_0
    df = 0.10 * fcen

    w_wg = 0.2400000000
    h_wg = 0.2200000000
    l_body = 2.0000000000

    clad_mat = mp.Medium(index=1.4490000000)
    core_mat = mp.Medium(index=3.4800000000)
    pcm_cr = mp.Medium(index=3.3000000000)
    pcm_am = mp.Medium(index=2.7000000000)

    cell_len = 2 * l_tip + l_body + 1.0000000000
    cell_w = 1.6000000000
    dpml = 0.5000000000
    pad_y = 0.4000000000
    pad_z = 0.4500000000

    sx = cell_len + 2 * dpml
    sy = cell_w + 2 * (dpml + pad_y)
    sz = h_wg + t_pcm + 2 * (dpml + pad_z)

    cell = mp.Vector3(sx, sy, sz)
    pml_layers = [mp.PML(dpml)]

    x_left = -sx / 2.0
    x0 = -cell_len / 2.0
    x_pcm_start = x0 + 0.5000000000
    x_body_start = x_pcm_start + l_tip
    x_body_end = x_body_start + l_body
    x_pcm_end = x_body_end + l_tip
    x_mon = x_pcm_end + 0.2500000000

    src_pos_x = x_left + dpml + 0.1500000000
    sources = [
        mp.EigenModeSource(
            src=mp.GaussianSource(fcen, fwidth=df),
            center=mp.Vector3(src_pos_x, 0, h_wg / 2.0),
            size=mp.Vector3(0, w_wg * 3.0, (h_wg + t_pcm) * 3.0),
            eig_band=1,
            direction=mp.X
        )
    ]

    # Reference
    geom_ref = [
        mp.Block(center=mp.Vector3(0, 0, -sz / 4.0), size=mp.Vector3(mp.inf, mp.inf, sz / 2.0), material=clad_mat),
        mp.Block(center=mp.Vector3(0, 0, h_wg / 2.0), size=mp.Vector3(mp.inf, w_wg, h_wg), material=core_mat)
    ]
    sim_ref = mp.Simulation(cell_size=cell, boundary_layers=pml_layers, geometry=geom_ref,
                            sources=sources, resolution=resolution, default_material=clad_mat)
    f_ref = sim_ref.add_flux(fcen, 0, 1, mp.FluxRegion(center=mp.Vector3(x_mon, 0, h_wg / 2.0), size=mp.Vector3(0, 0.6, 0.8)))
    sim_ref.run(until_after_sources=mp.stop_when_fields_decayed(20, mp.Hz, mp.Vector3(x_mon, 0, h_wg / 2.0), 1e-4))
    res_ref = sim_ref.get_eigenmode_coefficients(f_ref, [1], eig_parity=mp.NO_PARITY)
    p_ref = abs(res_ref.alpha[0, 0, 0])**2

    # Devices: test crystalline and amorphous
    for state, pcm_mat in [("crystalline", pcm_cr), ("amorphous", pcm_am)]:
        geom_dev = [
            mp.Block(center=mp.Vector3(0, 0, -sz / 4.0), size=mp.Vector3(mp.inf, mp.inf, sz / 2.0), material=clad_mat),
            mp.Block(center=mp.Vector3(0, 0, h_wg / 2.0), size=mp.Vector3(mp.inf, w_wg, h_wg), material=core_mat),
        ]
        poly_in = make_ccw_linear_taper(x_pcm_start, x_body_start, 0, 0.0150000000, w_wg)
        geom_dev.append(create_prism_exact(poly_in, t_pcm, h_wg, pcm_mat))
        geom_dev.append(mp.Block(center=mp.Vector3((x_body_start + x_body_end) / 2.0, 0, h_wg + t_pcm / 2.0),
                                 size=mp.Vector3(l_body, w_wg, t_pcm), material=pcm_mat))
        poly_out = make_ccw_linear_taper(x_body_end, x_pcm_end, 0, w_wg, 0.0150000000)
        geom_dev.append(create_prism_exact(poly_out, t_pcm, h_wg, pcm_mat))

        sim_dev = mp.Simulation(cell_size=cell, boundary_layers=pml_layers, geometry=geom_dev,
                                sources=sources, resolution=resolution, default_material=clad_mat)
        f_dev = sim_dev.add_flux(fcen, 0, 1, mp.FluxRegion(center=mp.Vector3(x_mon, 0, h_wg / 2.0), size=mp.Vector3(0, 0.6, 0.8)))
        sim_dev.run(until_after_sources=mp.stop_when_fields_decayed(20, mp.Hz, mp.Vector3(x_mon, 0, h_wg / 2.0), 1e-4))
        res_dev = sim_dev.get_eigenmode_coefficients(f_dev, [1], eig_parity=mp.NO_PARITY)
        p_dev = abs(res_dev.alpha[0, 0, 0])**2
        trans = p_dev / p_ref
        loss = 1.0 - trans
        print(f"L_tip = {l_tip:6.3f} um [{state.upper():<11}] Trans = {trans*100:15.10f}% (Loss = {loss*100:15.10f}%)")
        sys.stdout.flush()

if __name__ == "__main__":
    t0 = time.time()
    print("=" * 85)
    print("PCM TAPER TRANSMISSION VS TAPER LENGTH (t_pcm = 160 nm, Resolution = 20)")
    print("=" * 85)
    for lt in [1.000, 1.250, 1.500]:
        test_taper_length(l_tip=lt)
    print("=" * 85)
    print(f"Total runtime = {time.time() - t0:.2f} s")
