"""
MEASURE ISOLATED PATCH TRANSMISSION VS T_PCM IN 3D FDTD
======================================================
Tests:
- t_pcm in [120, 140, 160, 237.5] nm
- Measures fundamental mode transmission in both amorphous and crystalline states
"""

import sys, math, time
import numpy as np
import meep as mp

def make_ccw_linear_taper(x_start, x_end, y, w_start, w_end, n_pts=25):
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

def test_isolated_thickness(t_pcm_nm=120.0, l_tip=0.800, resolution=15):
    t_pcm = t_pcm_nm * 1e-3
    lambda_0 = 1.064
    fcen = 1.0 / lambda_0
    df = 0.10 * fcen

    w_wg = 0.240
    h_wg = 0.220
    l_body = 4.000

    clad_mat = mp.Medium(index=1.449)
    core_mat = mp.Medium(index=3.480)
    pcm_am = mp.Medium(index=2.70)
    pcm_cr = mp.Medium(index=3.30)

    cell_len = l_tip + l_body + l_tip + 1.0
    cell_w = 1.5
    dpml = 0.50
    pad = 0.35
    sx = cell_len + 2 * (dpml + pad)
    sy = cell_w + 2 * (dpml + pad)
    sz = h_wg + 0.250 + 2 * (dpml + pad)

    cell = mp.Vector3(sx, sy, sz)
    pml_layers = [mp.PML(dpml)]

    x_left = -sx / 2.0
    x0 = -cell_len / 2.0
    x_in_end = x0 + 0.5
    x_pcm_start = x_in_end
    x_pcm_body_start = x_pcm_start + l_tip
    x_pcm_body_end = x_pcm_body_start + l_body
    x_pcm_end = x_pcm_body_end + l_tip
    x_mon = x_pcm_end + 0.25

    src_pos_x = x_left + dpml + 0.10
    sources = [
        mp.EigenModeSource(
            src=mp.GaussianSource(fcen, fwidth=df),
            center=mp.Vector3(src_pos_x, 0, h_wg / 2.0),
            size=mp.Vector3(0, w_wg * 3.0, (h_wg + 0.250) * 3.0),
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
    f_ref = sim_ref.add_flux(fcen, 0, 1, mp.FluxRegion(center=mp.Vector3(x_mon, 0, h_wg / 2.0), size=mp.Vector3(0, 0.4, 0.8)))
    sim_ref.run(until_after_sources=mp.stop_when_fields_decayed(20, mp.Hz, mp.Vector3(x_mon, 0, h_wg / 2.0), 1e-4))
    ref_mode = max(abs(sim_ref.get_eigenmode_coefficients(f_ref, [1], eig_parity=mp.NO_PARITY).alpha[0, 0, 0])**2, 1e-15)

    results = {}
    for state, pcm_mat in [("amorphous", pcm_am), ("crystalline", pcm_cr)]:
        geom_dev = [
            mp.Block(center=mp.Vector3(0, 0, -sz / 4.0), size=mp.Vector3(mp.inf, mp.inf, sz / 2.0), material=clad_mat),
            mp.Block(center=mp.Vector3(0, 0, h_wg / 2.0), size=mp.Vector3(mp.inf, w_wg, h_wg), material=core_mat),
        ]
        poly_in = make_ccw_linear_taper(x_pcm_start, x_pcm_body_start, 0, 0.02, w_wg)
        geom_dev.append(create_prism_exact(poly_in, t_pcm, h_wg, pcm_mat))
        geom_dev.append(mp.Block(center=mp.Vector3((x_pcm_body_start + x_pcm_body_end) / 2.0, 0, h_wg + t_pcm / 2.0),
                                 size=mp.Vector3(l_body, w_wg, t_pcm), material=pcm_mat))
        poly_out = make_ccw_linear_taper(x_pcm_body_end, x_pcm_end, 0, w_wg, 0.02)
        geom_dev.append(create_prism_exact(poly_out, t_pcm, h_wg, pcm_mat))

        sim_dev = mp.Simulation(cell_size=cell, boundary_layers=pml_layers, geometry=geom_dev,
                                sources=sources, resolution=resolution, default_material=clad_mat)
        f_dev = sim_dev.add_flux(fcen, 0, 1, mp.FluxRegion(center=mp.Vector3(x_mon, 0, h_wg / 2.0), size=mp.Vector3(0, 0.4, 0.8)))
        sim_dev.run(until_after_sources=mp.stop_when_fields_decayed(20, mp.Hz, mp.Vector3(x_mon, 0, h_wg / 2.0), 1e-4))
        dev_mode = max(abs(sim_dev.get_eigenmode_coefficients(f_dev, [1], eig_parity=mp.NO_PARITY).alpha[0, 0, 0])**2, 1e-15)

        trans = dev_mode / ref_mode
        results[state] = trans
        print(f"t_pcm = {t_pcm_nm:5.1f}nm [{state.upper():<11}] Trans = {trans*100:10.6f}% (Loss = {(1-trans)*100:8.6f}%)")
        sys.stdout.flush()

    return results

if __name__ == "__main__":
    t0 = time.time()
    print("=" * 80)
    print("MEASURING ISOLATED PATCH TRANSMISSION VS PCM THICKNESS (l_tip = 800nm)")
    print("=" * 80)
    for t_nm in [120.0, 140.0, 160.0]:
        test_isolated_thickness(t_pcm_nm=t_nm, l_tip=0.800)
    print("=" * 80)
    print(f"Total elapsed time = {time.time() - t0:.2f} s")
