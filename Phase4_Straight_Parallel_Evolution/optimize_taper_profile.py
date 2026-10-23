"""
TAPER PROFILE OPTIMIZATION FOR ZERO RADIATION LOSS IN CRYSTALLINE PCM
======================================================================
Tests:
- Linear vs Quadratic vs Cubic Hermite taper profiles
- Tip lengths: 800 nm, 1000 nm
- Caches reference waveguide to minimize simulation time
"""

import sys, math, time
import numpy as np
import meep as mp

def make_ccw_taper(x_start, x_end, y, w_start, w_end, profile="hermite", n_pts=30):
    dx = x_end - x_start
    u = np.linspace(0, 1, n_pts)
    if profile == "linear":
        w = w_start + u * (w_end - w_start)
    elif profile == "quadratic":
        w = w_start + (u**2) * (w_end - w_start)
    elif profile == "hermite":
        w = w_start + (3 * u**2 - 2 * u**3) * (w_end - w_start)
    elif profile == "sqrt":
        w = w_start + np.sqrt(u) * (w_end - w_start)
    else:
        w = w_start + u * (w_end - w_start)

    x = x_start + u * dx
    bot = [mp.Vector3(xi, y - wi / 2.0, 0) for xi, wi in zip(x, w)]
    top = [mp.Vector3(xi, y + wi / 2.0, 0) for xi, wi in reversed(list(zip(x, w)))]
    return bot + top

def create_prism_exact(vertices, height, z_base, material):
    cx = float(np.mean([pt.x for pt in vertices]))
    cy = float(np.mean([pt.y for pt in vertices]))
    cz = float(z_base + height / 2.0)
    return mp.Prism(vertices=vertices, height=height, axis=mp.Vector3(0, 0, 1),
                    center=mp.Vector3(cx, cy, cz), material=material)

ref_cache = {}

def get_reference(cell_len, l_tip, l_body, t_pcm, resolution=15):
    if cell_len in ref_cache:
        return ref_cache[cell_len]

    lambda_0 = 1.064
    fcen = 1.0 / lambda_0
    df = 0.10 * fcen

    w_wg = 0.240
    h_wg = 0.220

    clad_mat = mp.Medium(index=1.449)
    core_mat = mp.Medium(index=3.480)

    dpml = 0.50
    pad = 0.35
    sx = cell_len + 2 * (dpml + pad)
    sy = 1.5 + 2 * (dpml + pad)
    sz = h_wg + t_pcm + 2 * (dpml + pad)

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
            size=mp.Vector3(0, w_wg * 3.0, (h_wg + t_pcm) * 3.0),
            eig_band=1,
            direction=mp.X
        )
    ]

    geom_ref = [
        mp.Block(center=mp.Vector3(0, 0, -sz / 4.0), size=mp.Vector3(mp.inf, mp.inf, sz / 2.0), material=clad_mat),
        mp.Block(center=mp.Vector3(0, 0, h_wg / 2.0), size=mp.Vector3(mp.inf, w_wg, h_wg), material=core_mat)
    ]
    sim_ref = mp.Simulation(cell_size=cell, boundary_layers=pml_layers, geometry=geom_ref,
                            sources=sources, resolution=resolution, default_material=clad_mat)
    f_ref = sim_ref.add_flux(fcen, 0, 1, mp.FluxRegion(center=mp.Vector3(x_mon, 0, h_wg / 2.0), size=mp.Vector3(0, 0.4, 0.8)))
    sim_ref.run(until_after_sources=mp.stop_when_fields_decayed(20, mp.Hz, mp.Vector3(x_mon, 0, h_wg / 2.0), 1e-4))
    ref_mode = max(abs(sim_ref.get_eigenmode_coefficients(f_ref, [1], eig_parity=mp.NO_PARITY).alpha[0, 0, 0])**2, 1e-15)
    ref_cache[cell_len] = ref_mode
    return ref_mode

def test_taper_profile(profile="hermite", l_tip=0.800, t_pcm=0.2375, state="crystalline", resolution=15):
    lambda_0 = 1.064
    fcen = 1.0 / lambda_0
    df = 0.10 * fcen

    w_wg = 0.240
    h_wg = 0.220
    l_body = 4.000

    clad_mat = mp.Medium(index=1.449)
    core_mat = mp.Medium(index=3.480)
    pcm_mat = mp.Medium(index=3.30 if state == "crystalline" else 2.70)

    cell_len = l_tip + l_body + l_tip + 1.0
    cell_w = 1.5
    dpml = 0.50
    pad = 0.35
    sx = cell_len + 2 * (dpml + pad)
    sy = cell_w + 2 * (dpml + pad)
    sz = h_wg + t_pcm + 2 * (dpml + pad)

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

    ref_mode = get_reference(cell_len, l_tip, l_body, t_pcm, resolution=resolution)

    src_pos_x = x_left + dpml + 0.10
    sources = [
        mp.EigenModeSource(
            src=mp.GaussianSource(fcen, fwidth=df),
            center=mp.Vector3(src_pos_x, 0, h_wg / 2.0),
            size=mp.Vector3(0, w_wg * 3.0, (h_wg + t_pcm) * 3.0),
            eig_band=1,
            direction=mp.X
        )
    ]

    geom_dev = [
        mp.Block(center=mp.Vector3(0, 0, -sz / 4.0), size=mp.Vector3(mp.inf, mp.inf, sz / 2.0), material=clad_mat),
        mp.Block(center=mp.Vector3(0, 0, h_wg / 2.0), size=mp.Vector3(mp.inf, w_wg, h_wg), material=core_mat),
    ]
    poly_in = make_ccw_taper(x_pcm_start, x_pcm_body_start, 0, 0.02, w_wg, profile=profile)
    geom_dev.append(create_prism_exact(poly_in, t_pcm, h_wg, pcm_mat))

    geom_dev.append(mp.Block(center=mp.Vector3((x_pcm_body_start + x_pcm_body_end) / 2.0, 0, h_wg + t_pcm / 2.0),
                             size=mp.Vector3(l_body, w_wg, t_pcm), material=pcm_mat))

    poly_out = make_ccw_taper(x_pcm_body_end, x_pcm_end, 0, w_wg, 0.02, profile=profile)
    geom_dev.append(create_prism_exact(poly_out, t_pcm, h_wg, pcm_mat))

    sim_dev = mp.Simulation(cell_size=cell, boundary_layers=pml_layers, geometry=geom_dev,
                            sources=sources, resolution=resolution, default_material=clad_mat)
    f_dev = sim_dev.add_flux(fcen, 0, 1, mp.FluxRegion(center=mp.Vector3(x_mon, 0, h_wg / 2.0), size=mp.Vector3(0, 0.4, 0.8)))
    sim_dev.run(until_after_sources=mp.stop_when_fields_decayed(20, mp.Hz, mp.Vector3(x_mon, 0, h_wg / 2.0), 1e-4))
    dev_mode = max(abs(sim_dev.get_eigenmode_coefficients(f_dev, [1], eig_parity=mp.NO_PARITY).alpha[0, 0, 0])**2, 1e-15)

    trans = dev_mode / ref_mode
    print(f"[{state.upper()}] profile={profile:<10} l_tip={l_tip*1e3:4.0f}nm: Trans = {trans*100:10.6f}% (Loss = {(1-trans)*100:8.6f}%)")
    sys.stdout.flush()
    return trans

if __name__ == "__main__":
    t0 = time.time()
    print("=" * 80)
    print("TESTING TAPER PROFILES FOR LOW-LOSS CRYSTALLINE PCM ADIABATIC TRANSITION")
    print("=" * 80)
    for l_tip in [0.800, 1.000]:
        for prof in ["hermite", "quadratic", "linear"]:
            test_taper_profile(profile=prof, l_tip=l_tip, state="crystalline")
    print("=" * 80)
    print(f"Total scan elapsed time = {time.time() - t0:.2f} s")
