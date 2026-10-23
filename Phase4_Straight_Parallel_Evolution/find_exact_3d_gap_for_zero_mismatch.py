"""
FIND EXACT 3D FDTD EVANESCENT GAP FOR ZERO-MISMATCH (Delta L = 0)
==================================================================
Runs 3D FDTD for Gap in [85, 95, 105] nm to find exact Delta L = 0 crossing.
For each gap:
  - Tracks State 0 power transfer vs x to find exact peak L_c0
  - Tracks State 1 power transfer vs x to find exact null L_c1
  - Computes Delta L(gap) = L_c0 - L_c1
"""

import sys, math, time
import numpy as np
import meep as mp

def measure_lc_for_gap(gap_nm=95.0, resolution=15):
    gap = gap_nm * 1e-3
    lambda_0 = 1.0640000000
    fcen = 1.0 / lambda_0
    df = 0.10 * fcen

    w_wg = 0.2400000000
    h_wg = 0.2200000000
    t_pcm = 0.1600000000
    l_coupler = 7.0000000000

    y1 = 0.0000000000
    y2 = -(w_wg + gap)
    y_mid = (y1 + y2) / 2.0

    clad_mat = mp.Medium(index=1.4490000000)
    core_mat = mp.Medium(index=3.4800000000)
    pcm_am = mp.Medium(index=2.7000000000)
    pcm_cr = mp.Medium(index=3.3000000000)

    l_lead = 0.3500000000
    total_len = l_lead + l_coupler + l_lead
    dpml = 0.5000000000
    pad = 0.4000000000

    sx = total_len + 2 * dpml
    sy = abs(y2) + w_wg + 2 * (dpml + pad)
    sz = h_wg + t_pcm + 2 * (dpml + pad)

    cell = mp.Vector3(sx, sy, sz)
    pml_layers = [mp.PML(dpml)]

    x_left = -sx / 2.0
    x_c_start = x_left + dpml + l_lead
    x_c_end = x_c_start + l_coupler
    src_pos_x = x_left + dpml + 0.1200000000

    sources = [
        mp.EigenModeSource(
            src=mp.GaussianSource(fcen, fwidth=df),
            center=mp.Vector3(src_pos_x, y1, h_wg / 2.0),
            size=mp.Vector3(0, w_wg * 3.0, (h_wg + t_pcm) * 3.0),
            eig_band=1,
            direction=mp.X
        )
    ]

    # Reference
    geom_ref = [
        mp.Block(center=mp.Vector3(0, 0, -sz / 4.0), size=mp.Vector3(mp.inf, mp.inf, sz / 2.0), material=clad_mat),
        mp.Block(center=mp.Vector3(0, y1, h_wg / 2.0), size=mp.Vector3(mp.inf, w_wg, h_wg), material=core_mat),
        mp.Block(center=mp.Vector3(0, y1, h_wg + t_pcm / 2.0), size=mp.Vector3(mp.inf, w_wg, t_pcm), material=pcm_am)
    ]
    sim_ref = mp.Simulation(cell_size=cell, boundary_layers=pml_layers, geometry=geom_ref,
                            sources=sources, resolution=resolution, default_material=clad_mat)
    f_ref = sim_ref.add_flux(fcen, 0, 1, mp.FluxRegion(center=mp.Vector3(x_c_start, y1, h_wg / 2.0), size=mp.Vector3(0, 1.2, 1.2)))
    sim_ref.run(until_after_sources=mp.stop_when_fields_decayed(20, mp.Hz, mp.Vector3(x_c_start, y1, h_wg / 2.0), 1e-4))
    p_ref = max(mp.get_fluxes(f_ref)[0], 1e-15)

    x_eval = np.linspace(3.5, 6.8, 34)
    mon_wg1_size = mp.Vector3(0, abs(y_mid - y1) * 2.0, (h_wg + t_pcm) * 3.0)
    mon_wg2_size = mp.Vector3(0, abs(y_mid - y2) * 2.0, (h_wg + t_pcm) * 3.0)

    # State 0 (amorphous)
    geom_s0 = [
        mp.Block(center=mp.Vector3(0, 0, -sz / 4.0), size=mp.Vector3(mp.inf, mp.inf, sz / 2.0), material=clad_mat),
        mp.Block(center=mp.Vector3(0, y1, h_wg / 2.0), size=mp.Vector3(mp.inf, w_wg, h_wg), material=core_mat),
        mp.Block(center=mp.Vector3((x_c_start + x_c_end) / 2.0, y1, h_wg + t_pcm / 2.0),
                 size=mp.Vector3(l_coupler, w_wg, t_pcm), material=pcm_am),
        mp.Block(center=mp.Vector3(0, y2, h_wg / 2.0), size=mp.Vector3(mp.inf, w_wg, h_wg), material=core_mat),
        mp.Block(center=mp.Vector3((x_c_start + x_c_end) / 2.0, y2, h_wg + t_pcm / 2.0),
                 size=mp.Vector3(l_coupler, w_wg, t_pcm), material=pcm_am),
    ]
    sim_s0 = mp.Simulation(cell_size=cell, boundary_layers=pml_layers, geometry=geom_s0,
                           sources=sources, resolution=resolution, default_material=clad_mat)
    flist_s0 = []
    for xp in x_eval:
        f1 = sim_s0.add_flux(fcen, 0, 1, mp.FluxRegion(center=mp.Vector3(x_c_start + xp, y1, (h_wg + t_pcm)/2.0), size=mon_wg1_size))
        f2 = sim_s0.add_flux(fcen, 0, 1, mp.FluxRegion(center=mp.Vector3(x_c_start + xp, y2, (h_wg + t_pcm)/2.0), size=mon_wg2_size))
        flist_s0.append((xp, f1, f2))
    sim_s0.run(until_after_sources=mp.stop_when_fields_decayed(20, mp.Hz, mp.Vector3(x_c_end, y_mid, h_wg/2.0), 1e-4))

    p2_frac_s0 = []
    for xp, f1, f2 in flist_s0:
        p1 = mp.get_fluxes(f1)[0]
        p2 = mp.get_fluxes(f2)[0]
        p2_frac_s0.append(p2 / (p1 + p2))

    # Fit quadratic near maximum to find exact peak L_c0
    idx_max0 = int(np.argmax(p2_frac_s0))
    idx_fit0 = slice(max(0, idx_max0 - 3), min(len(x_eval), idx_max0 + 4))
    poly0 = np.polyfit(x_eval[idx_fit0], p2_frac_s0[idx_fit0], 2)
    lc0 = -poly0[1] / (2.0 * poly0[0])

    # State 1 (crystalline)
    geom_s1 = [
        mp.Block(center=mp.Vector3(0, 0, -sz / 4.0), size=mp.Vector3(mp.inf, mp.inf, sz / 2.0), material=clad_mat),
        mp.Block(center=mp.Vector3(0, y1, h_wg / 2.0), size=mp.Vector3(mp.inf, w_wg, h_wg), material=core_mat),
        mp.Block(center=mp.Vector3((x_c_start + x_c_end) / 2.0, y1, h_wg + t_pcm / 2.0),
                 size=mp.Vector3(l_coupler, w_wg, t_pcm), material=pcm_cr),
        mp.Block(center=mp.Vector3(0, y2, h_wg / 2.0), size=mp.Vector3(mp.inf, w_wg, h_wg), material=core_mat),
        mp.Block(center=mp.Vector3((x_c_start + x_c_end) / 2.0, y2, h_wg + t_pcm / 2.0),
                 size=mp.Vector3(l_coupler, w_wg, t_pcm), material=pcm_am),
    ]
    sim_s1 = mp.Simulation(cell_size=cell, boundary_layers=pml_layers, geometry=geom_s1,
                           sources=sources, resolution=resolution, default_material=clad_mat)
    flist_s1 = []
    for xp in x_eval:
        f1 = sim_s1.add_flux(fcen, 0, 1, mp.FluxRegion(center=mp.Vector3(x_c_start + xp, y1, (h_wg + t_pcm)/2.0), size=mon_wg1_size))
        f2 = sim_s1.add_flux(fcen, 0, 1, mp.FluxRegion(center=mp.Vector3(x_c_start + xp, y2, (h_wg + t_pcm)/2.0), size=mon_wg2_size))
        flist_s1.append((xp, f1, f2))
    sim_s1.run(until_after_sources=mp.stop_when_fields_decayed(20, mp.Hz, mp.Vector3(x_c_end, y_mid, h_wg/2.0), 1e-4))

    p1_frac_s1 = []
    for xp, f1, f2 in flist_s1:
        p1 = mp.get_fluxes(f1)[0]
        p2 = mp.get_fluxes(f2)[0]
        p1_frac_s1.append(p1 / (p1 + p2))

    # Find peak of p1_frac_s1 (which is the null of cross-talk)
    # Search in range x > 4.5
    valid_mask = x_eval > 4.5
    valid_x = x_eval[valid_mask]
    valid_p1 = np.array(p1_frac_s1)[valid_mask]
    idx_max1 = int(np.argmax(valid_p1))
    fit_x = valid_x[max(0, idx_max1-3):min(len(valid_x), idx_max1+4)]
    fit_y = valid_p1[max(0, idx_max1-3):min(len(valid_p1), idx_max1+4)]
    poly1 = np.polyfit(fit_x, fit_y, 2)
    lc1 = -poly1[1] / (2.0 * poly1[0])

    delta_l = (lc0 - lc1) * 1e3 # in nm
    print(f"Gap = {gap_nm:6.2f} nm -> L_c0 = {lc0:8.5f} um, L_c1 = {lc1:8.5f} um, Delta_L = {delta_l:+8.2f} nm")
    sys.stdout.flush()
    return gap_nm, lc0, lc1, delta_l

if __name__ == "__main__":
    t0 = time.time()
    print("=" * 85)
    print("MAPPING 3D FDTD ZERO-MISMATCH GAP (Delta L vs Gap)")
    print("=" * 85)
    results = []
    for g in [85.0, 95.0, 105.51]:
        res = measure_lc_for_gap(gap_nm=g)
        results.append(res)
    print("=" * 85)
    print(f"{'Gap (nm)':<12} {'L_c0 (um)':<14} {'L_c1 (um)':<14} {'Delta L (nm)':<16}")
    print("-" * 85)
    for g, l0, l1, dl in results:
        print(f"{g:<12.2f} {l0:<14.5f} {l1:<14.5f} {dl:<+16.2f}")
    print("=" * 85)
    print(f"Total mapping time = {time.time() - t0:.2f} s")
