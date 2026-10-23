"""
VERIFY COINCIDENCE AT GAP = 98.0 NM IN 3D FDTD (Resolution = 20)
================================================================
"""

import sys, math, time
import numpy as np
import meep as mp

def test_gap(gap_nm=98.0000000000, resolution=20):
    gap = gap_nm * 1e-3
    lambda_0 = 1.0640000000
    fcen = 1.0 / lambda_0
    df = 0.10 * fcen

    w_wg = 0.2400000000
    h_wg = 0.2200000000
    t_pcm = 0.1600000000
    l_coupler = 6.8000000000

    y1 = 0.0000000000
    y2 = -(w_wg + gap)
    y_mid = (y1 + y2) / 2.0

    clad_mat = mp.Medium(index=1.4490000000)
    core_mat = mp.Medium(index=3.4800000000)
    pcm_am = mp.Medium(index=2.7000000000)
    pcm_cr = mp.Medium(index=3.3000000000)

    l_lead = 0.4000000000
    total_len = l_lead + l_coupler + l_lead
    dpml = 0.5000000000
    pad_y = 0.5000000000
    pad_z = 0.5000000000

    sx = total_len + 2 * dpml
    sy = abs(y2) + w_wg + 2 * (dpml + pad_y)
    sz = h_wg + t_pcm + 2 * (dpml + pad_z)

    cell = mp.Vector3(sx, sy, sz)
    pml_layers = [mp.PML(dpml)]

    x_left = -sx / 2.0
    x_c_start = x_left + dpml + l_lead
    x_c_end = x_c_start + l_coupler
    src_pos_x = x_left + dpml + 0.1500000000

    sources = [
        mp.EigenModeSource(
            src=mp.GaussianSource(fcen, fwidth=df),
            center=mp.Vector3(src_pos_x, y1, h_wg / 2.0),
            size=mp.Vector3(0, w_wg * 3.0, (h_wg + t_pcm) * 3.0),
            eig_band=1,
            direction=mp.X
        )
    ]

    x_points = np.arange(5.0, 6.21, 0.08)
    mon_wg1_size = mp.Vector3(0, abs(y_mid - y1) * 2.0, (h_wg + t_pcm) * 3.0)
    mon_wg2_size = mp.Vector3(0, abs(y_mid - y2) * 2.0, (h_wg + t_pcm) * 3.0)

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

    data = {}
    for state, mat1 in [("amorphous", pcm_am), ("crystalline", pcm_cr)]:
        geom_dev = [
            mp.Block(center=mp.Vector3(0, 0, -sz / 4.0), size=mp.Vector3(mp.inf, mp.inf, sz / 2.0), material=clad_mat),
            mp.Block(center=mp.Vector3(0, y1, h_wg / 2.0), size=mp.Vector3(mp.inf, w_wg, h_wg), material=core_mat),
            mp.Block(center=mp.Vector3((x_c_start + x_c_end) / 2.0, y1, h_wg + t_pcm / 2.0),
                     size=mp.Vector3(l_coupler, w_wg, t_pcm), material=mat1),
            mp.Block(center=mp.Vector3(0, y2, h_wg / 2.0), size=mp.Vector3(mp.inf, w_wg, h_wg), material=core_mat),
            mp.Block(center=mp.Vector3((x_c_start + x_c_end) / 2.0, y2, h_wg + t_pcm / 2.0),
                     size=mp.Vector3(l_coupler, w_wg, t_pcm), material=pcm_am),
        ]
        sim = mp.Simulation(cell_size=cell, boundary_layers=pml_layers, geometry=geom_dev,
                            sources=sources, resolution=resolution, default_material=clad_mat)
        f_list = []
        for xp in x_points:
            f1 = sim.add_flux(fcen, 0, 1, mp.FluxRegion(center=mp.Vector3(x_c_start + xp, y1, (h_wg + t_pcm)/2.0), size=mon_wg1_size))
            f2 = sim.add_flux(fcen, 0, 1, mp.FluxRegion(center=mp.Vector3(x_c_start + xp, y2, (h_wg + t_pcm)/2.0), size=mon_wg2_size))
            f_list.append((xp, f1, f2))
        sim.run(until_after_sources=mp.stop_when_fields_decayed(20, mp.Hz, mp.Vector3(x_c_end, y_mid, h_wg/2.0), 1e-4))

        results = []
        for xp, f1, f2 in f_list:
            p1 = mp.get_fluxes(f1)[0] / p_ref
            p2 = mp.get_fluxes(f2)[0] / p_ref
            frac = p2 / (p1 + p2) if state == "amorphous" else p1 / (p1 + p2)
            results.append((xp, p1*100, p2*100, (p1+p2)*100, frac*100))
        data[state] = results

    print(f"\n{'='*110}")
    print(f"[GAP = {gap_nm:.2f} nm] 3D FDTD COINCIDENCE SCAN (Resolution = {resolution})")
    print(f"{'Length x (um)':<14} {'S0 Cross Frac %':<18} {'S0 P2 Flux %':<16} {'S1 Bar Frac %':<18} {'S1 P1 Flux %':<16} {'S1 Cross Frac %':<16}")
    print("-" * 110)
    for (xp, p1_0, p2_0, tot_0, frac_0), (_, p1_1, p2_1, tot_1, frac_1) in zip(data["amorphous"], data["crystalline"]):
        s1_cross_frac = (100.0 - frac_1)
        print(f"{xp:<14.3f} {frac_0:<18.6f} {p2_0:<16.6f} {frac_1:<18.6f} {p1_1:<16.6f} {s1_cross_frac:<16.6f}")
    print(f"{'='*110}\n")
    sys.stdout.flush()

if __name__ == "__main__":
    t0 = time.time()
    test_gap(gap_nm=98.0000000000)
    print(f"Elapsed time = {time.time() - t0:.2f} s")
