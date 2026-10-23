"""
PROJECT JANUS: 3D FULL-WAVE MAXWELL MEEP FDTD SIMULATION
2x2 MMI 3dB POWER SPLITTER AT 1064 nm
=====================================
Target: 50:50 Power Split, Total Transmission > 97%
Dimensions:
  - Multimode width W_mmi = 1.20 um
  - Multimode length L_mmi = 2.30 um (sweep [2.0, 2.6] um)
  - Port separation: y_sep = 0.55 um (y = +0.275 um, -0.275 um)
  - Port width w_wg = 0.250 um, h_wg = 0.220 um
"""
import sys
import os
import math
import time
import numpy as np
import meep as mp

def run_mmi_test(l_mmi=2.30, w_mmi=1.20, y_port=0.275, resolution=20):
    lambda_0 = 1.064
    fcen = 1.0 / lambda_0
    df = 0.1 * fcen

    w_wg = 0.250
    h_wg = 0.220
    l_in = 0.60
    l_out = 0.60

    clad_mat = mp.Medium(index=1.449)
    core_mat = mp.Medium(index=3.480)

    cell_len = l_in + l_mmi + l_out
    cell_w = w_mmi + 0.60

    dpml = 0.50
    pad_x = 0.30
    pad_y = 0.30
    pad_z = 0.35
    sx = cell_len + 2 * (dpml + pad_x)
    sy = cell_w + 2 * (dpml + pad_y)
    sz = h_wg + 2 * (dpml + pad_z)

    cell = mp.Vector3(sx, sy, sz)
    pml_layers = [mp.PML(dpml)]

    x_left = -sx / 2.0
    x_right = sx / 2.0
    x_mmi_start = -l_mmi / 2.0
    x_mmi_end = l_mmi / 2.0

    sources = [
        mp.EigenModeSource(
            src=mp.GaussianSource(fcen, fwidth=df),
            center=mp.Vector3(x_left + dpml + 0.10, y_port, h_wg / 2.0),
            size=mp.Vector3(0, w_wg * 3.0, h_wg * 3.0),
            eig_band=1,
            direction=mp.X
        )
    ]

    # Reference straight waveguide
    x_mon = x_mmi_end + 0.30
    mon_w = w_wg * 2.5
    mon_h = h_wg * 2.5

    sim_ref = mp.Simulation(
        cell_size=cell, boundary_layers=pml_layers,
        geometry=[
            mp.Block(center=mp.Vector3(0, 0, -sz/4.0), size=mp.Vector3(mp.inf, mp.inf, sz/2.0), material=clad_mat),
            mp.Block(center=mp.Vector3(0, y_port, h_wg/2.0), size=mp.Vector3(mp.inf, w_wg, h_wg), material=core_mat)
        ],
        sources=sources, resolution=resolution, default_material=clad_mat
    )
    flux_ref = sim_ref.add_flux(fcen, 0, 1, mp.FluxRegion(
        center=mp.Vector3(x_mon, y_port, h_wg/2.0), size=mp.Vector3(0, mon_w, mon_h)
    ))
    sim_ref.run(until_after_sources=mp.stop_when_fields_decayed(20, mp.Hz, mp.Vector3(x_mon, y_port, h_wg/2.0), 1e-3))
    p_ref = max(mp.get_fluxes(flux_ref)[0], 1e-12)

    # 2x2 MMI Geometry
    geom = [
        mp.Block(center=mp.Vector3(0, 0, -sz/4.0), size=mp.Vector3(mp.inf, mp.inf, sz/2.0), material=clad_mat),
        # Input access waveguides
        mp.Block(center=mp.Vector3((x_left + x_mmi_start)/2.0, y_port, h_wg/2.0), size=mp.Vector3(x_mmi_start - x_left, w_wg, h_wg), material=core_mat),
        mp.Block(center=mp.Vector3((x_left + x_mmi_start)/2.0, -y_port, h_wg/2.0), size=mp.Vector3(x_mmi_start - x_left, w_wg, h_wg), material=core_mat),
        # MMI central multimode section
        mp.Block(center=mp.Vector3(0, 0, h_wg/2.0), size=mp.Vector3(l_mmi, w_mmi, h_wg), material=core_mat),
        # Output access waveguides
        mp.Block(center=mp.Vector3((x_mmi_end + x_right)/2.0, y_port, h_wg/2.0), size=mp.Vector3(x_right - x_mmi_end, w_wg, h_wg), material=core_mat),
        mp.Block(center=mp.Vector3((x_mmi_end + x_right)/2.0, -y_port, h_wg/2.0), size=mp.Vector3(x_right - x_mmi_end, w_wg, h_wg), material=core_mat),
    ]

    sim = mp.Simulation(
        cell_size=cell, boundary_layers=pml_layers, geometry=geom,
        sources=sources, resolution=resolution, default_material=clad_mat
    )
    flux_p1 = sim.add_flux(fcen, 0, 1, mp.FluxRegion(center=mp.Vector3(x_mon, y_port, h_wg/2.0), size=mp.Vector3(0, mon_w, mon_h)))
    flux_p2 = sim.add_flux(fcen, 0, 1, mp.FluxRegion(center=mp.Vector3(x_mon, -y_port, h_wg/2.0), size=mp.Vector3(0, mon_w, mon_h)))

    sim.run(until_after_sources=mp.stop_when_fields_decayed(20, mp.Hz, mp.Vector3(x_mon, y_port, h_wg/2.0), 1e-3))

    p1 = max(mp.get_fluxes(flux_p1)[0] / p_ref, 1e-12)
    p2 = max(mp.get_fluxes(flux_p2)[0] / p_ref, 1e-12)
    tot = p1 + p2

    return {
        "l_mmi": l_mmi, "w_mmi": w_mmi, "p1": p1, "p2": p2, "tot": tot,
        "p1_pct": p1 * 100, "p2_pct": p2 * 100, "tot_pct": tot * 100,
        "loss_db": -10.0 * math.log10(tot), "imbalance_db": abs(10.0 * math.log10(p1 / p2))
    }

if __name__ == "__main__":
    l = float(sys.argv[1]) if len(sys.argv) > 1 else 2.30
    res = run_mmi_test(l_mmi=l)
    print("\n" + "=" * 60)
    print(f"2x2 MMI RESULTS (L_mmi = {res['l_mmi']:.2f} um, W_mmi = {res['w_mmi']:.2f} um):")
    print(f"  Port 1 (Bar)   : {res['p1_pct']:.2f}%")
    print(f"  Port 2 (Cross) : {res['p2_pct']:.2f}%")
    print(f"  Total Trans    : {res['tot_pct']:.2f}%")
    print(f"  Excess Loss    : {res['loss_db']:.3f} dB")
    print(f"  Power Imbalance: {res['imbalance_db']:.2f} dB")
    print("=" * 60)
