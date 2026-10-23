"""
PROJECT JANUS: 3D FULL-WAVE VECTOR MAXWELL MEEP FDTD
PROTOTYPE: SINGLE-CAVITY 1x2 MMI SWITCH WITH ASYMMETRIC Sb2S3 PATCH
===================================================================
Tests geometry creation, eigenmode launch, and flux integration.
"""
import sys
import os
import math
import numpy as np
import meep as mp

def run_single_cavity_test(
    w_mmi=1.00,
    l_mmi=3.00,
    w_wg=0.250,
    h_wg=0.220,
    y_port=0.275,
    t_pcm=0.030,
    l_pcm=2.00,
    state="amorphous",
    resolution=20
):
    lambda_0 = 1.064
    fcen = 1.0 / lambda_0
    df = 0.1 * fcen

    l_in = 0.400
    l_out = 0.400

    clad_mat = mp.Medium(index=1.449)
    core_mat = mp.Medium(index=3.480)

    # Sb2S3 material properties
    k_am = 1.0e-5
    k_cr = 8.0e-5
    pcm_am = mp.Medium(index=2.700, D_conductivity=2 * math.pi * fcen * k_am)
    pcm_cr = mp.Medium(index=3.300, D_conductivity=2 * math.pi * fcen * k_cr)
    pcm_mat = pcm_cr if state == "crystalline" else pcm_am

    cell_len = l_in + l_mmi + l_out
    cell_w = max(w_mmi, 2 * y_port + w_wg) + 0.40
    cell_area = cell_len * cell_w

    dpml = 0.40
    pad_x = 0.20
    pad_y = 0.20
    pad_z = 0.30
    sx = cell_len + 2 * (dpml + pad_x)
    sy = cell_w + 2 * (dpml + pad_y)
    sz = h_wg + t_pcm + 2 * (dpml + pad_z)

    cell = mp.Vector3(sx, sy, sz)
    pml_layers = [mp.PML(dpml)]

    x_left = -sx / 2.0
    x_right = sx / 2.0
    x_mmi_start = -l_mmi / 2.0
    x_mmi_end = l_mmi / 2.0

    # Source at input
    x_src = x_left + dpml + 0.08
    sources = [
        mp.EigenModeSource(
            src=mp.GaussianSource(fcen, fwidth=df),
            center=mp.Vector3(x_src, y_port, h_wg / 2.0),
            size=mp.Vector3(0, w_wg * 3.0, h_wg * 3.0),
            eig_band=1,
            direction=mp.X
        )
    ]

    # Reference straight waveguide
    mon_w = 0.260
    mon_h = 0.280
    x_mon = x_mmi_end + 0.15

    geom_ref = [
        mp.Block(center=mp.Vector3(0, 0, -sz/4.0), size=mp.Vector3(mp.inf, mp.inf, sz/2.0), material=clad_mat),
        mp.Block(center=mp.Vector3(0, y_port, h_wg/2.0), size=mp.Vector3(mp.inf, w_wg, h_wg), material=core_mat)
    ]

    sim_ref = mp.Simulation(
        cell_size=cell, boundary_layers=pml_layers, geometry=geom_ref,
        sources=sources, resolution=resolution, default_material=clad_mat
    )
    flux_ref = sim_ref.add_flux(fcen, 0, 1, mp.FluxRegion(center=mp.Vector3(x_mon, y_port, h_wg/2.0), size=mp.Vector3(0, mon_w, mon_h)))
    sim_ref.run(until_after_sources=mp.stop_when_fields_decayed(15, mp.Hz, mp.Vector3(x_mon, y_port, h_wg/2.0), 1e-3))
    p_ref = max(mp.get_fluxes(flux_ref)[0], 1e-12)

    # 1x2 Single-Cavity MMI Geometry
    geom = [
        mp.Block(center=mp.Vector3(0, 0, -sz/4.0), size=mp.Vector3(mp.inf, mp.inf, sz/2.0), material=clad_mat),
        # 1 Single Input access lead (at +y_port)
        mp.Block(center=mp.Vector3((x_left + x_mmi_start)/2.0, y_port, h_wg/2.0), size=mp.Vector3(x_mmi_start - x_left, w_wg, h_wg), material=core_mat),
        # MMI Multimode central body
        mp.Block(center=mp.Vector3(0, 0, h_wg/2.0), size=mp.Vector3(l_mmi, w_mmi, h_wg), material=core_mat),
        # Asymmetric Sb2S3 Patch strictly covering upper half (y in [0, w_mmi/2])
        mp.Block(
            center=mp.Vector3(0, w_mmi / 4.0, h_wg + t_pcm / 2.0),
            size=mp.Vector3(l_pcm, w_mmi / 2.0, t_pcm),
            material=pcm_mat
        ),
        # Output access leads (Port 1 at +y_port, Port 2 at -y_port)
        mp.Block(center=mp.Vector3((x_mmi_end + x_right)/2.0, y_port, h_wg/2.0), size=mp.Vector3(x_right - x_mmi_end, w_wg, h_wg), material=core_mat),
        mp.Block(center=mp.Vector3((x_mmi_end + x_right)/2.0, -y_port, h_wg/2.0), size=mp.Vector3(x_right - x_mmi_end, w_wg, h_wg), material=core_mat),
    ]

    sim = mp.Simulation(
        cell_size=cell, boundary_layers=pml_layers, geometry=geom,
        sources=sources, resolution=resolution, default_material=clad_mat
    )
    flux_p1 = sim.add_flux(fcen, 0, 1, mp.FluxRegion(center=mp.Vector3(x_mon, y_port, h_wg/2.0), size=mp.Vector3(0, mon_w, mon_h)))
    flux_p2 = sim.add_flux(fcen, 0, 1, mp.FluxRegion(center=mp.Vector3(x_mon, -y_port, h_wg/2.0), size=mp.Vector3(0, mon_w, mon_h)))

    sim.run(until_after_sources=mp.stop_when_fields_decayed(15, mp.Hz, mp.Vector3(x_mon, y_port, h_wg/2.0), 1e-3))

    p1 = max(mp.get_fluxes(flux_p1)[0] / p_ref, 1e-12)
    p2 = max(mp.get_fluxes(flux_p2)[0] / p_ref, 1e-12)
    tot = p1 + p2

    return {
        "state": state, "w_mmi": w_mmi, "l_mmi": l_mmi, "l_pcm": l_pcm, "t_pcm": t_pcm,
        "area": cell_area, "p1": p1, "p2": p2, "tot": tot,
        "p1_pct": p1 * 100, "p2_pct": p2 * 100, "tot_pct": tot * 100,
        "loss_db": -10.0 * math.log10(tot),
        "xt_p1_db": 10.0 * math.log10(p2 / max(p1, 1e-12)),
        "xt_p2_db": 10.0 * math.log10(p1 / max(p2, 1e-12))
    }

if __name__ == "__main__":
    state = sys.argv[1] if len(sys.argv) > 1 else "amorphous"
    res = run_single_cavity_test(state=state)
    print("\n" + "=" * 65)
    print(f"SINGLE-CAVITY 1x2 MMI 3D RESULT ({res['state'].upper()}):")
    print(f"  Footprint Area   : {res['area']:.3f} um^2")
    print(f"  Port 1 (Bar)     : {res['p1_pct']:.2f}%")
    print(f"  Port 2 (Cross)   : {res['p2_pct']:.2f}%")
    print(f"  Total Trans      : {res['tot_pct']:.2f}%")
    print(f"  Excess Loss      : {res['loss_db']:.4f} dB")
    print(f"  Crosstalk (P1)   : {res['xt_p1_db']:.2f} dB")
    print(f"  Crosstalk (P2)   : {res['xt_p2_db']:.2f} dB")
    print("=" * 65)
