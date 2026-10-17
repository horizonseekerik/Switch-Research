"""
PROJECT JANUS: 3D FIELD MONITOR FOR EXACT COUPLING DYNAMICS
Tracks Poynting flux P1(x) and P2(x) along the coupler in 3D MEEP FDTD
with non-overlapping partition flux monitors.
"""
import sys
import os
import math
import numpy as np
import meep as mp

def monitor_coupler_evolution(gap=0.100, w=0.240, t_pcm=0.150, mode="symmetric_am", resolution=20):
    lambda_0 = 1.064
    fcen = 1.0 / lambda_0
    df = 0.1 * fcen
    h_wg = 0.220
    l_tot = 8.00  # 8 um long straight coupling section

    y_sep = w + gap
    y1 = y_sep / 2.0
    y2 = -y_sep / 2.0

    dpml = 0.50
    pad_x = 0.35
    pad_y = 0.35
    pad_z = 0.35
    sx = l_tot + 2 * (dpml + pad_x)
    sy = y_sep + w + 2 * (dpml + pad_y)
    sz = h_wg + t_pcm + 2 * (dpml + pad_z)

    cell = mp.Vector3(sx, sy, sz)
    pml_layers = [mp.PML(dpml)]

    x_left = -sx / 2.0
    x_right = sx / 2.0

    clad_mat = mp.Medium(index=1.449)
    core_mat = mp.Medium(index=3.480)

    pcm_am = mp.Medium(index=2.70, D_conductivity=2*math.pi*fcen*1e-5)
    pcm_cr = mp.Medium(index=3.30, D_conductivity=2*math.pi*fcen*8e-5)

    sources = [
        mp.EigenModeSource(
            src=mp.GaussianSource(fcen, fwidth=df),
            center=mp.Vector3(x_left + dpml + 0.10, y1, h_wg / 2.0),
            size=mp.Vector3(0, w * 3.0, h_wg * 3.0),
            eig_band=1,
            direction=mp.X
        )
    ]

    mon_w = y_sep       # Exact half-plane partition width
    mon_h = 0.50        # Large vertical capture
    x_ref_mon = -3.0

    # 1. Reference straight waveguide
    print(f"[{mode.upper()}] 1. Running Reference Straight Waveguide...")
    sim_ref = mp.Simulation(
        cell_size=cell, boundary_layers=pml_layers,
        geometry=[
            mp.Block(center=mp.Vector3(0, 0, -sz/4.0), size=mp.Vector3(mp.inf, mp.inf, sz/2.0), material=clad_mat),
            mp.Block(center=mp.Vector3(0, y1, h_wg/2.0), size=mp.Vector3(mp.inf, w, h_wg), material=core_mat)
        ],
        sources=sources, resolution=resolution, default_material=clad_mat
    )
    f_ref = sim_ref.add_flux(fcen, 0, 1, mp.FluxRegion(center=mp.Vector3(x_ref_mon, y1, h_wg/2.0), size=mp.Vector3(0, mon_w, mon_h)))
    sim_ref.run(until_after_sources=mp.stop_when_fields_decayed(20, mp.Hz, mp.Vector3(x_ref_mon, y1, h_wg/2.0), 1e-3))
    p_ref = max(mp.get_fluxes(f_ref)[0], 1e-12)
    print(f"[{mode.upper()}] Reference Flux = {p_ref:.6e}")

    # 2. Coupled Waveguides
    # Monitors every 0.25 um from x = -3.0 to +3.0 um
    x_positions = np.arange(-3.0, 3.25, 0.25)

    geom = [
        mp.Block(center=mp.Vector3(0, 0, -sz/4.0), size=mp.Vector3(mp.inf, mp.inf, sz/2.0), material=clad_mat),
        mp.Block(center=mp.Vector3(0, y1, h_wg/2.0), size=mp.Vector3(mp.inf, w, h_wg), material=core_mat),
        mp.Block(center=mp.Vector3(0, y2, h_wg/2.0), size=mp.Vector3(mp.inf, w, h_wg), material=core_mat),
    ]

    if mode == "symmetric_am":
        mat1, mat2 = pcm_am, pcm_am
    elif mode == "symmetric_cr":
        mat1, mat2 = pcm_cr, pcm_cr
    elif mode == "asymmetric_cr":
        mat1, mat2 = pcm_cr, pcm_am
    elif mode == "bare":
        mat1, mat2 = None, None
    else:
        raise ValueError("Unknown mode: " + mode)

    if mat1 is not None:
        geom.append(mp.Block(center=mp.Vector3(0, y1, h_wg + t_pcm/2.0), size=mp.Vector3(mp.inf, w, t_pcm), material=mat1))
    if mat2 is not None:
        geom.append(mp.Block(center=mp.Vector3(0, y2, h_wg + t_pcm/2.0), size=mp.Vector3(mp.inf, w, t_pcm), material=mat2))

    print(f"[{mode.upper()}] 2. Running Coupled Evolution Tracking...")
    sim = mp.Simulation(
        cell_size=cell, boundary_layers=pml_layers, geometry=geom,
        sources=sources, resolution=resolution, default_material=clad_mat
    )

    flux_p1 = [sim.add_flux(fcen, 0, 1, mp.FluxRegion(center=mp.Vector3(xi, y1, h_wg/2.0), size=mp.Vector3(0, mon_w, mon_h))) for xi in x_positions]
    flux_p2 = [sim.add_flux(fcen, 0, 1, mp.FluxRegion(center=mp.Vector3(xi, y2, h_wg/2.0), size=mp.Vector3(0, mon_w, mon_h))) for xi in x_positions]

    sim.run(until_after_sources=mp.stop_when_fields_decayed(20, mp.Hz, mp.Vector3(x_positions[-1], y1, h_wg/2.0), 1e-3))

    results = []
    for i, xi in enumerate(x_positions):
        p1 = mp.get_fluxes(flux_p1[i])[0] / p_ref
        p2 = mp.get_fluxes(flux_p2[i])[0] / p_ref
        dist = xi - x_positions[0]
        results.append((dist, p1, p2, p1 + p2))
    return results

if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else "symmetric_am"
    g = float(sys.argv[2]) if len(sys.argv) > 2 else 0.100
    res = monitor_coupler_evolution(gap=g, mode=mode)
    print("\n" + "=" * 70)
    print(f"POWER EVOLUTION ALONG COUPLER (mode={mode}, gap={g*1000:.0f} nm):")
    print(f"{'Dist (um)':<10} {'P1 (WG1)':<12} {'P2 (WG2)':<12} {'Total':<12} {'P2 Ratio':<12}")
    print("-" * 70)
    for dist, p1, p2, tot in res:
        ratio = (p2 / tot * 100.0) if tot > 0 else 0
        print(f"{dist:<10.2f} {p1*100:<12.2f}% {p2*100:<12.2f}% {tot*100:<12.2f}% {ratio:<12.2f}%")
