"""
PROJECT JANUS: STAGGERED-PCM DIRECTIONAL COUPLER SWEEP
======================================================
Systematic sweep over L_c and gap for the staggered-PCM buffer-zone switch architecture.
All quantities recorded and reported to >= 10 decimal places.
Pure TE0 modal decomposition via MPB eigenmode coefficients.
"""

import sys, os, math, time, argparse
import numpy as np
import meep as mp

from test_staggered_pcm_switch import simulate_staggered_pcm

def run_sweep(l_c_list, gap_list, l_buffer=0.350, resolution=15):
    print("="*100)
    print("STARTING STAGGERED-PCM 3D FDTD SWEEP")
    print(f"L_c values (um)   : {[f'{x:.10f}' for x in l_c_list]}")
    print(f"Gap values (um)   : {[f'{g:.10f}' for g in gap_list]}")
    print(f"Buffer length (um): {l_buffer:.10f}")
    print(f"Resolution        : {resolution}")
    print("="*100)

    # First get reference straight waveguide to reuse across all simulations
    print("\nComputing Master Reference Waveguide...")
    init_res = simulate_staggered_pcm(
        state="amorphous", resolution=resolution,
        l_c=l_c_list[0], gap=gap_list[0], l_buffer=l_buffer,
        verbose=False
    )
    p_ref_flux = init_res["p_ref_flux"]
    p_ref_mode = init_res["p_ref_mode"]
    print(f"Master p_ref_flux = {p_ref_flux:.10e}")
    print(f"Master p_ref_mode = {p_ref_mode:.10e}\n")

    results = []
    total_runs = len(l_c_list) * len(gap_list)
    run_idx = 0

    for gap in gap_list:
        for l_c in l_c_list:
            run_idx += 1
            print(f"\n--- Point {run_idx}/{total_runs}: L_c = {l_c:.10f} um, Gap = {gap*1000.0:.10f} nm ---")
            t0 = time.time()
            
            # State 0: Amorphous (cross port = Port 2)
            r0 = simulate_staggered_pcm(
                state="amorphous", resolution=resolution,
                l_c=l_c, gap=gap, l_buffer=l_buffer,
                p_ref_flux=p_ref_flux, p_ref_mode=p_ref_mode,
                verbose=False
            )
            
            # State 1: Crystalline (bar port = Port 1)
            r1 = simulate_staggered_pcm(
                state="crystalline", resolution=resolution,
                l_c=l_c, gap=gap, l_buffer=l_buffer,
                p_ref_flux=p_ref_flux, p_ref_mode=p_ref_mode,
                verbose=False
            )
            elapsed = time.time() - t0
            
            p2_s0 = r0["p2_mode"]
            p1_s0 = r0["p1_mode"]
            ptot_s0 = r0["ptot_mode"]
            il_s0 = r0["il_db"]
            xt_s0 = r0["xt_db"]

            p1_s1 = r1["p1_mode"]
            p2_s1 = r1["p2_mode"]
            ptot_s1 = r1["ptot_mode"]
            il_s1 = r1["il_db"]
            xt_s1 = r1["xt_db"]

            area = r0["cell_area"]
            
            point_data = {
                "l_c": l_c, "gap": gap, "l_buffer": l_buffer,
                "p2_s0": p2_s0, "p1_s0": p1_s0, "ptot_s0": ptot_s0, "il_s0": il_s0, "xt_s0": xt_s0,
                "p1_s1": p1_s1, "p2_s1": p2_s1, "ptot_s1": ptot_s1, "il_s1": il_s1, "xt_s1": xt_s1,
                "cell_area": area, "elapsed_s": elapsed
            }
            results.append(point_data)

            print(f"  [S0 Amorphous ] P2 (target) = {p2_s0*100.0:14.10f}%, P1 (leak) = {p1_s0*100.0:14.10f}%, Ptot = {ptot_s0*100.0:14.10f}%, IL = {il_s0:10.10f} dB, XT = {xt_s0:10.10f} dB")
            print(f"  [S1 Crystallin] P1 (target) = {p1_s1*100.0:14.10f}%, P2 (leak) = {p2_s1*100.0:14.10f}%, Ptot = {ptot_s1*100.0:14.10f}%, IL = {il_s1:10.10f} dB, XT = {xt_s1:10.10f} dB")
            print(f"  Footprint: {area:.10f} um^2 | Elapsed: {elapsed:.2f} s")

    print("\n" + "="*120)
    print("SWEEP SUMMARY TABLE (Sorted by min(P2_s0, P1_s1)):")
    print("="*120)
    print(f"{'L_c (um)':<14} | {'Gap (nm)':<12} | {'S0 P2 (%)':<18} | {'S1 P1 (%)':<18} | {'S0 Ptot (%)':<16} | {'S1 Ptot (%)':<16} | {'Min Target (%)':<18} | {'Area (um^2)':<12}")
    print("-" * 120)

    # Sort by minimum of target powers
    results.sort(key=lambda d: min(d["p2_s0"], d["p1_s1"]), reverse=True)

    for d in results:
        min_target = min(d["p2_s0"], d["p1_s1"]) * 100.0
        print(f"{d['l_c']:<14.10f} | {d['gap']*1000.0:<12.6f} | {d['p2_s0']*100.0:<18.10f} | {d['p1_s1']*100.0:<18.10f} | {d['ptot_s0']*100.0:<16.10f} | {d['ptot_s1']*100.0:<16.10f} | {min_target:<18.10f} | {d['cell_area']:<12.10f}")

    print("="*120)
    return results

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--l_c", nargs="+", type=float, default=[4.30, 4.40, 4.42, 4.50, 4.60])
    parser.add_argument("--gap", nargs="+", type=float, default=[0.096, 0.098, 0.100])
    parser.add_argument("--buffer", type=float, default=0.350)
    parser.add_argument("--res", type=int, default=15)
    args = parser.parse_args()

    run_sweep(args.l_c, args.gap, l_buffer=args.buffer, resolution=args.res)
