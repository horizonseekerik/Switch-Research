"""
PROJECT JANUS: CANDIDATE 3 GAP SCAN FOR SIMULTANEOUS STATE 0 & STATE 1 RESONANCE
================================================================================
Physics:
  State 0 half-beat scales as 1/kappa (moves ~100 nm per 1 nm gap change).
  State 1 Kogelnik null scales as 1/sqrt(kappa^2 + (dbeta/2)^2) where dbeta dominates
  (moves only ~9 nm per 1 nm gap change).
  
  Scanning gap from 100 to 109 nm will sweep State 0's peak across State 1's null
  at L_c ~ 4.15 - 4.30 um.
"""

import sys, os, math, time, argparse
import numpy as np
import meep as mp
from test_candidate3_optimal import simulate_candidate3

def run_gap_scan():
    t_pcm = 0.190
    res = 15
    
    # 2 gap values, 2 L_c values for focused initial scan
    gap_list = [0.102, 0.106]
    lc_list = [4.18, 4.30]

    print("=" * 125)
    print("  PROJECT JANUS: CANDIDATE 3 GAP SCAN (t_pcm = 190 nm, straight bus, low loss)")
    print(f"  Gap values (nm): {[f'{g*1e3:.1f}' for g in gap_list]}")
    print(f"  L_c values (um): {[f'{l:.2f}' for l in lc_list]}")
    print(f"  Resolution     : {res}")
    print("=" * 125)
    print(f"{'Gap(nm)':<10} {'L_c(um)':<10} {'S0 P2(Cross)%':<16} {'S0 P1(Leak)%':<15} {'S1 P1(Bar)%':<15} {'S1 P2(Leak)%':<16} {'S0 Tot%':<12} {'S1 Tot%':<12} {'Min Signal%':<14}")
    print("-" * 125)

    # Master reference waveguide
    res_ref = simulate_candidate3("amorphous", resolution=res, l_c=4.25, gap=0.103, t_pcm=t_pcm, verbose=False)
    p_ref_flux = res_ref["p_ref_flux"]
    p_ref_mode = res_ref["p_ref_mode"]

    results = []

    for gap in gap_list:
        for lc in lc_list:
            r0 = simulate_candidate3("amorphous", resolution=res, l_c=lc, gap=gap, t_pcm=t_pcm,
                                     p_ref_flux=p_ref_flux, p_ref_mode=p_ref_mode, verbose=False)
            r1 = simulate_candidate3("crystalline", resolution=res, l_c=lc, gap=gap, t_pcm=t_pcm,
                                     p_ref_flux=p_ref_flux, p_ref_mode=p_ref_mode, verbose=False)

            p2_s0 = r0["p2_mode"] * 100.0
            p1_s0 = r0["p1_mode"] * 100.0
            tot_s0 = r0["ptot_mode"] * 100.0

            p1_s1 = r1["p1_mode"] * 100.0
            p2_s1 = r1["p2_mode"] * 100.0
            tot_s1 = r1["ptot_mode"] * 100.0

            min_sig = min(p2_s0, p1_s1)

            entry = {
                "gap": gap, "lc": lc, "p2_s0": p2_s0, "p1_s0": p1_s0, "tot_s0": tot_s0,
                "p1_s1": p1_s1, "p2_s1": p2_s1, "tot_s1": tot_s1, "min_sig": min_sig
            }
            results.append(entry)

            print(f"{gap*1e3:<10.1f} {lc:<10.2f} {p2_s0:<16.10f} {p1_s0:<15.10f} {p1_s1:<15.10f} {p2_s1:<16.10f} {tot_s0:<12.6f} {tot_s1:<12.6f} {min_sig:<14.10f}")
            sys.stdout.flush()

    print("=" * 125)
    results.sort(key=lambda d: d["min_sig"], reverse=True)
    best = results[0]
    print(f"\nTOP CONFIGURATION: Gap = {best['gap']*1e3:.1f} nm, L_c = {best['lc']:.2f} um")
    print(f"  State 0 Cross (P2) = {best['p2_s0']:.10f}%, Leakage (P1) = {best['p1_s0']:.10f}%")
    print(f"  State 1 Bar   (P1) = {best['p1_s1']:.10f}%, Leakage (P2) = {best['p2_s1']:.10f}%")
    print(f"  Worst-Case Target   = {best['min_sig']:.10f}%")
    print("=" * 125)

if __name__ == '__main__':
    run_gap_scan()
