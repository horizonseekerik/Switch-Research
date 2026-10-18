"""
PROJECT JANUS: HIGH-SPEED DUAL-STATE RESONANCE SWEEP FOR PROTECTED SWITCH
========================================================================
Sweeps L_c across [2.80, 3.80] um with step 0.20 um (or user-specified)
Evaluates BOTH Amorphous (State 0) and Crystalline (State 1) states.
Finds the exact operating point where both states reach maximum transmission.
"""

import sys, os, math, time
import numpy as np
import meep as mp

from test_protected_sub6um_dc import simulate_protected_dc

def run_sweep(lc_list=[2.80, 3.00, 3.20, 3.40, 3.60, 3.80], gap=0.098, resolution=15):
    print("=" * 95)
    print("  PROJECT JANUS: DUAL-STATE PROTECTED COUPLER SWEEP")
    print(f"  Gap = {gap*1000:.1f} nm, Resolution = {resolution}")
    print("=" * 95)
    print(f"{'L_c(um)':<8} {'AM P2(Cross)%':<15} {'AM P1(Bar)%':<13} {'CR P1(Bar)%':<13} {'CR P2(Cross)%':<15} {'AM XT(dB)':<11} {'CR XT(dB)':<11} {'Area(um2)':<10}")
    print("-" * 95)

    # 1. Cache reference flux once using a standard cell
    res_ref = simulate_protected_dc("amorphous", resolution=resolution, l_c=3.40, gap=gap, verbose=False)
    p_ref = res_ref["p_ref"]

    best_lc = None
    best_score = -1e9
    results = []

    for lc in lc_list:
        res_am = simulate_protected_dc("amorphous", resolution=resolution, l_c=lc, gap=gap, p_ref=p_ref, verbose=False)
        res_cr = simulate_protected_dc("crystalline", resolution=resolution, l_c=lc, gap=gap, p_ref=p_ref, verbose=False)

        am_p2 = res_am["p2"] * 100.0
        am_p1 = res_am["p1"] * 100.0
        cr_p1 = res_cr["p1"] * 100.0
        cr_p2 = res_cr["p2"] * 100.0
        am_xt = res_am["xt_db"]
        cr_xt = res_cr["xt_db"]
        area = res_am["cell_area"]

        print(f"{lc:<8.2f} {am_p2:<15.2f} {am_p1:<13.2f} {cr_p1:<13.2f} {cr_p2:<15.2f} {am_xt:<11.2f} {cr_xt:<11.2f} {area:<10.3f}")
        sys.stdout.flush()

        # Score function: penalizes lower of the two signal transmissions
        score = min(am_p2, cr_p1)
        if score > best_score:
            best_score = score
            best_lc = lc

        results.append({
            "lc": lc, "am_p2": am_p2, "am_p1": am_p1,
            "cr_p1": cr_p1, "cr_p2": cr_p2,
            "am_xt": am_xt, "cr_xt": cr_xt, "area": area
        })

    print("=" * 95)
    print(f"Optimal Dual-State Operating Point: L_c = {best_lc:.2f} um (min signal = {best_score:.2f}%)")
    print("=" * 95)

if __name__ == '__main__':
    run_sweep()
