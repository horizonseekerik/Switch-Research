"""
PROJECT JANUS: PATH 2 - PARALLEL GAP (89 nm) WITH WIDENED WG2 (w2_end = 340 nm)
================================================================================
Investigates keeping the inner gap fixed at 89 nm while widening WG2 outward to 340 nm.
Eliminates premature coupling cutoff for State 1 while maintaining strong detuning for State 0.
"""

import sys, os, math, time
import numpy as np
import meep as mp

from test_path2_optimized_dc import simulate_path2_optimized

if __name__ == "__main__":
    t0 = time.time()
    # gap_s1 = 89 nm, gap_s2_end = 89 nm (parallel inner edge)
    print("Testing Path 2 with fixed inner gap = 89 nm and w2_end = 340 nm...")
    res_am = simulate_path2_optimized(state="amorphous", resolution=15,
                                      l_1=3.00, l_2=1.25, w2_end=0.340,
                                      gap_s1=0.089, gap_s2_end=0.089)
    res_cr = simulate_path2_optimized(state="crystalline", resolution=15,
                                      l_1=3.00, l_2=1.25, w2_end=0.340,
                                      gap_s1=0.089, gap_s2_end=0.089,
                                      p_ref_flux=res_am["p_ref_flux"], p_ref_mode=res_am["p_ref_mode"])
    print(f"\nElapsed run time = {time.time() - t0:.4f} s\n")
