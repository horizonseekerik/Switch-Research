"""
TEST STATE 1 (CRYSTALLINE) AT L_c = 3.60 um (GAP = 98.14 nm)
============================================================
Checks whether State 1 is already close to null at L_c = 3.60 um.
"""

import sys, os, math, time
import numpy as np

from test_symmetric_pcm_aligned_3d import simulate_symmetric_dc

if __name__ == "__main__":
    t0 = time.time()
    lc = 3.60
    gap = 0.09814
    print("=" * 105)
    print(f"TESTING STATE 1 AT L_c = {lc:.2f} um, Gap = {gap*1e3:.2f} nm")
    print("=" * 105)

    res_cr = simulate_symmetric_dc(state="crystalline", resolution=15, l_c=lc, gap=gap, verbose=True)
    print(f"\nElapsed time = {time.time() - t0:.2f} s")
