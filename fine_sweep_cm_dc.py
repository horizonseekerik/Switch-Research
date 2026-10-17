import math
import numpy as np
import meep as mp
from meep import mpb
import contextlib, os, sys

@contextlib.contextmanager
def nostdout():
    with open(os.devnull, "w") as devnull:
        old = sys.stdout; sys.stdout = devnull
        try: yield
        finally: sys.stdout = old

from test_symmetric_cm_dc import solve_symmetric_pcm_coupler

print("=" * 75)
print("FINE SWEEP FOR EXACT RATIO == 2.000 (SUB-4 um^2 CM-DC SWITCH)")
print("=" * 75)
print(f"{'w(nm)':<8} {'t_pcm(nm)':<10} {'gap(nm)':<9} {'L_c(am)':<10} {'L_c(cr)':<10} {'Ratio':<9} {'Area um^2':<10}")
print("-" * 75)

for w in [0.225, 0.230, 0.235, 0.240]:
    for tp in [0.120, 0.140, 0.150, 0.160]:
        for g in [0.020, 0.025, 0.030]:
            with nostdout():
                k1_am, k2_am = solve_symmetric_pcm_coupler(w=w, gap=g, t_pcm=tp, state="amorphous")
                k1_cr, k2_cr = solve_symmetric_pcm_coupler(w=w, gap=g, t_pcm=tp, state="crystalline")
            dn_am = k1_am - k2_am
            dn_cr = k1_cr - k2_cr
            lc_am = 1.064 / (2.0 * dn_am)
            lc_cr = 1.064 / (2.0 * dn_cr)
            ratio = lc_am / lc_cr
            area = (lc_am + 0.95) * 0.61
            if 1.95 <= ratio <= 2.05:
                print(f"{w*1000:<8.1f} {tp*1000:<10.0f} {g*1000:<9.0f} {lc_am:<10.2f} {lc_cr:<10.2f} {ratio:<9.3f} {area:<10.2f}  <-- EXACT MATCH")
