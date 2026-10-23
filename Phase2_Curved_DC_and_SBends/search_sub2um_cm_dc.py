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

print("=" * 80)
print("SEARCHING FOR SUB-2 um^2 / SUB-3 um^2 COUPLING-MODULATED SWITCHES")
print("Target: Find (w, gap, t_pcm) where Ratio == 2.000 AND L_c(am) < 3.5 um")
print("=" * 80)
print(f"{'w(nm)':<8} {'t_pcm(nm)':<10} {'gap(nm)':<9} {'L_c(am)':<10} {'L_c(cr)':<10} {'Ratio':<9} {'Area um^2':<10}")
print("-" * 80)

for w in [0.180, 0.200, 0.220]:
    for tp in [0.100, 0.120, 0.150]:
        for g in [0.020, 0.025, 0.030, 0.035, 0.040]:
            with nostdout():
                k1_am, k2_am = solve_symmetric_pcm_coupler(w=w, gap=g, t_pcm=tp, state="amorphous")
                k1_cr, k2_cr = solve_symmetric_pcm_coupler(w=w, gap=g, t_pcm=tp, state="crystalline")
            dn_am = k1_am - k2_am
            dn_cr = k1_cr - k2_cr
            lc_am = 1.064 / (2.0 * dn_am)
            lc_cr = 1.064 / (2.0 * dn_cr)
            ratio = lc_am / lc_cr
            # Estimated footprint area: (L_c + 0.25 + 0.45 + 0.25) * (0.45 + 0.16)
            area = (lc_am + 0.95) * 0.61
            if 1.80 <= ratio <= 2.20 or lc_am < 4.0:
                print(f"{w*1000:<8.0f} {tp*1000:<10.0f} {g*1000:<9.0f} {lc_am:<10.2f} {lc_cr:<10.2f} {ratio:<9.3f} {area:<10.2f}")
