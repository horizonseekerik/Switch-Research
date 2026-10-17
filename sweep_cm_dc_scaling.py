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
print("COUPLING-MODULATED 2:1 RATIO SWEEP (VARYING GAP & t_pcm)")
print("Target: Find gap & t_pcm where Ratio = L_c(am)/L_c(cr) == 2.000 AND L_c < 4 um")
print("=" * 75)
print(f"{'t_pcm(nm)':<10} {'gap(nm)':<10} {'L_c(am) um':<12} {'L_c(cr) um':<12} {'Ratio':<10} {'Area um^2':<10}")
print("-" * 75)

for tp in [0.080, 0.100, 0.120, 0.150]:
    for g in [0.015, 0.020, 0.025, 0.030, 0.035]:
        with nostdout():
            k1_am, k2_am = solve_symmetric_pcm_coupler(w=0.240, gap=g, t_pcm=tp, state="amorphous")
            k1_cr, k2_cr = solve_symmetric_pcm_coupler(w=0.240, gap=g, t_pcm=tp, state="crystalline")
        dn_am = k1_am - k2_am
        dn_cr = k1_cr - k2_cr
        lc_am = 1.064 / (2.0 * dn_am)
        lc_cr = 1.064 / (2.0 * dn_cr)
        ratio = lc_am / lc_cr
        # Area = (L_c + 0.35 + 0.50 + 0.35) * (0.45 + 0.16)
        area = (lc_am + 1.20) * 0.61
        print(f"{tp*1000:<10.0f} {g*1000:<10.0f} {lc_am:<12.2f} {lc_cr:<12.2f} {ratio:<10.3f} {area:<10.2f}")
