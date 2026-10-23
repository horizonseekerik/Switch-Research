import math
import contextlib
import os
import sys

@contextlib.contextmanager
def nostdout():
    with open(os.devnull, "w") as devnull:
        old = sys.stdout; sys.stdout = devnull
        try: yield
        finally: sys.stdout = old

from test_symmetric_cm_dc import solve_symmetric_pcm_coupler

print("=" * 60)
print("TESTING LARGER GAPS FOR CLEAN RESOLUTION (gap >= 50 nm)")
print("=" * 60)
print(f"{'gap (nm)':<10} {'L_c(am) um':<12} {'L_c(cr) um':<12} {'Ratio':<10}")
print("-" * 60)
for g in [0.040, 0.050, 0.060, 0.070, 0.080]:
    with nostdout():
        k1_am, k2_am = solve_symmetric_pcm_coupler(w=0.240, gap=g, t_pcm=0.150, state="amorphous")
        k1_cr, k2_cr = solve_symmetric_pcm_coupler(w=0.240, gap=g, t_pcm=0.150, state="crystalline")
    lc_am = 1.064 / (2.0 * (k1_am - k2_am))
    lc_cr = 1.064 / (2.0 * (k1_cr - k2_cr))
    print(f"{g*1000:<10.0f} {lc_am:<12.2f} {lc_cr:<12.2f} {lc_am/lc_cr:<10.3f}")
