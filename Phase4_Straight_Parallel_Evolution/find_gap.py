import os
import sys
import math
import contextlib
import numpy as np

@contextlib.contextmanager
def nostdout():
    with open(os.devnull, "w") as devnull:
        old = sys.stdout
        sys.stdout = devnull
        try:
            yield
        finally:
            sys.stdout = old

from solve_switch_modes import solve_coupler_modes

target_lc = 6.04
print(f"Searching for gap that gives L_c = {target_lc} um at w_avg = 212 nm, tpcm = 100 nm:")
for gap in [0.030, 0.045, 0.060, 0.080, 0.100, 0.120, 0.140, 0.160]:
    with nostdout():
        res = solve_coupler_modes(w_wg=0.212, h_wg=0.220, gap=gap, t_pcm=0.100, pcm_loc="top", pcm_state="amorphous")
    lc = res["l_c_um"]
    kap = res["kappa_um"]
    print(f"  gap = {gap*1000:3.0f} nm -> L_c = {lc:5.2f} um (kappa = {kap:6.4f} um^-1)")
