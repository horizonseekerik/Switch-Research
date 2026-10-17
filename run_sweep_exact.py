import numpy as np
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

from find_exact_switch_design import get_isolated_wg_neff

w2 = 0.240
with nostdout():
    neff_w2 = get_isolated_wg_neff(w2, t_pcm=0.0)
print(f"Bare WG2 ({w2*1000:.0f} nm): neff = {neff_w2:.5f}")

for tp in [0.080, 0.100, 0.120, 0.150, 0.180, 0.200]:
    best_w1 = None
    min_d = 999.0
    best_cr = None
    for w1 in np.linspace(0.100, 0.200, 26):
        with nostdout():
            n_cr = get_isolated_wg_neff(w1, t_pcm=tp, pcm_state="crystalline")
        d = abs(n_cr - neff_w2)
        if d < min_d:
            min_d = d
            best_w1 = w1
            best_cr = n_cr
    with nostdout():
        n_am = get_isolated_wg_neff(best_w1, t_pcm=tp, pcm_state="amorphous")
    d_beta = (2 * math.pi / 1.064) * abs(n_cr - n_am)
    kappa = d_beta / (2 * math.sqrt(3))
    lc = math.pi / (2 * kappa)
    area = (lc + 0.35 + 0.60 + 0.35) * (0.45 + 0.24)
    print(f"tp={tp*1000:3.0f}nm: w1={best_w1*1000:5.1f}nm (diff={min_d:.5f}), n_cr={best_cr:.4f}, n_am={n_am:.4f}, dbeta={d_beta:.4f}, lc={lc:.2f}um, Area={area:.2f}um^2")
