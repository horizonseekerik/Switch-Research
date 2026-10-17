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

from find_exact_switch_design import get_isolated_wg_neff

print("=" * 80)
print("PCM THICKNESS vs FOOTPRINT OPTIMIZATION (SWEEPING t_pcm & WIDTH ASYMMETRY)")
print("=" * 80)

for w2 in [0.220, 0.240, 0.250]:
    with nostdout():
        neff_w2 = get_isolated_wg_neff(w2, t_pcm=0.0)
    print(f"\n[Bare WG2 = {w2*1000:.0f} nm | n_eff = {neff_w2:.5f}]")
    for t_pcm in [0.040, 0.060, 0.080, 0.100, 0.120, 0.150]:
        min_diff = 999.0
        best_w1 = None
        best_neff_am = None
        for w1 in np.linspace(w2 - 0.080, w2, 21):
            with nostdout():
                neff_am = get_isolated_wg_neff(w1, t_pcm=t_pcm, pcm_state="amorphous")
            diff = abs(neff_am - neff_w2)
            if diff < min_diff:
                min_diff = diff
                best_w1 = w1
                best_neff_am = neff_am
        with nostdout():
            best_neff_cr = get_isolated_wg_neff(best_w1, t_pcm=t_pcm, pcm_state="crystalline")
        dn = best_neff_cr - best_neff_am
        dbeta = (2 * math.pi / 1.064) * dn
        kappa = dbeta / (2.0 * math.sqrt(3))
        lc = math.pi / (2.0 * kappa)
        cell_w = (w2 * 2 + 0.050) + 0.240 # pitch + 2*isolation
        cell_l = lc + 0.20 + 0.45 + 0.20 # Lin + Lc + Lbend + Lout
        area = cell_l * cell_w
        print(f"  t_pcm={t_pcm*1000:3.0f}nm | w1={best_w1*1000:5.1f}nm | dn={dn:6.4f} | Lc={lc:5.2f}um | Ltot={cell_l:5.2f}um | Area={area:5.2f}um^2")
