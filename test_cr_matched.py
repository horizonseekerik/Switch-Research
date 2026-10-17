import os
import sys
import math
import contextlib
import numpy as np

@contextlib.contextmanager
def nostdout():
    with open(os.devnull, "w") as devnull:
        old = sys.stdout; sys.stdout = devnull
        try: yield
        finally: sys.stdout = old

from find_exact_switch_design import get_isolated_wg_neff

print("=" * 80)
print("RESONANT PHASE-MATCHED SWITCH DESIGN (CR-MATCHED)")
print("Target: Amorphous = Detuned (Bar), Crystalline = Resonant (Cross)")
print("=" * 80)

w2 = 0.250
with nostdout():
    neff_w2 = get_isolated_wg_neff(w2, t_pcm=0.0)
print(f"Bare WG2 (w2 = {w2*1000:.0f} nm): n_eff = {neff_w2:.5f}")

for t_pcm in [0.040, 0.060, 0.080]:
    min_diff = 999.0
    best_w1 = None
    best_neff_cr = None
    for w1 in np.linspace(0.120, 0.220, 21):
        with nostdout():
            neff_cr = get_isolated_wg_neff(w1, t_pcm=t_pcm, pcm_state="crystalline")
        diff = abs(neff_cr - neff_w2)
        if diff < min_diff:
            min_diff = diff
            best_w1 = w1
            best_neff_cr = neff_cr
    with nostdout():
        neff_am = get_isolated_wg_neff(best_w1, t_pcm=t_pcm, pcm_state="amorphous")
    delta_beta_am = (2 * math.pi / 1.064) * abs(neff_am - neff_w2)
    print(f"\nt_pcm = {t_pcm*1000:.0f} nm:")
    print(f"  Matched w1 in Crystalline: {best_w1*1000:.1f} nm (n_eff_cr = {best_neff_cr:.5f}, diff = {min_diff:.5f})")
    print(f"  Amorphous n_eff          : {neff_am:.5f}")
    print(f"  Amorphous Detuning dBeta : {delta_beta_am:.4f} um^-1")
