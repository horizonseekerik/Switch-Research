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

from find_exact_switch_design import get_isolated_wg_neff

print("=" * 75)
print("EXACT 1x2 DIRECTIONAL COUPLER DESIGN (L_c = 4.60 um, t_pcm = 150 nm)")
print("Condition: Delta_beta / 2*kappa == sqrt(3) -> Extinction > 30 dB")
print("=" * 75)

# Target kappa for L_c = 4.60 um
target_lc = 4.60
target_kappa = math.pi / (2.0 * target_lc)
target_dbeta = 2.0 * math.sqrt(3) * target_kappa
target_dn = target_dbeta * 1.064 / (2.0 * math.pi)

print(f"Target L_c         : {target_lc:.2f} um")
print(f"Target kappa       : {target_kappa:.4f} um^-1")
print(f"Target Delta_beta  : {target_dbeta:.4f} um^-1")
print(f"Target Delta_n_eff : {target_dn:.5f}")

# Sweep w1 and w2 to find exact match
print("\nSearching for (w1, w2) matching Delta_n and phase matching in Amorphous...")
w2_candidates = [0.220, 0.230, 0.240, 0.250]
tp = 0.150

for w2 in w2_candidates:
    with nostdout():
        neff_w2 = get_isolated_wg_neff(w2, t_pcm=0.0)
    # Search w1 that matches neff_w2 in amorphous
    best_w1 = None
    min_diff = 999.0
    for w1 in np.linspace(0.120, 0.200, 33):
        with nostdout():
            neff_am = get_isolated_wg_neff(w1, t_pcm=tp, pcm_state="amorphous")
        d = abs(neff_am - neff_w2)
        if d < min_diff:
            min_diff = d
            best_w1 = w1
    with nostdout():
        neff_am = get_isolated_wg_neff(best_w1, t_pcm=tp, pcm_state="amorphous")
        neff_cr = get_isolated_wg_neff(best_w1, t_pcm=tp, pcm_state="crystalline")
    dn = neff_cr - neff_am
    ratio = (2.0 * math.pi / 1.064 * dn) / (2.0 * target_kappa)
    print(f"w2 = {w2*1000:.0f} nm: Matched w1 = {best_w1*1000:.1f} nm, n_am = {neff_am:.4f}, n_cr = {neff_cr:.4f}, Delta_n = {dn:.5f}, Ratio dbeta/2kappa = {ratio:.3f}")
