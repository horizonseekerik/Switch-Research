import sys, math
import numpy as np

# Let's inspect the actual FDTD extracted parameters:
# At gap = 89nm, L_0 = 2.9 um gives theta = pi/2 -> kappa_eff = 0.54 rad/um
# At gap = 89nm, L_1 = 4.3 um gives theta = pi -> S_eff = 0.73 rad/um
# Ratio S/k = 1.35, u = 0.91

# Now let's calculate what gap makes u_eff = sqrt(3) = 1.732:
# In silicon couplers, kappa(g) decays exponentially: kappa(g) = kappa_0 * exp(-g / g_0)
# where g_0 ~ 85 nm.

print("=" * 80)
print("  PHYSICAL COUPLING GAP ANALYSIS: SOLVING S/kappa = 2.000")
print("=" * 80)
print(f"{'Gap (nm)':<10} {'kappa_eff':<12} {'dbeta_eff':<12} {'u = dbeta/(2k)':<16} {'S/kappa':<12} {'Both >98% Possible?'}")
print("-" * 80)

# At gap = 89nm, kappa_eff = 0.54 rad/um, dbeta_eff = 2 * 0.54 * 0.906 = 0.978 rad/um
k_89 = 0.54
dbeta_eff = 0.978

for g_nm in [89, 95, 100, 105, 110, 115, 120, 125, 130, 135]:
    g = g_nm
    # exponential decay of kappa with gap
    k = k_89 * math.exp(-(g - 89) / 75.0)
    u = dbeta_eff / (2.0 * k)
    s_over_k = math.sqrt(1.0 + u**2)
    
    # Check if S/k is close to 2.0
    match = "YES! EXACT!" if abs(s_over_k - 2.0) < 0.05 else ("Near" if abs(s_over_k - 2.0) < 0.15 else "No")
    l_am = math.pi / (2.0 * k)
    l_cr = math.pi / (k * s_over_k)
    dl = abs(l_cr - l_am)
    
    print(f"{g_nm:<10} {k:<12.4f} {dbeta_eff:<12.4f} {u:<16.3f} {s_over_k:<12.3f} {match} (L_am={l_am:.2f}um, L_cr={l_cr:.2f}um, dL={dl:.2f}um)")

print("=" * 80)
