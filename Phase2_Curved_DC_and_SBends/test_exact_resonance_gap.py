import sys, math
import numpy as np

# Exact coupled mode solver at u = sqrt(3)
k_89 = 0.5400000000
dbeta_eff = 0.9780000000

print("=" * 105)
print("  EVALUATING EXACT RESONANCE GAP FOR S/kappa = 2.0000000000 (10 DECIMAL PLACES)")
print("=" * 105)
print(f"{'Gap (nm)':<12} {'L_c (um)':<16} {'AM P2%':<22} {'AM P1%':<22} {'CR P1%':<22} {'CR P2%':<22}")
print("-" * 105)

for g in np.arange(135.0, 142.0, 0.5):
    k = k_89 * math.exp(-(g - 89.0) / 75.0)
    u = dbeta_eff / (2.0 * k)
    s = k * math.sqrt(1.0 + u**2)
    
    # Optimum length where both states align:
    lc = math.pi / (2.0 * k)
    
    # State 0 (Amorphous):
    phi_am = k * lc
    p2_am = math.sin(phi_am)**2
    p1_am = math.cos(phi_am)**2
    
    # State 1 (Crystalline):
    # Cross power = 1/(1+u^2) * sin^2(s * lc)
    phi_cr = s * lc
    p2_cr = (1.0 / (1.0 + u**2)) * (math.sin(phi_cr)**2)
    p1_cr = 1.0 - p2_cr
    
    print(f"{g:<12.5f} {lc:<16.10f} {p2_am*100:<22.10f} {p1_am*100:<22.10f} {p1_cr*100:<22.10f} {p2_cr*100:<22.10f}")

print("=" * 105)
