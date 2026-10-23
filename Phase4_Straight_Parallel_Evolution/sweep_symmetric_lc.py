import sys, os, math
from test_symmetric_stack import simulate_switch

print("=" * 80)
print("  RESONANCE SWEEP FOR SYMMETRIC Sb2S3 STACK (L_c in [0.80, 1.50] um, gap = 45 nm)")
print("=" * 80)
print(f"{'L_c (um)':<10} {'P1 (Bar)%':<12} {'P2 (Cross)%':<14} {'Total%':<10} {'Rad Loss%':<12} {'XT (dB)':<10} {'Area (um2)':<10}")
print("-" * 80)

best_p2 = 0
best_lc = None

for lc in [1.40, 1.30, 1.20, 1.10, 1.00, 0.90, 0.80]:
    res = simulate_switch("amorphous", resolution=15, l_c=lc, symmetric_pcm=True)
    p1 = res["p1"] * 100
    p2 = res["p2"] * 100
    ptot = (res["p1"] + res["p2"]) * 100
    rad = (1.0 - (res["p1"] + res["p2"])) * 100
    xt = res["xt_db"]
    area = res["cell_area"]
    print(f"{lc:<10.2f} {p1:<12.2f} {p2:<14.2f} {ptot:<10.2f} {rad:<12.2f} {xt:<10.2f} {area:<10.3f}")
    sys.stdout.flush()
    if p2 > best_p2:
        best_p2 = p2
        best_lc = lc

print("=" * 80)
print(f"Optimal Symmetric L_c = {best_lc:.2f} um with P2 = {best_p2:.2f}%")
