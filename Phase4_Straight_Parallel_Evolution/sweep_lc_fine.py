import sys, os, math
from test_coupler_geometry_clean import simulate_switch

print("==================================================================")
print("  FINE RESONANCE SEARCH ACROSS L_c in [1.40, 1.90] um (gap = 45 nm)")
print("==================================================================")
print(f"{'L_c (um)':<10} {'P1 (Bar)%':<12} {'P2 (Cross)%':<14} {'Total%':<10} {'Rad Loss%':<12} {'XT (dB)':<10} {'Area (um2)':<10}")
print("-" * 80)

best_p2 = 0
best_lc = None

for lc in [1.40, 1.50, 1.60, 1.70, 1.80, 1.90]:
    res = simulate_switch("amorphous", "si", resolution=15, l_c_override=lc, gap_override=0.045)
    p1 = res["p1"] * 100
    p2 = res["p2"] * 100
    ptot = (res["p1"] + res["p2"]) * 100
    rad = (1.0 - (res["p1"] + res["p2"])) * 100
    xt = res["xt_db"]
    area = res["cell_area"]
    print(f"{lc:<10.2f} {p1:<12.2f} {p2:<14.2f} {ptot:<10.2f} {rad:<12.2f} {xt:<10.2f} {area:<10.3f}")
    if p2 > best_p2:
        best_p2 = p2
        best_lc = lc

print("=" * 80)
print(f"Optimal L_c = {best_lc:.2f} um with P2 = {best_p2:.2f}%")
