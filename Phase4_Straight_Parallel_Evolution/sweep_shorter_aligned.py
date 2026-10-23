import sys, os, time
from test_aligned_detuned_switch import simulate_aligned_switch

print("Checking L_c = [3.00, 3.20, 3.35, 3.45, 3.55] um at gap = 89 nm...")
gap = 0.089
t_pcm = 0.200
res = 15

res_ref = simulate_aligned_switch("amorphous", resolution=res, l_c=3.35, gap=gap, t_pcm=t_pcm, verbose=False)
p_ref_flux = res_ref["p_ref_flux"]
p_ref_mode = res_ref["p_ref_mode"]

print(f"{'L_c(um)':<8} {'AM P2(Mod)%':<14} {'AM P1(Mod)%':<13} {'CR P1(Mod)%':<13} {'CR P2(Mod)%':<14} {'AM Sig/G%':<12} {'CR Sig/G%':<12}")
print("-" * 90)

for lc in [3.00, 3.15, 3.30, 3.45]:
    res_am = simulate_aligned_switch("amorphous", resolution=res, l_c=lc, gap=gap, t_pcm=t_pcm,
                                     p_ref_flux=p_ref_flux, p_ref_mode=p_ref_mode, verbose=False)
    res_cr = simulate_aligned_switch("crystalline", resolution=res, l_c=lc, gap=gap, t_pcm=t_pcm,
                                     p_ref_flux=p_ref_flux, p_ref_mode=p_ref_mode, verbose=False)

    am_p2 = res_am["p2_mode"] * 100.0
    am_p1 = res_am["p1_mode"] * 100.0
    cr_p1 = res_cr["p1_mode"] * 100.0
    cr_p2 = res_cr["p2_mode"] * 100.0
    am_ratio = res_am["ratio_guided"]
    cr_ratio = res_cr["ratio_guided"]

    print(f"{lc:<8.2f} {am_p2:<14.2f} {am_p1:<13.2f} {cr_p1:<13.2f} {cr_p2:<14.2f} {am_ratio:<12.2f} {cr_ratio:<12.2f}")
    sys.stdout.flush()
