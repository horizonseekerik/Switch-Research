from test_aligned_detuned_switch import simulate_aligned_switch

print("Running single point check at L_c = 4.40 um, gap = 98.2 nm with verbose=True...")
res = simulate_aligned_switch("amorphous", resolution=15, l_c=4.40, gap=0.0982, t_pcm=0.200,
                              y_out=0.420, l_bend=0.900, w_bend=0.320, verbose=True)

res_cr = simulate_aligned_switch("crystalline", resolution=15, l_c=4.40, gap=0.0982, t_pcm=0.200,
                                 y_out=0.420, l_bend=0.900, w_bend=0.320,
                                 p_ref_flux=res["p_ref_flux"], p_ref_mode=res["p_ref_mode"], verbose=True)
