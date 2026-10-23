import sys
from optimize_1x2_switch_meep import simulate_1x2_switch

st = sys.argv[1] if len(sys.argv) > 1 else "amorphous"
res = simulate_1x2_switch(w1=0.240, w2=0.240, gap=0.035, l_c=3.38, t_pcm=0.150, state=st, resolution=20)
print("\n" + "=" * 60)
print(f"OPTIMIZED 1x2 SWITCH ({st.upper()}, L_c = 3.38 um, gap = 35 nm):")
print(f"  Footprint Area : {res['area_um2']:.3f} um^2 ({res['cell_len']:.2f} um x {res['cell_w']:.2f} um)")
print(f"  Port 1 (Bar)   : {res['p1_bar_pct']:.2f}%")
print(f"  Port 2 (Cross) : {res['p2_cross_pct']:.2f}%")
print(f"  Total Trans    : {res['tot_trans_pct']:.2f}%")
print(f"  Insertion Loss : {res['il_db']:.3f} dB")
print(f"  Crosstalk      : {res['xt_db']:.2f} dB")
print("=" * 60)
