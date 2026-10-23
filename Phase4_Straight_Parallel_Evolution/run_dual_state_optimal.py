"""
Run both states of the optimal 3D directional coupler switch
at gap = 105 nm, L_c = 4.35 um
"""
import sys, os
from test_dc_optimal_3d import run_dc_switch

if __name__ == "__main__":
    lc = float(sys.argv[1]) if len(sys.argv) > 1 else 4.35
    gp = float(sys.argv[2]) if len(sys.argv) > 2 else 0.105

    print(f"\n{'='*75}")
    print(f"RUNNING 3D FULL-WAVE VERIFICATION: L_c = {lc:.2f} um, gap = {gp*1000:.0f} nm")
    print(f"{'='*75}\n")

    res_am = run_dc_switch(state="amorphous", l_c=lc, gap=gp)
    res_cr = run_dc_switch(state="crystalline", l_c=lc, gap=gp)

    print("\n" + "#"*75)
    print(f"FINAL DUAL-STATE VERIFICATION SUMMARY (L_c = {lc:.2f} um, gap = {gp*1000:.0f} nm):")
    print("#"*75)
    print(f"{'Metric':<25} {'Amorphous (State 0)':<25} {'Crystalline (State 1)':<25}")
    print("-" * 75)
    print(f"{'Port 1 (Bar)':<25} {res_am['p1']*100:.2f}%{'':<19} {res_cr['p1']*100:.2f}%")
    print(f"{'Port 2 (Cross)':<25} {res_am['p2']*100:.2f}%{'':<19} {res_cr['p2']*100:.2f}%")
    print(f"{'Total Trans':<25} {res_am['tot']*100:.2f}%{'':<19} {res_cr['tot']*100:.2f}%")
    print(f"{'Insertion Loss':<25} {res_am['il']:.3f} dB{'':<18} {res_cr['il']:.3f} dB")
    print(f"{'Crosstalk':<25} {res_am['xt']:.2f} dB{'':<18} {res_cr['xt']:.2f} dB")
    print(f"{'Active Footprint':<25} {res_am['area']:.3f} um^2{'':<18} {res_cr['area']:.3f} um^2")
    fabric_area = 4096 * res_am['area'] * 1e-6
    print(f"{'4096-Switch Fabric':<25} {fabric_area:.2f} mm^2{'':<18} {fabric_area:.2f} mm^2")
    print("#"*75 + "\n")
