import sys, os, json
from run_optimized_1x2_switch_meep import run_1x2_switch_3d

if __name__ == '__main__':
    lc = float(sys.argv[1]) if len(sys.argv) > 1 else 4.65
    gap = float(sys.argv[2]) if len(sys.argv) > 2 else 0.040

    print("=" * 75)
    print("PROJECT JANUS: FULL 3D MEEP FDTD VERIFICATION OF BOTH SWITCH STATES")
    print("=" * 75)

    print("\n--- 1. RUNNING AMORPHOUS STATE (Target: Port 2 Cross) ---")
    res_am = run_1x2_switch_3d(state="amorphous", l_c=lc, gap=gap)

    print("\n--- 2. RUNNING CRYSTALLINE STATE (Target: Port 1 Bar) ---")
    res_cr = run_1x2_switch_3d(state="crystalline", l_c=lc, gap=gap)

    summary = {
        "device_architecture": "1x2 Non-Volatile Directional Coupler Switch (Zero-Bend Bar WG)",
        "simulation_mode": "3D Full-Wave Vector Maxwell FDTD (MIT MEEP)",
        "wavelength_nm": 1064.0,
        "parameters": {
            "w_wg_nm": 240,
            "h_wg_nm": 220,
            "t_pcm_nm": 150,
            "gap_nm": gap * 1000,
            "l_c_um": lc,
            "l_bend_um": 1.50,
            "footprint_length_um": res_am["cell_len"],
            "footprint_width_um": res_am["cell_w"],
            "footprint_area_um2": res_am["area_um2"]
        },
        "amorphous_state": res_am,
        "crystalline_state": res_cr
    }

    out_file = "final_1x2_optimized_results.json"
    with open(out_file, "w") as f:
        json.dump(summary, f, indent=2)

    print("\n" + "=" * 75)
    print("FINAL 1x2 OPTICAL SWITCH 3D FDTD PERFORMANCE SUMMARY:")
    print("=" * 75)
    print(f"Active Footprint Area : {res_am['area_um2']:.3f} um^2 ({res_am['cell_len']:.2f} um x {res_am['cell_w']:.2f} um)")
    print(f"\nAMORPHOUS STATE (Target: Port 2 Cross):")
    print(f"  Port 1 (Bar)   : {res_am['p1_bar_pct']:.2f}%")
    print(f"  Port 2 (Cross) : {res_am['p2_cross_pct']:.2f}% (Signal)")
    print(f"  Total Trans    : {res_am['tot_trans_pct']:.2f}%")
    print(f"  Insertion Loss : {res_am['il_db']:.3f} dB (Target: <= 0.08 dB)")
    print(f"  Crosstalk      : {res_am['xt_db']:.2f} dB (Target: <= -20.0 dB)")
    print(f"\nCRYSTALLINE STATE (Target: Port 1 Bar):")
    print(f"  Port 1 (Bar)   : {res_cr['p1_bar_pct']:.2f}% (Signal)")
    print(f"  Port 2 (Cross) : {res_cr['p2_cross_pct']:.2f}%")
    print(f"  Total Trans    : {res_cr['tot_trans_pct']:.2f}%")
    print(f"  Insertion Loss : {res_cr['il_db']:.3f} dB (Target: <= 0.08 dB)")
    print(f"  Crosstalk      : {res_cr['xt_db']:.2f} dB (Target: <= -20.0 dB)")
    print("=" * 75)
