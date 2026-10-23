"""
PROJECT JANUS: DUAL-STATE 3D FDTD VERIFICATION
===============================================
Rapidly verifies both Amorphous and Crystalline states at a given
(L_c, gap) and reports against all design targets.

Usage:
  python verify_both_states_3d.py                         # defaults L_c=3.55, gap=95nm
  python verify_both_states_3d.py --lc 3.55 --gap 0.095
  python verify_both_states_3d.py --lc 3.60 --gap 0.090
"""
import argparse
from master_directional_coupler_3d import run_both_states

TARGETS = {
    "signal_port_pct": (">= 98.00%", lambda x: x >= 98.0),
    "leakage_port_pct": ("<= 1.00%",  lambda x: x <= 1.0),
    "il_db":            ("<= 0.088 dB", lambda x: x <= 0.088),
    "xt_db":            ("<= -20.0 dB", lambda x: x <= -20.0),
    "area_um2":         ("< 5.0 um^2",  lambda x: x < 5.0),
}

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--lc",         type=float, default=3.55)
    parser.add_argument("--gap",        type=float, default=0.095)
    parser.add_argument("--resolution", type=int,   default=20)
    args = parser.parse_args()

    r_am, r_cr = run_both_states(l_c=args.lc, gap=args.gap, resolution=args.resolution)

    print("\n" + "="*70)
    print("SPEC COMPLIANCE CHECK")
    print("="*70)
    all_ok = True
    for key, (spec, fn) in TARGETS.items():
        va = r_am[key]
        vc = r_cr[key]
        pa = "✓ PASS" if fn(va) else "✗ FAIL"
        pc = "✓ PASS" if fn(vc) else "✗ FAIL"
        if not (fn(va) and fn(vc)):
            all_ok = False
        print(f"  {key:<22} Target {spec:<15} | Am: {va:>8.3f} {pa} | Cr: {vc:>8.3f} {pc}")

    print("="*70)
    if all_ok:
        print("  *** ALL DESIGN TARGETS MET FOR BOTH STATES! ***")
    else:
        print("  Some targets not yet met — continue parameter tuning.")
    print("="*70 + "\n")
