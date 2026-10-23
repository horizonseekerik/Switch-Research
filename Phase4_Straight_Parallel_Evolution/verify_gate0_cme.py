"""
PROJECT JANUS: GATE 0 COUPLED-MODE EQUATIONS (CME) SOLVER & CO-DESIGN
====================================================================
Physics Principle:
  Solves non-autonomous coupled-mode equations (CME) along the longitudinal 
  coordinate x in [0, L_total] incorporating:
    1. Input fanout coupling integral
    2. Central interaction zone with candidate modulation profile g(x)
    3. Output Hermite S-bend fanout coupling integral
    4. Crystalline detuning profile Delta_beta_1(x) and kappa_1(x)
    5. Amorphous phase-matching Delta_beta_0(x) = 0 and kappa_0(x)

Outputs all values to >= 10 decimal places.
Proves whether a profile achieves theoretical extinction > 99.0000000000%
prior to running expensive 3D Maxwell FDTD.
"""

import math
import numpy as np
from scipy.integrate import solve_ivp

# --- Empirical MPB Database Constants at lambda_0 = 1064 nm (w = 240 nm, h = 220 nm) ---
GAP_REF = 0.09814            # 98.14 nm
KAPPA_0_REF = 0.3654531482   # um^-1 (L_c0 = 4.298177 um)
GAMMA_0 = 12.8532104521      # um^-1 (decay constant d ln(kappa)/dg)

KAPPA_1_REF = 0.3337034185   # um^-1
GAMMA_1 = 13.1245892105      # um^-1
DELTA_BETA_1 = 1.2676834192  # um^-1 (2*pi*fcen * (n_cr - n_am))

def get_gap(x, l_in, l_c, l_bend, g0, delta_g, y_pitch=0.450, w_wg=0.240):
    """Calculates local edge-to-edge gap g(x) in um across all three regions."""
    if x < l_in:
        # Input bare-Si inverse taper / straight access: parallel or wide
        return g0
    elif x <= l_in + l_c:
        # Central interaction zone: smooth profile g(x)
        u = (x - l_in) / l_c
        # Gentle symmetric cosine elevation
        return g0 + delta_g * math.sin(math.pi * u)**2
    else:
        # Output C^2 Hermite S-bend fanout
        u = (x - (l_in + l_c)) / l_bend
        u = min(max(u, 0.0), 1.0)
        # S-bend curve: y_shift(u) = dy * (10*u^3 - 15*u^4 + 6*u^5)
        dy = y_pitch - (w_wg + g0)
        y_bend = dy * (10.0 * u**3 - 15.0 * u**4 + 6.0 * u**5)
        return g0 + y_bend

def cme_amorphous(x, y, l_in, l_c, l_bend, g0, delta_g):
    # y = [Re(a1), Im(a1), Re(a2), Im(a2)]
    g = get_gap(x, l_in, l_c, l_bend, g0, delta_g)
    kappa = KAPPA_0_REF * math.exp(-GAMMA_0 * (g - GAP_REF))
    # da1/dx = -j kappa a2
    # da2/dx = -j kappa a1
    # a1 = u1 + j v1 -> du1/dx = kappa v2, dv1/dx = -kappa u2
    # a2 = u2 + j v2 -> du2/dx = kappa v1, dv2/dx = -kappa u1
    u1, v1, u2, v2 = y
    du1 = kappa * v2
    dv1 = -kappa * u2
    du2 = kappa * v1
    dv2 = -kappa * u1
    return [du1, dv1, du2, dv2]

def cme_crystalline(x, y, l_in, l_c, l_bend, g0, delta_g):
    # y = [Re(a1), Im(a1), Re(a2), Im(a2)]
    g = get_gap(x, l_in, l_c, l_bend, g0, delta_g)
    kappa = KAPPA_1_REF * math.exp(-GAMMA_1 * (g - GAP_REF))
    
    # Delta_beta is active only where PCM is on WG1: x in [l_in, l_in + l_c + l_hold]
    # For conservative analysis, hold length l_hold = 0.35 um into fanout
    l_hold = 0.350
    if l_in <= x <= (l_in + l_c + l_hold):
        dbeta = DELTA_BETA_1
    else:
        dbeta = 0.0  # bare silicon on both waveguides
        
    delta_half = dbeta / 2.0
    u1, v1, u2, v2 = y
    # da1/dx = -j delta_half a1 - j kappa a2
    # da2/dx = +j delta_half a2 - j kappa a1
    du1 = delta_half * v1 + kappa * v2
    dv1 = -delta_half * u1 - kappa * u2
    du2 = -delta_half * v2 + kappa * v1
    dv2 = delta_half * u2 - kappa * u1
    return [du1, dv1, du2, dv2]

def evaluate_profile(l_c=4.425, g0=0.09814, delta_g=0.000, l_in=0.350, l_bend=1.000):
    x_span = [0.0, l_in + l_c + l_bend]
    y0 = [1.0, 0.0, 0.0, 0.0]  # Launch 100% into Waveguide 1
    
    # Solve Amorphous (State 0)
    sol_am = solve_ivp(cme_amorphous, x_span, y0, args=(l_in, l_c, l_bend, g0, delta_g),
                       rtol=1e-11, atol=1e-13, max_step=0.01)
    u1_0, v1_0, u2_0, v2_0 = sol_am.y[:, -1]
    p1_0 = u1_0**2 + v1_0**2
    p2_0 = u2_0**2 + v2_0**2
    
    # Solve Crystalline (State 1)
    sol_cr = solve_ivp(cme_crystalline, x_span, y0, args=(l_in, l_c, l_bend, g0, delta_g),
                       rtol=1e-11, atol=1e-13, max_step=0.01)
    u1_1, v1_1, u2_1, v2_1 = sol_cr.y[:, -1]
    p1_1 = u1_1**2 + v1_1**2
    p2_1 = u2_1**2 + v2_1**2
    
    return {
        "p2_am": p2_0 * 100.0, "p1_am": p1_0 * 100.0,
        "p1_cr": p1_1 * 100.0, "p2_cr": p2_1 * 100.0,
        "xt_am": 10.0 * math.log10(max(p1_0 / p2_0, 1e-15)),
        "xt_cr": 10.0 * math.log10(max(p2_1 / p1_1, 1e-15))
    }

if __name__ == "__main__":
    print("=" * 115)
    print("GATE 0: COUPLED-MODE EQUATIONS (CME) VERIFICATION & PROFILE OPTIMIZATION")
    print("=" * 115)
    print(f"{'L_c (um)':<10} {'g0 (nm)':<10} {'dg (nm)':<10} {'State 0 P2 %':<16} {'State 0 P1 %':<16} {'State 1 P1 %':<16} {'State 1 P2 %':<16} {'State 1 XT (dB)':<16}")
    print("-" * 115)
    
    # 1. Baseline Uniform at L_c = 4.425 um
    res_base = evaluate_profile(l_c=4.425, g0=0.09814, delta_g=0.000)
    print(f"{4.425:<10.3f} {98.14:<10.2f} {0.0:<10.2f} {res_base['p2_am']:<16.10f} {res_base['p1_am']:<16.10f} {res_base['p1_cr']:<16.10f} {res_base['p2_cr']:<16.10f} {res_base['xt_cr']:<16.10f}")
    
    # 2. Sweep delta_g (elevation) at L_c = 4.425 um to find exact State 0 lock
    best_entry = None
    min_error = 1e9
    for dg_nm in np.linspace(0.0, 14.0, 15):
        dg = dg_nm * 1e-3
        r = evaluate_profile(l_c=4.425, g0=0.09814, delta_g=dg)
        err = (100.0 - r['p2_am'])**2 + (100.0 - r['p1_cr'])**2
        if err < min_error:
            min_error = err
            best_entry = (dg_nm, r)
        print(f"{4.425:<10.3f} {98.14:<10.2f} {dg_nm:<10.2f} {r['p2_am']:<16.10f} {r['p1_am']:<16.10f} {r['p1_cr']:<16.10f} {r['p2_cr']:<16.10f} {r['xt_cr']:<16.10f}")
        
    print("=" * 115)
    dg_opt, r_opt = best_entry
    print(f"OPTIMAL GATE 0 POINT: delta_g = {dg_opt:.2f} nm")
    print(f"  State 0 (Amorphous Cross P2) = {r_opt['p2_am']:.10f}%, Bar Leakage = {r_opt['p1_am']:.10f}%")
    print(f"  State 1 (Crystalline Bar P1) = {r_opt['p1_cr']:.10f}%, Cross Leakage = {r_opt['p2_cr']:.10f}%")
    print(f"  State 1 Optical Crosstalk    = {r_opt['xt_cr']:.10f} dB")
    print("=" * 115)
