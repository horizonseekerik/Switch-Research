import sys, math
import numpy as np

def simulate_segmented_switch():
    # Section 1: Length L1 = 2.91 um, gap = 89 nm
    # State 0 is perfectly phase-matched: Delta_beta_1,am = 0
    # State 1 is detuned: Delta_beta_1,cr = 0.978 rad/um
    k1 = 0.5400000000
    L1 = math.pi / (2.0 * k1)  # 2.9088820867 um
    
    # In State 0:
    phi1_am = k1 * L1  # pi / 2 = 1.5707963268 rad
    # Section 1 transfer matrix for State 0:
    # WG1 -> 0, WG2 -> 1
    # State 1:
    dbeta1 = 0.9780000000
    s1 = math.sqrt(k1**2 + (dbeta1 / 2.0)**2)
    c1 = math.cos(s1 * L1) - 1j * (dbeta1 / (2.0 * s1)) * math.sin(s1 * L1)
    s1_val = -1j * (k1 / s1) * math.sin(s1 * L1)
    M1_cr = np.array([[c1, s1_val], [s1_val, np.conj(c1)]])
    state1_after_sec1 = M1_cr @ np.array([1.0, 0.0])
    
    # Section 2: Length L2.
    # In Section 2, Waveguide 2 has core width offset such that:
    # State 1 is phase-matched: Delta_beta_2,cr = 0
    # State 0 is strongly detuned: Delta_beta_2,am = -0.978 rad/um
    # We want Section 2 to rotate State 1 back to [1.0, 0.0]
    # Let k2 be the coupling in Section 2.
    k2 = 0.4000000000
    dbeta2_am = 1.8500000000  # Strong detuning in State 0
    dbeta2_cr = 0.0000000000  # Exact phase match in State 1
    
    print("=" * 100)
    print("  APPROACH 1: SEGMENTED DUAL-SECTION SWITCH (10 DECIMAL PLACES ACCURACY)")
    print("=" * 100)
    print(f"Section 1: Length L1 = {L1:.10f} um, kappa_1 = {k1:.10f} rad/um")
    print(f"State 0 Power after Sec 1: P2 = {(math.sin(phi1_am)**2)*100:.10f}%, P1 = {(math.cos(phi1_am)**2)*100:.10f}%")
    print(f"State 1 Amplitude after Sec 1: |A1|^2 = {abs(state1_after_sec1[0])**2*100:.10f}%, |A2|^2 = {abs(state1_after_sec1[1])**2*100:.10f}%")
    print("-" * 100)
    print(f"{'L2 (um)':<14} {'State 0 P2%':<22} {'State 0 P1%':<22} {'State 1 P1%':<22} {'State 1 P2%':<22} {'Both >98%?'}")
    print("-" * 100)

    best_res = None
    for L2 in np.arange(0.50, 2.50, 0.05):
        # Section 2 for State 0 (Detuned):
        s2_am = math.sqrt(k2**2 + (dbeta2_am / 2.0)**2)
        c2_am = math.cos(s2_am * L2) - 1j * (dbeta2_am / (2.0 * s2_am)) * math.sin(s2_am * L2)
        s2_am_val = -1j * (k2 / s2_am) * math.sin(s2_am * L2)
        M2_am = np.array([[c2_am, s2_am_val], [s2_am_val, np.conj(c2_am)]])
        
        # State 0 enters Section 2 with [0.0, 1.0]:
        out_am = M2_am @ np.array([0.0, 1.0])
        p1_am = abs(out_am[0])**2
        p2_am = abs(out_am[1])**2
        
        # Section 2 for State 1 (Phase-matched Delta_beta = 0):
        # M2_cr is pure rotation:
        phi2_cr = k2 * L2
        c2_cr = math.cos(phi2_cr)
        s2_cr_val = -1j * math.sin(phi2_cr)
        M2_cr = np.array([[c2_cr, s2_cr_val], [s2_cr_val, c2_cr]])
        
        out_cr = M2_cr @ state1_after_sec1
        p1_cr = abs(out_cr[0])**2
        p2_cr = abs(out_cr[1])**2
        
        pass_str = "YES! BOTH >98%!" if (p2_am >= 0.98 and p1_cr >= 0.98) else "No"
        if p2_am >= 0.95 or p1_cr >= 0.95:
            print(f"{L2:<14.10f} {p2_am*100:<22.10f} {p1_am*100:<22.10f} {p1_cr*100:<22.10f} {p2_cr*100:<22.10f} {pass_str}")
        
        if p2_am >= 0.98 and p1_cr >= 0.98:
            if best_res is None or (p2_am + p1_cr) > (best_res[1] + best_res[2]):
                best_res = (L2, p2_am, p1_cr, p1_am, p2_cr)

    print("=" * 100)
    if best_res:
        print(f"Optimal Section 2 Length: L2 = {best_res[0]:.10f} um")
        print(f"Total Cell Length L_total = {L1 + best_res[0]:.10f} um")
        print(f"State 0 (Amorphous   -> Port 2 Cross): P2 = {best_res[1]*100:.10f}%, P1 (Leak) = {best_res[3]*100:.10f}%")
        print(f"State 1 (Crystalline -> Port 1 Bar)  : P1 = {best_res[2]*100:.10f}%, P2 (Leak) = {best_res[4]*100:.10f}%")
    print("=" * 100)

if __name__ == '__main__':
    simulate_segmented_switch()
