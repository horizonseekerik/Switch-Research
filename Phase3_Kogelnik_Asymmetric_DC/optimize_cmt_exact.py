import sys, os, math
import numpy as np
from scipy.integrate import solve_ivp

# Data from MPB for bare Si and PCM loaded Si (w = 240 nm, h = 220 nm)
# Gap in nm -> kappa in rad/um
# Bare Si kappa(gap) fitted from MPB:
# gap=90nm: 0.2043, 150nm: 0.1563, 210nm: 0.0780, 300nm: 0.0308, 400nm: 0.0091
# Exponential decay: kappa(g) = 0.38 * exp(-g / 0.085) for PCM loaded
# For bare Si: kappa_bare(g) = 0.35 * exp(-g / 0.100)

def kappa_pcm(gap_um):
    # MPB data: gap=0.09 -> 0.387, gap=0.10 -> 0.340, gap=0.11 -> 0.301, gap=0.12 -> 0.265
    return 0.387 * math.exp(-(gap_um - 0.090) / 0.080)

def kappa_bare(gap_um):
    return 0.204 * math.exp(-(gap_um - 0.090) / 0.095)

def simulate_coupler_cmode(lc, gap, t_pcm=0.200, l_in=0.35, l_bend=1.0, dy=0.12):
    # delta_beta for crystalline vs amorphous PCM (from MPB)
    dbeta_pcm = 1.371  # rad/um at t_pcm = 200 nm
    
    # Entrance: x in [0, l_in]. WG2 enters. Average gap is ~0.15 um, bare Si
    # Transition coupling in entrance:
    phi_in = 0.5 * kappa_bare(gap + 0.03) * l_in
    
    # Central coupling: x in [0, lc]. Gap = gap.
    k_pcm = kappa_pcm(gap)
    
    # Exit: x in [0, l_bend]. S-bend separates waveguides from gap to gap + dy
    # Integrate kappa along S-bend:
    u = np.linspace(0, 1, 50)
    g_bend = gap + dy * (3 * u**2 - 2 * u**3)
    k_bend = np.array([kappa_bare(g) for g in g_bend])
    phi_exit = np.trapz(k_bend, u * l_bend)
    
    # Total phase in Amorphous state (Delta_beta = 0):
    phi_am = phi_in + k_pcm * lc + phi_exit
    p2_am = math.sin(phi_am)**2
    p1_am = 1.0 - p2_am
    
    # In Crystalline state:
    # In central section: Delta_beta = dbeta_pcm
    # Effective S = sqrt(k_pcm^2 + (dbeta_pcm/2)^2)
    s_pcm = math.sqrt(k_pcm**2 + (dbeta_pcm / 2.0)**2)
    u_param = dbeta_pcm / (2.0 * k_pcm)
    f_factor = 1.0 / (1.0 + u_param**2)
    
    # State 1 transfer matrix:
    # Entrance bare Si transfer matrix:
    # Bare Si has Delta_beta = 0, so rotation by phi_in
    c_in = math.cos(phi_in)
    s_in = -1j * math.sin(phi_in)
    M_in = np.array([[c_in, s_in], [s_in, c_in]])
    
    # Central section transfer matrix with detuning:
    c_c = math.cos(s_pcm * lc) - 1j * (dbeta_pcm / (2.0 * s_pcm)) * math.sin(s_pcm * lc)
    s_c = -1j * (k_pcm / s_pcm) * math.sin(s_pcm * lc)
    M_c = np.array([[c_c, s_c], [s_c, np.conj(c_c)]])
    
    # Exit bare Si transfer matrix:
    c_ex = math.cos(phi_exit)
    s_ex = -1j * math.sin(phi_exit)
    M_ex = np.array([[c_ex, s_ex], [s_ex, c_ex]])
    
    # Total transfer matrix:
    M_tot = M_ex @ M_c @ M_in
    a_out = M_tot @ np.array([1.0, 0.0])
    
    p1_cr = abs(a_out[0])**2
    p2_cr = abs(a_out[1])**2
    
    return p2_am, p1_cr, phi_in + phi_exit, k_pcm, s_pcm

print("=" * 80)
print("  ANALYTICAL CMT GRID SCAN FOR >98% IN BOTH STATES")
print("=" * 80)
print(f"{'Gap(nm)':<10} {'L_c(um)':<10} {'AM P2(Cross)%':<16} {'CR P1(Bar)%':<15} {'Both >98%?':<15}")
print("-" * 80)

best_score = 0
best_params = None

for gap_nm in range(85, 125, 2):
    gap = gap_nm / 1000.0
    for lc in np.arange(3.0, 4.8, 0.02):
        p2_am, p1_cr, phi_ends, k, s = simulate_coupler_cmode(lc, gap)
        score = min(p2_am, p1_cr)
        if p2_am >= 0.95 and p1_cr >= 0.95:
            match = "YES! PASS!" if (p2_am >= 0.98 and p1_cr >= 0.98) else "Near"
            print(f"{gap_nm:<10} {lc:<10.2f} {p2_am*100:<16.2f} {p1_cr*100:<15.2f} {match:<15}")
        if score > best_score:
            best_score = score
            best_params = (gap_nm, lc, p2_am, p1_cr)

print("-" * 80)
print(f"Optimal Configuration: Gap = {best_params[0]} nm, L_c = {best_params[1]:.2f} um")
print(f"Predicted Performance: AM P2 = {best_params[2]*100:.2f}%, CR P1 = {best_params[3]*100:.2f}%")
print("=" * 80)
