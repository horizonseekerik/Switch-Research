import numpy as np
from scipy.integrate import solve_ivp
import math

lambda_0 = 1.064
k_0 = 2.0 * math.pi / lambda_0

# MPB refractive indices:
n_bare = 2.253
n_am = 2.54110
n_cr = 2.73897

# Evanescent decay constants in cladding:
n_clad = 1.449
gamma_bare = k_0 * math.sqrt(n_bare**2 - n_clad**2)
gamma_am = k_0 * math.sqrt(n_am**2 - n_clad**2)
gamma_cr = k_0 * math.sqrt(n_cr**2 - n_clad**2)

# Coupling at 100nm gap:
# L_c ~ 4.40 um in amorphous state => kappa_am_0 = pi / (2 * 4.40)
kappa_am_0 = math.pi / (2.0 * 4.40)
# In crystalline state, mode is more tightly confined, so field in gap is slightly smaller:
# Mode confinement factor ~ n_core/n_clad ratio or overlap integral
kappa_cr_0 = kappa_am_0 * 0.92  # ~8% reduction due to tighter confinement

def solve_switch(arch='old', state='crystalline', l_c=4.40, l_tip=0.40, gap_0=0.100, y_pitch=0.560, l_bend=1.050, l_ext=0.0):
    # Lengths
    w_wg = 0.240
    
    # Define x domain:
    # Let x=0 be the start of the straight coupling region
    if arch == 'old':
        x_min = -0.350  # inverse taper of WG2 starts at -0.350 at gap_0!
        x_max = l_c + 2*l_tip + l_bend
    elif arch == 'curve_in':
        x_min = -l_bend # WG2 curves in from y_pitch
        x_max = l_c + 2*l_tip + l_bend
    elif arch == 'pure_uniform':
        x_min = -l_bend
        x_max = l_c + l_bend

    def geom_props(x):
        # returns (n_eff1, n_eff2, gap, in_wg2)
        # WG1
        if arch in ['old', 'curve_in']:
            # WG1 has PCM patch from 0 to l_c + 2*l_tip (+ l_ext)
            l_pcm1 = l_c + 2*l_tip + l_ext
            if x < 0:
                n1 = n_bare
            elif 0 <= x < l_tip:
                u = x / l_tip
                n_peak = n_cr if state == 'crystalline' else n_am
                n1 = n_bare + (n_peak - n_bare) * 0.5 * (1.0 - math.cos(math.pi * u))
            elif l_tip <= x <= l_pcm1 - l_tip:
                n1 = n_cr if state == 'crystalline' else n_am
            elif l_pcm1 - l_tip < x <= l_pcm1:
                u = (l_pcm1 - x) / l_tip
                n_peak = n_cr if state == 'crystalline' else n_am
                n1 = n_bare + (n_peak - n_bare) * 0.5 * (1.0 - math.cos(math.pi * u))
            else:
                n1 = n_bare
                
            # WG2
            l_pcm2 = l_c + 2*l_tip
            if x < 0:
                n2 = n_bare
            elif 0 <= x < l_tip:
                u = x / l_tip
                n2 = n_bare + (n_am - n_bare) * 0.5 * (1.0 - math.cos(math.pi * u))
            elif l_tip <= x <= l_pcm2 - l_tip:
                n2 = n_am
            elif l_pcm2 - l_tip < x <= l_pcm2:
                u = (l_pcm2 - x) / l_tip
                n2 = n_bare + (n_am - n_bare) * 0.5 * (1.0 - math.cos(math.pi * u))
            else:
                n2 = n_bare
        else: # pure_uniform
            # WG1 and WG2 uniform across l_c
            if 0 <= x <= l_c:
                n1 = n_cr if state == 'crystalline' else n_am
                n2 = n_am
            else:
                n1 = n_bare
                n2 = n_bare

        # Gap calculation:
        if arch == 'old':
            if x < -0.350:
                gap = 10.0 # not present
                in_wg2 = False
            elif -0.350 <= x < 0:
                gap = gap_0
                in_wg2 = True
            elif 0 <= x <= l_c + 2*l_tip:
                gap = gap_0
                in_wg2 = True
            else:
                u = (x - (l_c + 2*l_tip)) / l_bend
                u = min(max(u, 0.0), 1.0)
                # Cosine S-bend
                dy = (y_pitch - (w_wg + gap_0)) * 0.5 * (1.0 - math.cos(math.pi * u))
                gap = gap_0 + dy
                in_wg2 = True
        elif arch in ['curve_in', 'pure_uniform']:
            x_str_end = (l_c + 2*l_tip) if arch == 'curve_in' else l_c
            if x < 0:
                u = (0.0 - x) / l_bend
                u = min(max(u, 0.0), 1.0)
                dy = (y_pitch - (w_wg + gap_0)) * 0.5 * (1.0 - math.cos(math.pi * u))
                gap = gap_0 + dy
                in_wg2 = True
            elif 0 <= x <= x_str_end:
                gap = gap_0
                in_wg2 = True
            else:
                u = (x - x_str_end) / l_bend
                u = min(max(u, 0.0), 1.0)
                dy = (y_pitch - (w_wg + gap_0)) * 0.5 * (1.0 - math.cos(math.pi * u))
                gap = gap_0 + dy
                in_wg2 = True

        return n1, n2, gap, in_wg2

    # Differential equations for A1, A2:
    # A1 = y[0] + 1j*y[1]
    # A2 = y[2] + 1j*y[3]
    def ode(x, y):
        n1, n2, gap, in_wg2 = geom_props(x)
        if not in_wg2 or gap > 0.400:
            kappa = 0.0
        else:
            gamma = gamma_cr if state == 'crystalline' else gamma_am
            kappa_0 = kappa_cr_0 if state == 'crystalline' else kappa_am_0
            kappa = kappa_0 * math.exp(-gamma * (gap - gap_0))
            
        beta1 = k_0 * n1
        beta2 = k_0 * n2
        
        A1 = y[0] + 1j*y[1]
        A2 = y[2] + 1j*y[3]
        
        dA1 = -1j * beta1 * A1 - 1j * kappa * A2
        dA2 = -1j * beta2 * A2 - 1j * kappa * A1
        
        return [dA1.real, dA1.imag, dA2.real, dA2.imag]

    # Initial condition: all light in WG1 at x_min
    y0 = [1.0, 0.0, 0.0, 0.0]
    sol = solve_ivp(ode, [x_min, x_max], y0, max_step=0.01, rtol=1e-6, atol=1e-8)
    
    A1_end = sol.y[0][-1] + 1j*sol.y[1][-1]
    A2_end = sol.y[2][-1] + 1j*sol.y[3][-1]
    P1 = abs(A1_end)**2
    P2 = abs(A2_end)**2
    
    sig = P2 if state == 'amorphous' else P1
    leak = P1 if state == 'amorphous' else P2
    il_db = -10.0 * math.log10(max(sig, 1e-9))
    xt_db = 10.0 * math.log10(max(leak / max(sig, 1e-9), 1e-9))
    return P1, P2, sig, leak, il_db, xt_db

print("==============================================================================")
print("CMT PREDICTION: OLD GEOMETRY (with bare Si inverse taper at entrance):")
print("==============================================================================")
for state in ['amorphous', 'crystalline']:
    p1, p2, sig, leak, il, xt = solve_switch('old', state, l_c=4.40, l_tip=0.40, l_ext=0.35)
    print(f"State: {state.upper():12s} | P1 (Through)={p1*100:6.2f}% | P2 (Cross)={p2*100:6.2f}% | IL={il:6.3f}dB | XT={xt:6.2f}dB")

print("\n==============================================================================")
print("CMT PREDICTION: CURVE-IN GEOMETRY (WG2 curves in smoothly from y_pitch):")
print("==============================================================================")
for state in ['amorphous', 'crystalline']:
    p1, p2, sig, leak, il, xt = solve_switch('curve_in', state, l_c=4.40, l_tip=0.40, l_ext=0.0)
    print(f"State: {state.upper():12s} | P1 (Through)={p1*100:6.2f}% | P2 (Cross)={p2*100:6.2f}% | IL={il:6.3f}dB | XT={xt:6.2f}dB")
