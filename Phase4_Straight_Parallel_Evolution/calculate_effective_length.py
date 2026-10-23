import math
import numpy as np

# MPB results at Gap = 103.7165 nm, w = 240 nm:
# lc_0 = 4.722 um -> kappa_0 = pi / (2 * lc_0) = pi / 9.444 um^-1
lc_0 = 4.7220000000
kappa_0 = math.pi / (2.0 * lc_0) # um^-1

# Evanescent decay constant in SiO2 cladding:
# gamma = sqrt(beta^2 - k0^2 * n_clad^2)
lambda_0 = 1.064
fcen = 1.0 / lambda_0
n_eff = 2.5472041900
n_clad = 1.449
gamma = 2.0 * math.pi * fcen * math.sqrt(n_eff**2 - n_clad**2)

print(f"kappa_0 = {kappa_0:.10f} um^-1")
print(f"gamma   = {gamma:.10f} um^-1")

# In S-bend fanout: y(x) = y0 + dy * (10*u^3 - 15*u^4 + 6*u^5)
# Extra gap Delta_g(x) = dy * (10*u^3 - 15*u^4 + 6*u^5)
# kappa(x) = kappa_0 * exp(-gamma * Delta_g(x))
# We integrate kappa(x) over the S-bend:
l_bend = 1.0000000000
y_pitch = 0.4500000000
w = 0.2400000000
gap = 0.1037165000
dy = y_pitch - (w + gap)

u = np.linspace(0, 1, 1000)
delta_g = dy * (10.0 * u**3 - 15.0 * u**4 + 6.0 * u**5)
kappa_bend = kappa_0 * np.exp(-gamma * delta_g)
integrated_phase_bend = np.trapz(kappa_bend, u * l_bend)
l_bend_eff = integrated_phase_bend / kappa_0

print(f"dy                 = {dy*1e3:.10f} nm")
print(f"l_bend_eff         = {l_bend_eff:.10f} um")
print(f"Total phase in bend = {integrated_phase_bend:.10f} rad (vs required pi/2 = {math.pi/2:.10f})")

# Entrance inverse taper on WG2:
# WG2 core tapers from 20 to 240 nm over l_in = 0.350 um
# Coupling is only significant over the second half where w > 160 nm
# Approximate effective length of entrance:
l_in = 0.3500000000
l_in_eff = 0.0800000000 # um

l_c_phys = lc_0 - l_bend_eff - l_in_eff
print(f"Required physical parallel length L_c_phys = {l_c_phys:.10f} um")
