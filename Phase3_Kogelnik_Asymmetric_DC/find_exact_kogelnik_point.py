"""
FIND EXACT KOGELNIK GAP WHERE u = Delta_beta / (2*kappa) = sqrt(3) = 1.7320508
=============================================================================
Calculates the exact gap between 89 nm and 98 nm where:
  L_c(amorphous) == L_null(crystalline)
"""

import math
import numpy as np

# From 3D FDTD measurements:
# At gap = 89 nm: L_c(am) = 3.00 um -> kappa = pi / (2 * 3.00) = 0.52360 rad/um
# At gap = 98 nm: L_c(am) = 4.28 um -> kappa = pi / (2 * 4.28) = 0.36701 rad/um

gap_1 = 89.0
kappa_1 = math.pi / (2.0 * 3.00)
gap_2 = 98.0
kappa_2 = math.pi / (2.0 * 4.28)

# Evanescent coupling decay: kappa(g) = kappa_0 * exp(-gamma * g)
gamma = math.log(kappa_1 / kappa_2) / (gap_2 - gap_1)
kappa_0 = kappa_1 * math.exp(gamma * gap_1)

print(f"Evanescent decay rate gamma = {gamma:.6f} nm^-1")
print(f"Target u = sqrt(3) = {math.sqrt(3):.10f}")

# At gap = 98 nm, L_null(cr) was 4.28 um -> S(98nm) = pi / 4.28 = 0.7340 rad/um
# S^2 = kappa^2 + (Delta_beta/2)^2 -> (Delta_beta/2)^2 = S^2 - kappa^2
S_2 = math.pi / 4.28
kappa_2_val = kappa_2
delta_beta_half = math.sqrt(S_2**2 - kappa_2_val**2)
delta_beta = 2.0 * delta_beta_half
print(f"Measured Delta_beta = {delta_beta:.6f} rad/um")

# For u = sqrt(3), we require kappa_target = Delta_beta / (2 * sqrt(3))
kappa_target = delta_beta / (2.0 * math.sqrt(3))
L_c_target = math.pi / (2.0 * kappa_target)
gap_target = (math.log(kappa_0 / kappa_target)) / gamma

print(f"\nEXACT ANALYTICAL SOLUTION:")
print(f"  Target kappa       = {kappa_target:.10f} rad/um")
print(f"  Target Gap         = {gap_target:.10f} nm")
print(f"  Target Length L_c  = {L_c_target:.10f} um")
print(f"  At this point: L_c(am) = {math.pi/(2*kappa_target):.10f} um, L_null(cr) = {math.pi/(2*kappa_target):.10f} um!")
