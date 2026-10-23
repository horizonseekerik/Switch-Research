import meep as mp
from meep import mpb

lambda_0 = 1.064
fcen = 1.0 / lambda_0
w = 0.240
h = 0.220
gap = 0.10551

# MPB simulation for bare silicon coupled waveguides
geometry_lattice = mp.Lattice(size=mp.Vector3(0, 2.0, 1.5))
default_material = mp.Medium(index=1.449)

geom_bare = [
    mp.Block(center=mp.Vector3(0, (w+gap)/2.0, 0), size=mp.Vector3(mp.inf, w, h), material=mp.Medium(index=3.48)),
    mp.Block(center=mp.Vector3(0, -(w+gap)/2.0, 0), size=mp.Vector3(mp.inf, w, h), material=mp.Medium(index=3.48))
]

ms = mpb.ModeSolver(
    geometry_lattice=geometry_lattice,
    geometry=geom_bare,
    resolution=32,
    default_material=default_material
)

# Find even and odd modes at fcen
k_guess = fcen * 2.0
k_even = ms.find_k(mp.NO_PARITY, fcen, 1, 1, mp.Vector3(1, 0, 0), 1e-4, k_guess, 0.1, 4.0)
k_odd = ms.find_k(mp.NO_PARITY, fcen, 2, 2, mp.Vector3(1, 0, 0), 1e-4, k_guess, 0.1, 4.0)

ke = k_even[0]
ko = k_odd[0]
delta_k = abs(ke - ko)
Lc_bare = 0.5 / delta_k

print(f"Bare Si: k_even = {ke:.6f}, k_odd = {ko:.6f}, Delta_k = {delta_k:.6f} um^-1, Lc = {Lc_bare:.6f} um")
