import meep as mp
import numpy as np

# Define a 1x1 square at x in [2, 3], y in [1, 2]
pts = [mp.Vector3(2, 1, 0), mp.Vector3(3, 1, 0), mp.Vector3(3, 2, 0), mp.Vector3(2, 2, 0)]
centroid_x = np.mean([p.x for p in pts])
centroid_y = np.mean([p.y for p in pts])
print(f"Original vertices centroid: ({centroid_x}, {centroid_y})")

# Let's create a simulation and check where epsilon is placed!
cell = mp.Vector3(6, 6, 2)
p1 = mp.Prism(vertices=pts, height=0.5, axis=mp.Vector3(0, 0, 1), material=mp.Medium(index=2.0))
sim1 = mp.Simulation(cell_size=cell, geometry=[p1], resolution=10)
sim1.init_sim()
eps_at_orig = sim1.get_epsilon_point(mp.Vector3(2.5, 1.5, 0.25))
eps_at_zero = sim1.get_epsilon_point(mp.Vector3(0.0, 0.0, 0.25))
print(f"Case 1 (omitted center):")
print(f"  eps at (2.5, 1.5, 0.25) [original location]: {eps_at_orig}")
print(f"  eps at (0.0, 0.0, 0.25) [shifted to center]: {eps_at_zero}")

# Case 2: pass center = mp.Vector3(centroid_x, centroid_y, 0.25)
p2 = mp.Prism(vertices=pts, height=0.5, axis=mp.Vector3(0, 0, 1), center=mp.Vector3(centroid_x, centroid_y, 0.25), material=mp.Medium(index=2.0))
sim2 = mp.Simulation(cell_size=cell, geometry=[p2], resolution=10)
sim2.init_sim()
eps2_at_orig = sim2.get_epsilon_point(mp.Vector3(2.5, 1.5, 0.25))
eps2_at_zero = sim2.get_epsilon_point(mp.Vector3(0.0, 0.0, 0.25))
print(f"\nCase 2 (center = (centroid_x, centroid_y, z_center)):")
print(f"  eps at (2.5, 1.5, 0.25) [original location]: {eps2_at_orig}")
print(f"  eps at (0.0, 0.0, 0.25) [shifted to center]: {eps2_at_zero}")

# Case 3: pass center = mp.Vector3(0, 0, 0.25)
p3 = mp.Prism(vertices=pts, height=0.5, axis=mp.Vector3(0, 0, 1), center=mp.Vector3(0, 0, 0.25), material=mp.Medium(index=2.0))
sim3 = mp.Simulation(cell_size=cell, geometry=[p3], resolution=10)
sim3.init_sim()
eps3_at_orig = sim3.get_epsilon_point(mp.Vector3(2.5, 1.5, 0.25))
eps3_at_zero = sim3.get_epsilon_point(mp.Vector3(0.0, 0.0, 0.25))
print(f"\nCase 3 (center = (0, 0, z_center)):")
print(f"  eps at (2.5, 1.5, 0.25) [original location]: {eps3_at_orig}")
print(f"  eps at (0.0, 0.0, 0.25) [shifted to center]: {eps3_at_zero}")
