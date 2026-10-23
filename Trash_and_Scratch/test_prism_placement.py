import meep as mp

v = [mp.Vector3(1, 1), mp.Vector3(2, 1), mp.Vector3(2, 2), mp.Vector3(1, 2)]
p1 = mp.Prism(vertices=v, height=1)
p2 = mp.Prism(vertices=v, height=1, center=mp.Vector3(5, 5, 0))

cell = mp.Vector3(10, 10, 2)
sim1 = mp.Simulation(cell_size=cell, geometry=[p1], resolution=10)
sim1.init_sim()
print("p1 at (1.5, 1.5, 0): eps =", sim1.get_epsilon_point(mp.Vector3(1.5, 1.5, 0)))

sim2 = mp.Simulation(cell_size=cell, geometry=[p2], resolution=10)
sim2.init_sim()
print("p2 at (1.5, 1.5, 0): eps =", sim2.get_epsilon_point(mp.Vector3(1.5, 1.5, 0)))
print("p2 at (5.0, 5.0, 0): eps =", sim2.get_epsilon_point(mp.Vector3(5.0, 5.0, 0)))
print("p2 at (6.5, 6.5, 0): eps =", sim2.get_epsilon_point(mp.Vector3(6.5, 6.5, 0)))
