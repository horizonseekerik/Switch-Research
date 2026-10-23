import meep as mp
p = mp.Prism(vertices=[mp.Vector3(1, 1), mp.Vector3(2, 1), mp.Vector3(2, 2)], height=1, center=mp.Vector3(10, 10, 0))
print("center:", p.center)
for v in p.vertices:
    print("v:", v)
sim = mp.Simulation(cell_size=mp.Vector3(20, 20, 2), resolution=1, geometry=[p])
sim.init_sim()
print("Epsilon at (1.5, 1.2, 0):", sim.get_epsilon_point(mp.Vector3(1.5, 1.2, 0)))
print("Epsilon at (11.5, 11.2, 0):", sim.get_epsilon_point(mp.Vector3(11.5, 11.2, 0)))
