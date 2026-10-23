import math
import numpy as np
import meep as mp

def test_bend_lengths():
    w_wg = 0.240
    h_wg = 0.220
    dy = 0.200
    lambda_0 = 1.064
    fcen = 1.0 / lambda_0
    df = 0.1 * fcen
    clad = mp.Medium(index=1.449)
    core = mp.Medium(index=3.480)

    for l_bend in [1.50, 2.00, 2.50]:
        l_in, l_out = 0.5, 0.5
        sx = l_in + l_bend + l_out + 1.0
        sy = dy + 1.4
        sz = h_wg + 1.0
        cell = mp.Vector3(sx, sy, sz)
        pml = [mp.PML(0.4)]

        x_start = -l_bend / 2.0
        x_end = l_bend / 2.0
        y_start = 0.0
        y_end = dy

        # Cosine S-bend
        def make_cosine(x0, x1, y0, y1, w, n_pts=50):
            dx = x1 - x0
            dy_ = y1 - y0
            u = np.linspace(0, 1, n_pts)
            x = x0 + u * dx
            y = y0 + dy_ * (0.5 * (1.0 - np.cos(np.pi * u)))
            top = [mp.Vector3(xi, yi + w / 2.0, 0) for xi, yi in zip(x, y)]
            bot = [mp.Vector3(xi, yi - w / 2.0, 0) for xi, yi in reversed(list(zip(x, y)))]
            return top + bot

        src = [mp.EigenModeSource(src=mp.GaussianSource(fcen, fwidth=df),
                                  center=mp.Vector3(-sx/2 + 0.5, 0, h_wg/2),
                                  size=mp.Vector3(0, w_wg*3, h_wg*3),
                                  eig_band=1, direction=mp.X)]

        # Ref straight
        sim_ref = mp.Simulation(cell_size=cell, boundary_layers=pml,
                                geometry=[mp.Block(center=mp.Vector3(0,0,h_wg/2), size=mp.Vector3(mp.inf, w_wg, h_wg), material=core)],
                                sources=src, resolution=20, default_material=clad)
        f_ref = sim_ref.add_flux(fcen, 0, 1, mp.FluxRegion(center=mp.Vector3(sx/2 - 0.6, 0, h_wg/2), size=mp.Vector3(0, 0.4, 0.4)))
        sim_ref.run(until_after_sources=mp.stop_when_fields_decayed(20, mp.Hz, mp.Vector3(sx/2 - 0.6, 0, h_wg/2), 1e-3))
        p_ref = max(mp.get_fluxes(f_ref)[0], 1e-12)

        # Bend
        geom = [
            mp.Block(center=mp.Vector3((-sx/2 + x_start)/2, y_start, h_wg/2), size=mp.Vector3(x_start - (-sx/2), w_wg, h_wg), material=core),
            mp.Prism(vertices=make_cosine(x_start, x_end, y_start, y_end, w_wg), height=h_wg, axis=mp.Vector3(0,0,1), material=core),
            mp.Block(center=mp.Vector3((x_end + sx/2)/2, y_end, h_wg/2), size=mp.Vector3(sx/2 - x_end, w_wg, h_wg), material=core)
        ]
        sim_b = mp.Simulation(cell_size=cell, boundary_layers=pml, geometry=geom, sources=src, resolution=20, default_material=clad)
        f_b = sim_b.add_flux(fcen, 0, 1, mp.FluxRegion(center=mp.Vector3(sx/2 - 0.6, y_end, h_wg/2), size=mp.Vector3(0, 0.4, 0.4)))
        sim_b.run(until_after_sources=mp.stop_when_fields_decayed(20, mp.Hz, mp.Vector3(sx/2 - 0.6, y_end, h_wg/2), 1e-3))
        p_b = max(mp.get_fluxes(f_b)[0], 1e-12)

        t = p_b / p_ref
        il = -10 * math.log10(max(t, 1e-6))
        print(f"[RESULT] L_bend = {l_bend:.2f} um, dy = {dy:.3f} um -> Trans = {t*100:.2f}%, IL = {il:.3f} dB")

if __name__ == '__main__':
    test_bend_lengths()
