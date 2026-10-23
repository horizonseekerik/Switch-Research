import sys, os, math
import meep as mp
from meep import mpb

def check_janus_baseline():
    lambda_0 = 1.064
    fcen = 1.0 / lambda_0
    w = 0.280
    h = 0.220
    gap = 0.080
    t_pcm = 0.020 # 20 nm!

    clad = mp.Medium(index=1.449)
    core = mp.Medium(index=3.480)
    pcm_am = mp.Medium(index=2.70)
    pcm_cr = mp.Medium(index=3.30)

    # 1. Single waveguide with 20nm PCM
    sy, sz = 3.0, 3.0
    cell = mp.Vector3(0, sy, sz)
    geom_am = [
        mp.Block(center=mp.Vector3(0, 0, -sz/4.0), size=mp.Vector3(mp.inf, mp.inf, sz/2.0), material=clad),
        mp.Block(center=mp.Vector3(0, 0, h/2.0), size=mp.Vector3(mp.inf, w, h), material=core),
        mp.Block(center=mp.Vector3(0, 0, h + t_pcm/2.0), size=mp.Vector3(mp.inf, w, t_pcm), material=pcm_am),
    ]
    ms = mpb.ModeSolver(geometry=geom_am, geometry_lattice=mp.Lattice(size=cell),
                        k_points=[mp.Vector3(2.5*fcen, 0, 0)], resolution=40,
                        default_material=clad, num_bands=2)
    old_stdout = sys.stdout
    sys.stdout = open(os.devnull, 'w')
    k_res_am = ms.find_k(mp.NO_PARITY, fcen, 1, 1, mp.Vector3(1, 0, 0), 1e-4, 2.5*fcen, 0.5*fcen, 3.48*fcen)
    
    geom_cr = [
        mp.Block(center=mp.Vector3(0, 0, -sz/4.0), size=mp.Vector3(mp.inf, mp.inf, sz/2.0), material=clad),
        mp.Block(center=mp.Vector3(0, 0, h/2.0), size=mp.Vector3(mp.inf, w, h), material=core),
        mp.Block(center=mp.Vector3(0, 0, h + t_pcm/2.0), size=mp.Vector3(mp.inf, w, t_pcm), material=pcm_cr),
    ]
    ms.geometry = geom_cr
    k_res_cr = ms.find_k(mp.NO_PARITY, fcen, 1, 1, mp.Vector3(1, 0, 0), 1e-4, 2.5*fcen, 0.5*fcen, 3.48*fcen)

    # 2. Coupled waveguides
    y_sep = w + gap
    y1, y2 = y_sep / 2.0, -y_sep / 2.0
    geom_dc = [
        mp.Block(center=mp.Vector3(0, 0, -sz/4.0), size=mp.Vector3(mp.inf, mp.inf, sz/2.0), material=clad),
        mp.Block(center=mp.Vector3(0, y1, h/2.0), size=mp.Vector3(mp.inf, w, h), material=core),
        mp.Block(center=mp.Vector3(0, y2, h/2.0), size=mp.Vector3(mp.inf, w, h), material=core),
        mp.Block(center=mp.Vector3(0, y1, h + t_pcm/2.0), size=mp.Vector3(mp.inf, w, t_pcm), material=pcm_am),
        mp.Block(center=mp.Vector3(0, y2, h + t_pcm/2.0), size=mp.Vector3(mp.inf, w, t_pcm), material=pcm_am),
    ]
    ms_dc = mpb.ModeSolver(geometry=geom_dc, geometry_lattice=mp.Lattice(size=cell),
                           k_points=[mp.Vector3(2.5*fcen, 0, 0)], resolution=40,
                           default_material=clad, num_bands=4)
    k_res_dc = ms_dc.find_k(mp.NO_PARITY, fcen, 1, 2, mp.Vector3(1, 0, 0), 1e-4, 2.5*fcen, 0.5*fcen, 3.48*fcen)
    sys.stdout.close()
    sys.stdout = old_stdout

    nam = k_res_am[0] / fcen
    ncr = k_res_cr[0] / fcen
    ne = max(k_res_dc[0], k_res_dc[1]) / fcen
    no = min(k_res_dc[0], k_res_dc[1]) / fcen
    split = ne - no
    lc = lambda_0 / (2.0 * split)
    print(f"Janus Baseline (w=280nm, h=220nm, gap=80nm, t_pcm=20nm):")
    print(f"  n_am = {nam:.5f}, n_cr = {ncr:.5f}, dn = {ncr - nam:.5f}")
    print(f"  n_even = {ne:.5f}, n_odd = {no:.5f}, dn_split = {split:.5f}")
    print(f"  Coupling Length L_c = {lc:.3f} um")

check_janus_baseline()
