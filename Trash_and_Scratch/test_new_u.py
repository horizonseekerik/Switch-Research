import math
import meep as mp
from meep import mpb

lambda_0 = 1.064
fcen = 1.0 / lambda_0
clad = mp.Medium(index=1.449)
core = mp.Medium(index=3.480)
w_wg = 0.240
h_wg = 0.220
cell = mp.Vector3(0, 4.0, 3.0)

for t_pcm in [0.180, 0.200]:
    # 1. Delta_beta
    geom_am = [
        mp.Block(center=mp.Vector3(0, 0, -0.75), size=mp.Vector3(mp.inf, mp.inf, 1.5), material=clad),
        mp.Block(center=mp.Vector3(0, 0, h_wg/2.0), size=mp.Vector3(mp.inf, w_wg, h_wg), material=core),
        mp.Block(center=mp.Vector3(0, 0, h_wg + t_pcm/2.0), size=mp.Vector3(mp.inf, w_wg, t_pcm), material=mp.Medium(index=2.70)),
    ]
    ms = mpb.ModeSolver(geometry=geom_am, geometry_lattice=mp.Lattice(size=cell),
                        k_points=[mp.Vector3(2.5*fcen, 0, 0)], resolution=30,
                        default_material=clad, num_bands=1)
    k_am = ms.find_k(mp.NO_PARITY, fcen, 1, 1, mp.Vector3(1, 0, 0), 1e-4, 2.5*fcen, 0.5*fcen, 3.48*fcen)[0]
    
    geom_cr = [
        mp.Block(center=mp.Vector3(0, 0, -0.75), size=mp.Vector3(mp.inf, mp.inf, 1.5), material=clad),
        mp.Block(center=mp.Vector3(0, 0, h_wg/2.0), size=mp.Vector3(mp.inf, w_wg, h_wg), material=core),
        mp.Block(center=mp.Vector3(0, 0, h_wg + t_pcm/2.0), size=mp.Vector3(mp.inf, w_wg, t_pcm), material=mp.Medium(index=3.30)),
    ]
    ms.geometry = geom_cr
    k_cr = ms.find_k(mp.NO_PARITY, fcen, 1, 1, mp.Vector3(1, 0, 0), 1e-4, 2.5*fcen, 0.5*fcen, 3.48*fcen)[0]
    
    nam = k_am / fcen
    ncr = k_cr / fcen
    dbeta = 2 * math.pi * fcen * (ncr - nam)
    print(f"\nt_pcm = {t_pcm*1000:.0f} nm: nam={nam:.4f}, ncr={ncr:.4f}, Delta_beta={dbeta:.4f} rad/um")

    # 2. Kappa for gaps
    for gap_nm in [98, 105, 110]:
        gap = gap_nm / 1000.0
        y_sep = w_wg + gap
        y1, y2 = y_sep / 2.0, -y_sep / 2.0
        geom_coupler = [
            mp.Block(center=mp.Vector3(0, 0, -0.75), size=mp.Vector3(mp.inf, mp.inf, 1.5), material=clad),
            mp.Block(center=mp.Vector3(0, y1, h_wg/2.0), size=mp.Vector3(mp.inf, w_wg, h_wg), material=core),
            mp.Block(center=mp.Vector3(0, y2, h_wg/2.0), size=mp.Vector3(mp.inf, w_wg, h_wg), material=core),
            mp.Block(center=mp.Vector3(0, y1, h_wg + t_pcm/2.0), size=mp.Vector3(mp.inf, w_wg, t_pcm), material=mp.Medium(index=2.70)),
            mp.Block(center=mp.Vector3(0, y2, h_wg + t_pcm/2.0), size=mp.Vector3(mp.inf, w_wg, t_pcm), material=mp.Medium(index=2.70)),
        ]
        ms.geometry = geom_coupler
        ms.num_bands = 4
        k_sym = ms.find_k(mp.NO_PARITY, fcen, 1, 2, mp.Vector3(1, 0, 0), 1e-4, 2.5*fcen, 0.5*fcen, 3.48*fcen)
        ne = max(k_sym[0], k_sym[1]) / fcen
        no = min(k_sym[0], k_sym[1]) / fcen
        kappa = (math.pi / lambda_0) * (ne - no)
        u = dbeta / (2 * kappa)
        lc_am = (math.pi / 2) / kappa
        print(f"  gap={gap_nm} nm: kappa={kappa:.4f} rad/um, u={u:.3f}, L_pi/2={lc_am:.3f} um")
