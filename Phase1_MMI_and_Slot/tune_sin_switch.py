import sys, os, math
import meep as mp
from meep import mpb

def check_sin_width_scaling():
    lambda_0 = 1.064
    fcen = 1.0 / lambda_0
    clad = mp.Medium(index=1.449)
    core = mp.Medium(index=2.016)
    h = 0.300
    t_pcm = 0.080
    pcm_am = mp.Medium(index=2.70)
    pcm_cr = mp.Medium(index=3.30)

    print(f"{'W (nm)':<8} {'Gap (nm)':<10} {'dn_eff':<10} {'dn_split':<10} {'L_c (um)':<10} {'u':<8} {'Area (um2)':<12} {'XT(Cr) dB':<12}")
    print("-" * 85)

    for w in [0.350, 0.380, 0.400, 0.450]:
        # Single mode
        sy, sz = 3.5, 3.0
        cell = mp.Vector3(0, sy, sz)
        geom_am = [
            mp.Block(center=mp.Vector3(0, 0, -sz/4.0), size=mp.Vector3(mp.inf, mp.inf, sz/2.0), material=clad),
            mp.Block(center=mp.Vector3(0, 0, h/2.0), size=mp.Vector3(mp.inf, w, h), material=core),
            mp.Block(center=mp.Vector3(0, 0, h + t_pcm/2.0), size=mp.Vector3(mp.inf, w, t_pcm), material=pcm_am),
        ]
        ms = mpb.ModeSolver(geometry=geom_am, geometry_lattice=mp.Lattice(size=cell),
                            k_points=[mp.Vector3(1.75*fcen, 0, 0)], resolution=35,
                            default_material=clad, num_bands=2)
        old_stdout = sys.stdout
        sys.stdout = open(os.devnull, 'w')
        k_am = ms.find_k(mp.NO_PARITY, fcen, 1, 1, mp.Vector3(1, 0, 0), 1e-4, 1.75*fcen, 0.5*fcen, 2.5*fcen)
        geom_cr = [
            mp.Block(center=mp.Vector3(0, 0, -sz/4.0), size=mp.Vector3(mp.inf, mp.inf, sz/2.0), material=clad),
            mp.Block(center=mp.Vector3(0, 0, h/2.0), size=mp.Vector3(mp.inf, w, h), material=core),
            mp.Block(center=mp.Vector3(0, 0, h + t_pcm/2.0), size=mp.Vector3(mp.inf, w, t_pcm), material=pcm_cr),
        ]
        ms.geometry = geom_cr
        k_cr = ms.find_k(mp.NO_PARITY, fcen, 1, 1, mp.Vector3(1, 0, 0), 1e-4, 1.75*fcen, 0.5*fcen, 2.5*fcen)
        sys.stdout.close()
        sys.stdout = old_stdout

        nam = k_am[0] / fcen
        ncr = k_cr[0] / fcen
        dn_eff = ncr - nam

        for gap in [0.050, 0.080, 0.100]:
            y_sep = w + gap
            geom_dc = [
                mp.Block(center=mp.Vector3(0, 0, -sz/4.0), size=mp.Vector3(mp.inf, mp.inf, sz/2.0), material=clad),
                mp.Block(center=mp.Vector3(0, y_sep/2.0, h/2.0), size=mp.Vector3(mp.inf, w, h), material=core),
                mp.Block(center=mp.Vector3(0, -y_sep/2.0, h/2.0), size=mp.Vector3(mp.inf, w, h), material=core),
                mp.Block(center=mp.Vector3(0, y_sep/2.0, h + t_pcm/2.0), size=mp.Vector3(mp.inf, w, t_pcm), material=pcm_am),
                mp.Block(center=mp.Vector3(0, -y_sep/2.0, h + t_pcm/2.0), size=mp.Vector3(mp.inf, w, t_pcm), material=pcm_am),
            ]
            ms_dc = mpb.ModeSolver(geometry=geom_dc, geometry_lattice=mp.Lattice(size=cell),
                                   k_points=[mp.Vector3(1.75*fcen, 0, 0)], resolution=35,
                                   default_material=clad, num_bands=4)
            sys.stdout = open(os.devnull, 'w')
            k_dc = ms_dc.find_k(mp.NO_PARITY, fcen, 1, 2, mp.Vector3(1, 0, 0), 1e-4, 1.75*fcen, 0.5*fcen, 2.5*fcen)
            sys.stdout.close()
            sys.stdout = old_stdout
            ne = max(k_dc[0], k_dc[1]) / fcen
            no = min(k_dc[0], k_dc[1]) / fcen
            split = ne - no
            lc = lambda_0 / (2.0 * split)
            u = dn_eff / split
            arg = 0.5 * math.pi * math.sqrt(1.0 + u**2)
            p_cross_cr = (1.0 / (1.0 + u**2)) * (math.sin(arg)**2)
            p_bar_cr = 1.0 - p_cross_cr
            xt_cr = 10.0 * math.log10(max(p_cross_cr / max(p_bar_cr, 1e-9), 1e-9))
            
            # Total cell footprint:
            cell_l = 0.35 + lc + 0.85 + 0.35
            cell_w = 0.70 + 2 * 0.20
            area = cell_l * cell_w
            print(f"{w*1000:<8.0f} {gap*1000:<10.0f} {dn_eff:<10.4f} {split:<10.4f} {lc:<10.3f} {u:<8.3f} {area:<12.3f} {xt_cr:<12.2f}")

if __name__ == '__main__':
    check_sin_width_scaling()
