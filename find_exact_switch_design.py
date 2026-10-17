"""
PROJECT JANUS: EXACT SWITCH DESIGN SOLVER
Finds the exact (w1, w2, gap, L_c, t_pcm) for:
  1. Perfect Phase Matching in Amorphous State: n_eff1(am) = n_eff2(bare) -> 100% Cross
  2. Exact Delta_beta / 2*kappa = sqrt(3) -> 100% Bar in Crystalline State
"""
import math
import numpy as np
import meep as mp
from meep import mpb

def get_isolated_wg_neff(w, h=0.220, t_pcm=0.0, pcm_state="none", resolution=32):
    lambda_0 = 1.064
    fcen = 1.0 / lambda_0
    sy, sz = 2.0, 2.0
    cell = mp.Vector3(0, sy, sz)
    clad_mat = mp.Medium(index=1.449)
    core_mat = mp.Medium(index=3.480)

    geom = [
        mp.Block(center=mp.Vector3(0, 0, -sz/4.0), size=mp.Vector3(mp.inf, mp.inf, sz/2.0), material=clad_mat),
        mp.Block(center=mp.Vector3(0, 0, h/2.0), size=mp.Vector3(mp.inf, w, h), material=core_mat)
    ]
    if t_pcm > 0:
        n_pcm = 2.70 if pcm_state == "amorphous" else 3.30
        geom.append(mp.Block(
            center=mp.Vector3(0, 0, h + t_pcm/2.0),
            size=mp.Vector3(mp.inf, w, t_pcm),
            material=mp.Medium(index=n_pcm)
        ))

    k_guess = 2.45 * fcen
    ms = mpb.ModeSolver(
        geometry=geom,
        geometry_lattice=mp.Lattice(size=cell),
        k_points=[mp.Vector3(k_guess, 0, 0)],
        resolution=resolution,
        default_material=clad_mat,
        num_bands=2
    )
    k_res = ms.find_k(mp.NO_PARITY, fcen, 1, 1, mp.Vector3(1, 0, 0), 1e-4, k_guess, 0.5*fcen, 3.48*fcen)
    return k_res[0] / fcen

if __name__ == "__main__":
    import contextlib, os, sys
    @contextlib.contextmanager
    def nostdout():
        with open(os.devnull, 'w') as devnull:
            old = sys.stdout
            sys.stdout = devnull
            try: yield
            finally: sys.stdout = old

    print("=" * 70)
    print("MAPPING SINGLE-WAVEGUIDE EFFECTIVE INDICES vs WIDTH & PCM")
    print("=" * 70)

    # Let bare WG2 be w2 = 250 nm
    w2 = 0.250
    with nostdout():
        neff_w2 = get_isolated_wg_neff(w2, t_pcm=0.0)
    print(f"Bare WG2 (w2 = {w2*1000:.0f} nm): n_eff = {neff_w2:.5f}")

    # For different t_pcm (20nm, 40nm, 60nm, 80nm), find w1 that matches neff_w2 in amorphous state
    for t_pcm in [0.020, 0.040, 0.060, 0.080]:
        # Search w1 from 0.180 to 0.250
        print(f"\n--- Testing t_pcm = {t_pcm*1000:.0f} nm ---")
        best_w1 = None
        min_diff = 999.0
        best_neff_am = None
        best_neff_cr = None

        for w1 in np.linspace(0.180, 0.250, 15):
            with nostdout():
                neff_am = get_isolated_wg_neff(w1, t_pcm=t_pcm, pcm_state="amorphous")
            diff = abs(neff_am - neff_w2)
            if diff < min_diff:
                min_diff = diff
                best_w1 = w1
                best_neff_am = neff_am

        with nostdout():
            best_neff_cr = get_isolated_wg_neff(best_w1, t_pcm=t_pcm, pcm_state="crystalline")

        delta_n_cr = best_neff_cr - best_neff_am
        delta_beta = (2 * math.pi / 1.064) * delta_n_cr
        target_kappa = delta_beta / (2.0 * math.sqrt(3))
        ideal_lc = math.pi / (2.0 * target_kappa)

        print(f"  Phase-Matched w1 = {best_w1*1000:.1f} nm (n_eff_am = {best_neff_am:.5f}, diff = {min_diff:.5f})")
        print(f"  Crystalline n_eff= {best_neff_cr:.5f} -> Delta_n = {delta_n_cr:.5f}")
        print(f"  Delta_beta       = {delta_beta:.4f} um^-1")
        print(f"  Target kappa     = {target_kappa:.4f} um^-1")
        print(f"  Ideal L_c        = {ideal_lc:.2f} um")
        # Approximate footprint area
        cell_w = 0.600 + 2 * 0.150
        cell_l = ideal_lc + 0.35 + 0.80 + 0.35
        print(f"  Estimated Area   = {cell_l * cell_w:.2f} um^2 (L_tot = {cell_l:.2f} um)")
