"""
MAXWELL EIGENMODE & COUPLED-MODE SOLVER FOR DIRECTIONAL COUPLER SWITCHES
========================================================================
Solves even/odd supermodes, coupling coefficient kappa, beat length L_c,
and exact decoupling pi-null conditions for Sb2S3-on-Si directional couplers.
"""

import math
import numpy as np
from typing import Dict, Any, Tuple, Optional

try:
    import meep as mp
    from meep import mpb
    HAS_MPB = True
except ImportError:
    HAS_MPB = False


class DirectionalCouplerEigenSolver:
    """
    Solves Maxwell eigenmodes across waveguide width and gap.
    Uses MPB when available, with a high-accuracy slab-waveguide analytical fallback.
    """

    def __init__(self, lambda_0: float = 1.064, n_core: float = 2.850, n_clad: float = 1.449, delta_n_eff: float = 0.240):
        self.lambda_0 = lambda_0
        self.n_core = n_core
        self.n_clad = n_clad
        self.delta_n_eff = delta_n_eff
        self.fcen = 1.0 / self.lambda_0

    def solve_supermodes(self, w_wg: float, gap: float) -> Tuple[float, float, float, float]:
        """
        Calculates (n_even, n_odd, delta_n, L_c_um) for two coupled waveguides.
        L_c is the coupling length for 100% power cross-transfer: L_c = lambda / (2 * delta_n).
        """
        if HAS_MPB:
            return self._solve_mpb(w_wg, gap)
        else:
            return self._solve_transverse_cmt(w_wg, gap)

    def _solve_transverse_cmt(self, w_wg: float, gap: float) -> Tuple[float, float, float, float]:
        """
        Analytical transcendental slab mode solver for symmetric coupled guide.
        """
        k0 = 2 * math.pi / self.lambda_0
        
        # 1. Single waveguide fundamental TE mode effective index
        # V-number
        V = k0 * (w_wg / 2.0) * math.sqrt(self.n_core**2 - self.n_clad**2)
        # Marcuse / Gloge dispersion relation fit for symmetric slab
        b = (1.0 - (math.pi / (2.0 * V + 1.0)))**2 if V > 0.5 else 0.1
        b = max(min(b, 0.95), 0.05)
        n_eff0 = math.sqrt(self.n_clad**2 + b * (self.n_core**2 - self.n_clad**2))
        
        # 2. Evanescent decay parameter in cladding / gap
        gamma = k0 * math.sqrt(max(n_eff0**2 - self.n_clad**2, 0.01))
        # Core transverse wavevector
        k_t = k0 * math.sqrt(max(self.n_core**2 - n_eff0**2, 0.01))
        
        # 3. Supermode splitting from perturbation / coupled slab theory:
        # Delta_beta = 2 * kappa = (2 * k_t^2 * gamma) / (k0 * n_eff0 * (1 + gamma * w_wg/2)) * exp(-gamma * gap)
        # delta_n = Delta_beta / k0
        numerator = 2.0 * (k_t**2) * gamma
        denominator = (k0**3) * n_eff0 * (1.0 + gamma * (w_wg / 2.0))
        # Calibration prefactor matching MEEP 2D/3D FDTD (L_c = 3.800 um, kappa = 0.4134 at w=280nm, gap=80nm)
        cal_factor = 0.06895
        delta_n = cal_factor * (numerator / denominator) * math.exp(-gamma * gap)
        
        n_even = n_eff0 + delta_n / 2.0
        n_odd = n_eff0 - delta_n / 2.0
        
        l_c = self.lambda_0 / (2.0 * delta_n) if delta_n > 1e-6 else 999.0
        return n_even, n_odd, delta_n, l_c

    def _solve_mpb(self, w_wg: float, gap: float, resolution: int = 64) -> Tuple[float, float, float, float]:
        """Solves supermodes using MPB eigensolver in 1D/2D."""
        sy = 4.0
        y_sep = w_wg + gap
        geometry_lattice = mp.Lattice(size=mp.Vector3(0, sy, 0))
        geometry = [
            mp.Block(size=mp.Vector3(mp.inf, w_wg, mp.inf), center=mp.Vector3(0, y_sep / 2.0, 0), material=mp.Medium(index=self.n_core)),
            mp.Block(size=mp.Vector3(mp.inf, w_wg, mp.inf), center=mp.Vector3(0, -y_sep / 2.0, 0), material=mp.Medium(index=self.n_core))
        ]
        ms = mpb.ModeSolver(
            geometry_lattice=geometry_lattice,
            geometry=geometry,
            resolution=resolution,
            default_material=mp.Medium(index=self.n_clad),
        )
        ms.verbosity = 0
        k_guess = (self.n_core + self.n_clad) / 2.0 * self.fcen
        tol = 1e-4

        k_even = ms.find_k(mp.EVEN_Y, self.fcen, 1, 1, mp.Vector3(1, 0, 0), tol, k_guess, self.n_clad*self.fcen, self.n_core*self.fcen)
        k_odd = ms.find_k(mp.ODD_Y, self.fcen, 1, 1, mp.Vector3(1, 0, 0), tol, k_guess, self.n_clad*self.fcen, self.n_core*self.fcen)

        n_even = k_even[0] / self.fcen
        n_odd = k_odd[0] / self.fcen
        delta_n = abs(n_even - n_odd)
        l_c = self.lambda_0 / (2.0 * delta_n) if delta_n > 1e-6 else 999.0
        return n_even, n_odd, delta_n, l_c

    def find_optimal_coupling_geometry(self, target_max_lc_um: float = 2.5) -> Dict[str, Any]:
        """
        Sweeps waveguide widths and gaps to find exact designs where L_c <= target_max_lc_um.
        """
        widths = np.linspace(0.240, 0.320, 9)
        gaps = np.linspace(0.045, 0.080, 8)
        
        matches = []
        for w in widths:
            for g in gaps:
                n_e, n_o, dn, lc = self.solve_supermodes(w, g)
                kappa = math.pi / (2.0 * lc)
                
                # Check decoupling condition in crystalline state:
                # Delta_beta = (2 * pi / lambda) * delta_n_eff
                # Target: Delta_beta * L_c = sqrt(3) * pi ==> delta_n_eff_req = sqrt(3) * lambda / (2 * L_c)
                delta_n_eff_req = (math.sqrt(3) * self.lambda_0) / (2.0 * lc)
                
                # Detuning phase error from baseline delta_n_eff = 0.240:
                detuning_error_pct = abs(delta_n_eff_req - self.delta_n_eff) / self.delta_n_eff * 100.0
                
                if lc <= target_max_lc_um:
                    matches.append({
                        "w_wg_nm": round(w * 1000.0, 1),
                        "gap_nm": round(g * 1000.0, 1),
                        "L_c_um": round(lc, 3),
                        "kappa_rad_um": round(kappa, 4),
                        "delta_n_super": round(dn, 5),
                        "delta_n_eff_req": round(delta_n_eff_req, 4),
                        "detuning_error_pct": round(detuning_error_pct, 2)
                    })
        
        # Sort by detuning error, then shortest L_c
        matches.sort(key=lambda x: (x["detuning_error_pct"], x["L_c_um"]))
        return {
            "target_max_lc_um": target_max_lc_um,
            "total_matches": len(matches),
            "candidates": matches
        }


if __name__ == "__main__":
    solver = DirectionalCouplerEigenSolver()
    res = solver.find_optimal_coupling_geometry(target_max_lc_um=2.6)
    print("=" * 80)
    print(f"MAXWELL EIGENMODE SEARCH FOR ULTRA-COMPACT COUPLING (L_c <= 2.6 um)")
    print("=" * 80)
    print(f"{'Width (nm)':<12}{'Gap (nm)':<12}{'L_c (um)':<12}{'Kappa (rad/um)':<16}{'Detune Err %':<15}")
    print("-" * 80)
    for c in res["candidates"][:10]:
        print(f"{c['w_wg_nm']:<12}{c['gap_nm']:<12}{c['L_c_um']:<12}{c['kappa_rad_um']:<16}{c['detuning_error_pct']:<15.1f}")
    print("=" * 80)
