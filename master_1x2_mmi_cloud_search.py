"""
PROJECT JANUS: MASTER CLOUD SEARCH FOR ULTRA-COMPACT 1x2 MMI OPTICAL SWITCH
=============================================================================
Optimized for Multi-Core Cloud Compute (Azure 4-vCPU or HPC Clusters)

Architecture:
  - Device: Single-Cavity 1x2 Multimode Interference (MMI) Optical Switch
  - Input: 1 single optical port (TE0 fundamental mode at lambda_0 = 1064 nm)
  - Multimode cavity: Standard 220 nm SOI platform (Si core, SiO2 cladding)
  - Active Region: Asymmetric Sb2S3 phase-change thin film atop multimode cavity
  - Outputs: 2 isolated output ports (Port 1 = Bar, Port 2 = Cross)

Target Performance Constraints:
  1. Routing: Exactly 1 designated port destination per state (State 0 -> Port 2 Cross, State 1 -> Port 1 Bar)
  2. Crosstalk: <= -20.0 dB (Unselected port power leakage < 1.0%)
  3. Insertion Loss: <= 0.08 dB (Total Transmission >= 98.18%) for BOTH states
  4. Active Footprint: Aggressively minimized (< 5.0 um^2, prioritizing Sub-2.0 um^2)
  5. 100% Truth: Direct Poynting flux integration in MIT MEEP FDTD

Usage:
  python master_1x2_mmi_cloud_search.py --cpus 4 --mode both --output-dir .
"""

import sys
import os

# Optimize parallel CPU execution across multiple processes
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["NUMEXPR_NUM_THREADS"] = "1"

import time
import math
import argparse
import multiprocessing as mp_proc
import json
import csv
from typing import Dict, Any, List, Tuple
import numpy as np

try:
    import meep as mp
    HAS_MEEP = True
except ImportError:
    HAS_MEEP = False


class SingleCavity1x2MMI:
    """
    Parametric geometry and full 3D Vector Maxwell FDTD solver for a
    Single-Cavity 1x2 MMI Optical Switch with asymmetric Sb2S3 patch.
    """
    def __init__(
        self,
        w_mmi: float = 0.88,       # um (MMI cavity width, strictly TMI regime)
        l_mmi: float = 2.70,       # um (MMI cavity length)
        w_wg: float = 0.250,       # um (Access waveguide width)
        h_wg: float = 0.220,       # um (Silicon core height)
        w_tap: float = 0.380,      # um (Adiabatic taper width at MMI junction)
        l_tap: float = 0.300,      # um (Adiabatic taper length)
        y_port: float = 0.260,     # um (Output port center offset, +/- y_port)
        t_pcm: float = 0.150,      # um (Sb2S3 thickness, 100-180 nm)
        l_pcm: float = 2.50,       # um (Sb2S3 patch length)
        w_pcm: float = 0.440,      # um (Sb2S3 patch width)
        y_pcm_center: float = 0.220, # um (Sb2S3 patch lateral center, y > 0 for asymmetry)
        l_in: float = 0.350,       # um (Input access lead length)
        l_out: float = 0.350,      # um (Output access lead length)
        input_type: str = "offset" # "offset" (at +y_port) or "center" (at y = 0)
    ):
        self.lambda_0 = 1.064
        self.fcen = 1.0 / self.lambda_0
        self.df = 0.1 * self.fcen

        self.w_mmi = w_mmi
        self.l_mmi = l_mmi
        self.w_wg = w_wg
        self.h_wg = h_wg
        self.w_tap = w_tap
        self.l_tap = l_tap
        self.y_port = y_port
        self.t_pcm = t_pcm
        self.l_pcm = l_pcm
        self.w_pcm = w_pcm
        self.y_pcm_center = y_pcm_center
        self.l_in = l_in
        self.l_out = l_out
        self.input_type = input_type

        # Materials
        self.n_si = 3.480
        self.n_sio2 = 1.449
        self.n_pcm_am = 2.700
        self.k_pcm_am = 1.0e-5
        self.n_pcm_cr = 3.300
        self.k_pcm_cr = 8.0e-5

        # Footprint calculations
        self.total_length = self.l_in + 2.0 * self.l_tap + self.l_mmi + self.l_out
        self.total_width = max(self.w_mmi, 2.0 * self.y_port + self.w_wg) + 0.30
        self.footprint_area = self.total_length * self.total_width
        self.mmi_core_area = self.l_mmi * self.w_mmi

    def estimate_performance(self) -> Dict[str, Any]:
        """
        Fast analytical two-mode interference (TMI) estimation for initial screening.
        Calculates beat length L_pi, modal propagation constants, and phase accumulation.
        """
        n_eff_slab = 2.85
        sigma = 0  # TE polarization
        delta_w = (self.lambda_0 / np.pi) * (n_eff_slab**2 - self.n_sio2**2)**(-0.5)
        w_eff = self.w_mmi + delta_w

        # Fundamental beat length between mode 0 and mode 1
        l_pi = (4.0 * n_eff_slab * (w_eff**2)) / (3.0 * self.lambda_0)
        l_half_beat = l_pi / 2.0

        # Phase accumulation without PCM
        delta_beta_bare = math.pi / l_half_beat if l_half_beat > 0 else 0.0
        phase_bare = delta_beta_bare * self.l_mmi

        # Optical mode confinement factor in Sb2S3 thin film (up to 30% for 150 nm)
        gamma_pcm = min(0.35, 0.06 + 0.24 * (self.t_pcm / 0.150))
        delta_n_pcm = (self.n_pcm_cr - self.n_pcm_am) * gamma_pcm
        delta_phase_cr = (2.0 * math.pi / self.lambda_0) * delta_n_pcm * self.l_pcm

        # State 0: Amorphous
        phase_am = phase_bare
        p2_am_ideal = math.sin(phase_am / 2.0)**2
        p1_am_ideal = math.cos(phase_am / 2.0)**2

        # State 1: Crystalline
        phase_cr = phase_bare + delta_phase_cr
        p1_cr_ideal = math.cos(phase_cr / 2.0)**2
        p2_cr_ideal = math.sin(phase_cr / 2.0)**2

        # Estimated losses
        loss_scatter = 0.015  # dB modal transition
        loss_abs_am = (4 * math.pi * self.k_pcm_am / self.lambda_0) * self.l_pcm * 4.343 * gamma_pcm
        loss_abs_cr = (4 * math.pi * self.k_pcm_cr / self.lambda_0) * self.l_pcm * 4.343 * gamma_pcm

        il_am = -10.0 * math.log10(max(p2_am_ideal, 1e-4)) + loss_scatter + loss_abs_am
        il_cr = -10.0 * math.log10(max(p1_cr_ideal, 1e-4)) + loss_scatter + loss_abs_cr

        xt_am = 10.0 * math.log10(max(p1_am_ideal / max(p2_am_ideal, 1e-4), 1e-5))
        xt_cr = 10.0 * math.log10(max(p2_cr_ideal / max(p1_cr_ideal, 1e-4), 1e-5))

        promising = (self.footprint_area < 5.0 and il_am <= 0.35 and il_cr <= 0.35 and xt_am <= -10.0 and xt_cr <= -10.0)

        return {
            "l_pi": l_pi,
            "l_half_beat": l_half_beat,
            "footprint_area": self.footprint_area,
            "il_am_est": il_am,
            "il_cr_est": il_cr,
            "xt_am_est": xt_am,
            "xt_cr_est": xt_cr,
            "promising": promising
        }

    def run_meep_3d(self, state: str = "amorphous", resolution: int = 20) -> Dict[str, Any]:
        """
        Executes full-wave 3D Vector Maxwell MEEP FDTD simulation.
        Measures exact Poynting flux transmitted to Port 1 (Bar) and Port 2 (Cross).
        """
        if not HAS_MEEP:
            raise RuntimeError("MEEP is not installed in the active Python environment.")

        state_norm = state.strip().lower()
        is_crystalline = state_norm in ("crystalline", "cr", "1", "bar")

        clad_mat = mp.Medium(index=self.n_sio2)
        core_mat = mp.Medium(index=self.n_si)

        pcm_cond = 2 * math.pi * self.fcen * (self.k_pcm_cr if is_crystalline else self.k_pcm_am)
        pcm_idx = self.n_pcm_cr if is_crystalline else self.n_pcm_am
        pcm_mat = mp.Medium(index=pcm_idx, D_conductivity=pcm_cond)

        dpml = 0.40
        pad_x = 0.20
        pad_y = 0.20
        pad_z = 0.25
        sx = self.total_length + 2 * (dpml + pad_x)
        sy = self.total_width + 2 * (dpml + pad_y)
        sz = self.h_wg + self.t_pcm + 2 * (dpml + pad_z)

        cell = mp.Vector3(sx, sy, sz)
        pml_layers = [mp.PML(dpml)]

        x_left = -sx / 2.0
        x_right = sx / 2.0
        x_mmi_start = -self.l_mmi / 2.0
        x_mmi_end = self.l_mmi / 2.0
        x_tap_in = x_mmi_start - self.l_tap
        x_tap_out = x_mmi_end + self.l_tap

        y_in = self.y_port if self.input_type == "offset" else 0.0

        # Source placed in input access lead
        x_src = x_left + dpml + 0.08
        sources = [
            mp.EigenModeSource(
                src=mp.GaussianSource(self.fcen, fwidth=self.df),
                center=mp.Vector3(x_src, y_in, self.h_wg / 2.0),
                size=mp.Vector3(0, self.w_wg * 3.0, self.h_wg * 3.0),
                eig_band=1,
                direction=mp.X
            )
        ]

        # Flux monitor positions
        x_mon = x_tap_out + 0.15
        mon_w = min(0.260, 2.0 * self.y_port - 0.040)
        mon_h = 0.280

        # 1. Reference straight waveguide for normalization
        geom_ref = [
            mp.Block(center=mp.Vector3(0, 0, -sz/4.0), size=mp.Vector3(mp.inf, mp.inf, sz/2.0), material=clad_mat),
            mp.Block(center=mp.Vector3(0, y_in, self.h_wg/2.0), size=mp.Vector3(mp.inf, self.w_wg, self.h_wg), material=core_mat)
        ]

        sim_ref = mp.Simulation(
            cell_size=cell, boundary_layers=pml_layers, geometry=geom_ref,
            sources=sources, resolution=resolution, default_material=clad_mat
        )
        flux_ref = sim_ref.add_flux(self.fcen, 0, 1, mp.FluxRegion(center=mp.Vector3(x_mon, y_in, self.h_wg/2.0), size=mp.Vector3(0, mon_w, mon_h)))
        sim_ref.run(until_after_sources=mp.stop_when_fields_decayed(15, mp.Hz, mp.Vector3(x_mon, y_in, self.h_wg/2.0), 1e-3))
        p_ref = max(mp.get_fluxes(flux_ref)[0], 1e-12)

        # 2. Single-Cavity 1x2 TMI Switch Structure with Adiabatic Tapers
        geom = [
            mp.Block(center=mp.Vector3(0, 0, -sz/4.0), size=mp.Vector3(mp.inf, mp.inf, sz/2.0), material=clad_mat),
            # Input access lead (w_wg)
            mp.Block(
                center=mp.Vector3((x_left + x_tap_in)/2.0, y_in, self.h_wg/2.0),
                size=mp.Vector3(x_tap_in - x_left, self.w_wg, self.h_wg),
                material=core_mat
            ),
            # Input adiabatic taper (w_wg -> w_tap)
            mp.Prism(
                vertices=[
                    mp.Vector3(x_tap_in, y_in - self.w_wg / 2.0),
                    mp.Vector3(x_mmi_start, y_in - self.w_tap / 2.0),
                    mp.Vector3(x_mmi_start, y_in + self.w_tap / 2.0),
                    mp.Vector3(x_tap_in, y_in + self.w_wg / 2.0)
                ],
                height=self.h_wg,
                axis=mp.Vector3(0, 0, 1),
                material=core_mat
            ),
            # Multimode Cavity Body (Two-Mode Interference)
            mp.Block(
                center=mp.Vector3(0, 0, self.h_wg/2.0),
                size=mp.Vector3(self.l_mmi, self.w_mmi, self.h_wg),
                material=core_mat
            ),
            # Asymmetric Sb2S3 Patch atop MMI upper half
            mp.Block(
                center=mp.Vector3(0, self.y_pcm_center, self.h_wg + self.t_pcm / 2.0),
                size=mp.Vector3(self.l_pcm, self.w_pcm, self.t_pcm),
                material=pcm_mat
            ),
            # Output Port 1 Taper (w_tap -> w_wg at +y_port)
            mp.Prism(
                vertices=[
                    mp.Vector3(x_mmi_end, self.y_port - self.w_tap / 2.0),
                    mp.Vector3(x_tap_out, self.y_port - self.w_wg / 2.0),
                    mp.Vector3(x_tap_out, self.y_port + self.w_wg / 2.0),
                    mp.Vector3(x_mmi_end, self.y_port + self.w_tap / 2.0)
                ],
                height=self.h_wg,
                axis=mp.Vector3(0, 0, 1),
                material=core_mat
            ),
            # Output Port 2 Taper (w_tap -> w_wg at -y_port)
            mp.Prism(
                vertices=[
                    mp.Vector3(x_mmi_end, -self.y_port - self.w_tap / 2.0),
                    mp.Vector3(x_tap_out, -self.y_port - self.w_wg / 2.0),
                    mp.Vector3(x_tap_out, -self.y_port + self.w_wg / 2.0),
                    mp.Vector3(x_mmi_end, -self.y_port + self.w_tap / 2.0)
                ],
                height=self.h_wg,
                axis=mp.Vector3(0, 0, 1),
                material=core_mat
            ),
            # Output Port 1 Lead (+y_port, Bar)
            mp.Block(
                center=mp.Vector3((x_tap_out + x_right)/2.0, self.y_port, self.h_wg/2.0),
                size=mp.Vector3(x_right - x_tap_out, self.w_wg, self.h_wg),
                material=core_mat
            ),
            # Output Port 2 Lead (-y_port, Cross)
            mp.Block(
                center=mp.Vector3((x_tap_out + x_right)/2.0, -self.y_port, self.h_wg/2.0),
                size=mp.Vector3(x_right - x_tap_out, self.w_wg, self.h_wg),
                material=core_mat
            )
        ]

        sim = mp.Simulation(
            cell_size=cell, boundary_layers=pml_layers, geometry=geom,
            sources=sources, resolution=resolution, default_material=clad_mat
        )
        flux_p1 = sim.add_flux(self.fcen, 0, 1, mp.FluxRegion(center=mp.Vector3(x_mon, self.y_port, self.h_wg/2.0), size=mp.Vector3(0, mon_w, mon_h)))
        flux_p2 = sim.add_flux(self.fcen, 0, 1, mp.FluxRegion(center=mp.Vector3(x_mon, -self.y_port, self.h_wg/2.0), size=mp.Vector3(0, mon_w, mon_h)))

        sim.run(until_after_sources=mp.stop_when_fields_decayed(15, mp.Hz, mp.Vector3(x_mon, self.y_port, self.h_wg/2.0), 1e-3))

        p1 = max(mp.get_fluxes(flux_p1)[0] / p_ref, 1e-12)
        p2 = max(mp.get_fluxes(flux_p2)[0] / p_ref, 1e-12)
        tot = p1 + p2

        # Target port definition
        # State 0 (Amorphous) -> Port 2 (Cross)
        # State 1 (Crystalline) -> Port 1 (Bar)
        if not is_crystalline:
            p_target = p2
            p_unsel = p1
            dest_port = "Port 2 (Cross)"
        else:
            p_target = p1
            p_unsel = p2
            dest_port = "Port 1 (Bar)"

        il_db = -10.0 * math.log10(max(p_target, 1e-6))
        xt_db = 10.0 * math.log10(max(p_unsel / max(p_target, 1e-6), 1e-6))
        er_db = -xt_db

        return {
            "state": "crystalline" if is_crystalline else "amorphous",
            "destination_port": dest_port,
            "p1": p1,
            "p2": p2,
            "p_target": p_target,
            "p_unsel": p_unsel,
            "p1_pct": p1 * 100.0,
            "p2_pct": p2 * 100.0,
            "target_trans_pct": p_target * 100.0,
            "unsel_leak_pct": p_unsel * 100.0,
            "total_trans_pct": tot * 100.0,
            "insertion_loss_db": il_db,
            "crosstalk_db": xt_db,
            "extinction_ratio_db": er_db,
            "meets_il": il_db <= 0.08,
            "meets_xt": xt_db <= -20.0
        }


def evaluate_single_config_worker(args_tuple) -> Dict[str, Any]:
    """
    Multiprocessing worker function. Evaluates a single parameter combination.
    """
    params, resolution, run_full3d = args_tuple
    switch = SingleCavity1x2MMI(
        w_mmi=params["w_mmi"],
        l_mmi=params["l_mmi"],
        w_wg=params.get("w_wg", 0.250),
        y_port=params["y_port"],
        t_pcm=params["t_pcm"],
        l_pcm=params["l_pcm"],
        w_pcm=params.get("w_pcm", params["w_mmi"] / 2.0),
        y_pcm_center=params.get("y_pcm_center", params["w_mmi"] / 4.0),
        input_type=params.get("input_type", "offset")
    )

    screening = switch.estimate_performance()

    result = {
        "params": params,
        "footprint_area": switch.footprint_area,
        "mmi_core_area": switch.mmi_core_area,
        "screening": screening
    }

    if not run_full3d:
        result["sim_evaluated"] = False
        result["all_passed"] = screening["promising"]
        return result

    # Run 3D MEEP simulations for both states
    try:
        res_am = switch.run_meep_3d(state="amorphous", resolution=resolution)
        res_cr = switch.run_meep_3d(state="crystalline", resolution=resolution)

        all_passed = (
            res_am["meets_il"] and res_cr["meets_il"] and
            res_am["meets_xt"] and res_cr["meets_xt"] and
            switch.footprint_area < 5.0
        )

        result.update({
            "sim_evaluated": True,
            "all_passed": all_passed,
            "res_amorphous": res_am,
            "res_crystalline": res_cr,
            "il_am_db": res_am["insertion_loss_db"],
            "il_cr_db": res_cr["insertion_loss_db"],
            "xt_am_db": res_am["crosstalk_db"],
            "xt_cr_db": res_cr["crosstalk_db"],
            "leak_am_pct": res_am["unsel_leak_pct"],
            "leak_cr_pct": res_cr["unsel_leak_pct"],
            "trans_am_pct": res_am["target_trans_pct"],
            "trans_cr_pct": res_cr["target_trans_pct"],
            "fabric_area_mm2": (3_932_160 * switch.footprint_area) / 1e6
        })
    except Exception as e:
        result.update({
            "sim_evaluated": False,
            "error": str(e),
            "all_passed": False
        })

    return result


def generate_search_space(quick: bool = False) -> List[Dict[str, Any]]:
    """
    Generates the discrete grid of candidate parameters for the Single-Cavity 1x2 TMI Switch.
    """
    if quick:
        w_mmi_list = [0.85, 0.88]
        l_mmi_list = [2.60, 2.80]
        y_port_list = [0.260]
        t_pcm_list = [0.150]
        l_pcm_list = [2.40, 2.60]
    else:
        # Calibrated TMI exploratory space (Strictly 2-mode cavity)
        w_mmi_list = [0.82, 0.85, 0.88, 0.92, 0.96]
        l_mmi_list = [2.40, 2.60, 2.80, 3.00, 3.20]
        y_port_list = [0.240, 0.260, 0.280]
        t_pcm_list = [0.120, 0.150, 0.180]
        l_pcm_list = [2.20, 2.40, 2.60, 2.80]

    grid = []
    for w in w_mmi_list:
        for l in l_mmi_list:
            for y_p in y_port_list:
                if y_p >= w / 2.0 - 0.05:
                    continue  # Port must be safely inside the MMI boundary with cladding margin
                for t in t_pcm_list:
                    for l_p in l_pcm_list:
                        if l_p >= l:
                            continue  # Patch must fit on the MMI cavity
                        grid.append({
                            "w_mmi": w,
                            "l_mmi": l,
                            "y_port": y_p,
                            "t_pcm": t,
                            "l_pcm": l_p,
                            "w_pcm": w / 2.0,
                            "y_pcm_center": w / 4.0,
                            "w_wg": 0.250,
                            "w_tap": 0.380,
                            "l_tap": 0.300,
                            "input_type": "offset"
                        })
    return grid


def main():
    parser = argparse.ArgumentParser(description="Master Cloud Search for 1x2 MMI Optical Switch")
    parser.add_argument("--cpus", type=int, default=max(4, mp_proc.cpu_count()), help="Number of parallel CPU workers")
    parser.add_argument("--resolution", type=int, default=20, help="MEEP spatial grid resolution")
    parser.add_argument("--output-dir", type=str, default=".", help="Directory to save output files")
    parser.add_argument("--quick", action="store_true", help="Quick mode with reduced grid size")
    parser.add_argument("--mode", type=str, choices=["screen", "full3d", "both"], default="both",
                        help="Execution mode: 'screen' (fast eigenmode), 'full3d' (direct 3D FDTD), 'both' (screen then 3D)")
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)
    csv_path = os.path.join(args.output_dir, "single_cavity_1x2_mmi_results.csv")
    json_path = os.path.join(args.output_dir, "top_1x2_mmi_candidates.json")
    report_path = os.path.join(args.output_dir, "1X2_MMI_OPTIMIZATION_REPORT.md")

    print("=" * 80)
    print("  PROJECT JANUS: MASTER CLOUD SEARCH FOR SINGLE-CAVITY 1x2 MMI OPTICAL SWITCH")
    print(f"  Target: Footprint < 5.0 um^2 (Sub-2 um^2 prioritized) | IL <= 0.08 dB | XT <= -20 dB")
    print(f"  CPU Workers: {args.cpus} cores | Resolution: {args.resolution} | Mode: {args.mode.upper()}")
    print("=" * 80)

    grid = generate_search_space(quick=args.quick)
    print(f"[INFO] Generated search space: {len(grid)} candidate configurations.")

    t_start = time.time()

    # Phase 1: Fast Screening
    print("\n--- PHASE 1: HIGH-THROUGHPUT ANALYTICAL SCREENING ---")
    screen_tasks = [(p, args.resolution, False) for p in grid]
    with mp_proc.Pool(processes=args.cpus) as pool:
        screen_results = pool.map(evaluate_single_config_worker, screen_tasks)

    # Sort all screened results by best estimated balance and footprint
    screen_results.sort(key=lambda x: (
        not x["screening"].get("promising", False),
        x["footprint_area"],
        x["screening"]["il_am_est"] + x["screening"]["il_cr_est"]
    ))

    promising = [r for r in screen_results if r["screening"].get("promising", False)]
    if not promising:
        promising = screen_results[:20]  # Fallback to top 20 screened

    print(f"[INFO] Screening complete in {time.time() - t_start:.2f} s. {len(promising)} top configurations selected.")

    if args.mode == "screen":
        print(f"[INFO] Screening mode completed. Top candidate area: {promising[0]['footprint_area']:.3f} um^2")
        return

    # Phase 2: Full 3D Vector Maxwell MEEP FDTD Verification
    candidates_to_simulate = promising[:20] if args.mode == "both" else [r for r in screen_results]
    print(f"\n--- PHASE 2: 3D VECTOR MAXWELL MEEP FDTD VERIFICATION ({len(candidates_to_simulate)} candidates) ---")

    sim_tasks = [(c["params"], args.resolution, True) for c in candidates_to_simulate]
    sim_results = []

    # Run in parallel across available CPUs
    with mp_proc.Pool(processes=args.cpus) as pool:
        for idx, res in enumerate(pool.imap_unordered(evaluate_single_config_worker, sim_tasks), 1):
            sim_results.append(res)
            p = res["params"]
            if res.get("sim_evaluated"):
                passed_str = "PASS" if res.get("all_passed") else "FAIL"
                print(f"[{idx}/{len(sim_tasks)}] Area={res['footprint_area']:.3f} um^2 | "
                      f"Am: IL={res.get('il_am_db', 0):.3f}dB (XT={res.get('xt_am_db', 0):.1f}dB) | "
                      f"Cr: IL={res.get('il_cr_db', 0):.3f}dB (XT={res.get('xt_cr_db', 0):.1f}dB) -> [{passed_str}]")

    # Sort results
    sim_results.sort(key=lambda x: (
        not x.get("all_passed", False),
        x.get("footprint_area", 999),
        x.get("il_am_db", 99) + x.get("il_cr_db", 99)
    ))

    # Save to CSV
    with open(csv_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow([
            "w_mmi_um", "l_mmi_um", "y_port_um", "t_pcm_nm", "l_pcm_um",
            "footprint_area_um2", "fabric_area_mm2", "all_passed",
            "il_am_db", "il_cr_db", "xt_am_db", "xt_cr_db",
            "trans_am_pct", "trans_cr_pct", "leak_am_pct", "leak_cr_pct"
        ])
        for r in sim_results:
            p = r["params"]
            writer.writerow([
                p["w_mmi"], p["l_mmi"], p["y_port"], int(p["t_pcm"] * 1000), p["l_pcm"],
                round(r["footprint_area"], 4), round(r.get("fabric_area_mm2", 0), 4),
                r.get("all_passed", False),
                round(r.get("il_am_db", 0), 4), round(r.get("il_cr_db", 0), 4),
                round(r.get("xt_am_db", 0), 2), round(r.get("xt_cr_db", 0), 2),
                round(r.get("trans_am_pct", 0), 2), round(r.get("trans_cr_pct", 0), 2),
                round(r.get("leak_am_pct", 0), 2), round(r.get("leak_cr_pct", 0), 2)
            ])

    # Save top candidates to JSON
    top_candidates = [r for r in sim_results if r.get("all_passed", False)]
    with open(json_path, "w") as f:
        json.dump(top_candidates[:10], f, indent=2)

    total_time = time.time() - t_start
    print("\n" + "=" * 80)
    print("  SEARCH COMPLETE")
    print(f"  Total Duration: {total_time:.2f} s ({total_time / 60.0:.2f} min)")
    print(f"  Candidates Verified in 3D: {len(sim_results)}")
    print(f"  Winning Candidates Passing All Targets: {len(top_candidates)}")
    if top_candidates:
        win = top_candidates[0]
        print(f"\n  [TOP WINNER]")
        print(f"    Footprint Area   : {win['footprint_area']:.3f} um^2 (Core: {win['mmi_core_area']:.3f} um^2)")
        print(f"    Dimensions       : W_mmi = {win['params']['w_mmi']:.2f} um, L_mmi = {win['params']['l_mmi']:.2f} um")
        print(f"    Amorphous State  : IL = {win['il_am_db']:.4f} dB, Trans = {win['trans_am_pct']:.2f}%, XT = {win['xt_am_db']:.2f} dB (Leak = {win['leak_am_pct']:.2f}%)")
        print(f"    Crystalline State: IL = {win['il_cr_db']:.4f} dB, Trans = {win['trans_cr_pct']:.2f}%, XT = {win['xt_cr_db']:.2f} dB (Leak = {win['leak_cr_pct']:.2f}%)")
        print(f"    Fabric Area      : {win['fabric_area_mm2']:.2f} mm^2 (Baseline: 47.34 mm^2)")
    print("=" * 80)


if __name__ == "__main__":
    main()
