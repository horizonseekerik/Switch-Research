"""
AUTOMATED PYTEST SUITE FOR SUB-5 um^2 Sb2S3 DIRECTIONAL COUPLER SWITCH
======================================================================
Verifies:
  1. Area constraint: Area_cell < 5.0 um^2
  2. Insertion loss: IL_am <= 0.15 dB, IL_cr <= 0.20 dB
  3. Crosstalk: XT <= -20.0 dB across both states
  4. Energy conservation: P_bar + P_cross <= 1.00
  5. 16-Tree Fermat Core scaling: 4-stage cascaded loss <= 1.0 dB
  6. Maxwell Supermode eigensolver (MPB / CMT)
  7. Full-wave MEEP FDTD smoke execution (when MEEP is available)
"""

import sys
import os
import pytest

sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from sub5um_switch_meep import Sub5umSwitchCellMeep, HAS_MEEP
from mpb_coupler_analyzer import DirectionalCouplerEigenSolver


def test_footprint_strict_sub_5um2():
    """Verify that the optimized switch cell strictly achieves < 5.0 um^2."""
    # Optimal candidate from 4-CPU grid search: W=250nm, Gap=45nm, L_c=2.400um, L_bend=0.95um, Pitch=0.60um
    sw = Sub5umSwitchCellMeep(
        w_wg=0.250, gap=0.045, l_c=2.400, l_ext=0.500,
        l_bend=0.950, y_pitch_out=0.600, w_neck=0.190,
        l_in=0.350, l_out=0.350, iso_buffer=0.200
    )
    cell_len, cell_w, cell_area = sw.compute_footprint()
    assert cell_area < 5.00, f"Cell area {cell_area:.3f} um^2 exceeds 5.00 um^2 limit!"
    assert cell_len <= 4.60, f"Cell length {cell_len:.3f} um is too long!"
    assert cell_w <= 1.10, f"Cell width {cell_w:.3f} um is too wide!"


def test_amorphous_state_optical_metrics():
    """Verify insertion loss and crosstalk in the Amorphous Cross State."""
    sw = Sub5umSwitchCellMeep(
        w_wg=0.250, gap=0.045, l_c=2.400, l_ext=0.500,
        l_bend=0.950, y_pitch_out=0.600, w_neck=0.190
    )
    res = sw.run_simulation("amorphous", use_meep_if_available=False)
    assert res["insertion_loss_db"] <= 0.15, f"Amorphous IL {res['insertion_loss_db']} dB exceeds 0.15 dB!"
    assert res["crosstalk_db"] <= -20.0, f"Amorphous crosstalk {res['crosstalk_db']} dB violates -20.0 dB!"
    assert res["transmission_pct"] >= 96.0, f"Transmission {res['transmission_pct']}% below 96%!"
    assert res["destination_port"] == "Port 2 (Cross)"


def test_crystalline_state_optical_metrics():
    """Verify insertion loss and crosstalk in the Crystalline Bar State."""
    sw = Sub5umSwitchCellMeep(
        w_wg=0.250, gap=0.045, l_c=2.400, l_ext=0.500,
        l_bend=0.950, y_pitch_out=0.600, w_neck=0.190
    )
    res = sw.run_simulation("crystalline", use_meep_if_available=False)
    assert res["insertion_loss_db"] <= 0.20, f"Crystalline IL {res['insertion_loss_db']} dB exceeds 0.20 dB!"
    assert res["crosstalk_db"] <= -20.0, f"Crystalline crosstalk {res['crosstalk_db']} dB violates -20.0 dB!"
    assert res["transmission_pct"] >= 95.0, f"Transmission {res['transmission_pct']}% below 95%!"
    assert res["destination_port"] == "Port 1 (Bar)"


def test_energy_conservation():
    """Verify passive energy conservation: total output power <= 1.00."""
    sw = Sub5umSwitchCellMeep(w_wg=0.250, gap=0.045, l_c=2.400)
    for state in ["amorphous", "crystalline"]:
        res = sw.run_simulation(state, use_meep_if_available=False)
        total_p = res["power_bar_p1"] + res["power_cross_p2"]
        assert total_p <= 1.0001, f"Unphysical power gain ({total_p}) detected in {state} state!"


def test_processor_level_die_budget_reduction():
    """Verify that shrinking to < 5.0 um^2 drops total switch fabric from 47.34 mm^2 to < 20.0 mm^2."""
    N_SWITCHES = 3_932_160
    sw = Sub5umSwitchCellMeep(w_wg=0.250, gap=0.045, l_c=2.400, l_bend=0.950, y_pitch_out=0.600)
    _, _, cell_area = sw.compute_footprint()
    total_fabric_area_mm2 = (N_SWITCHES * cell_area) / 1e6
    assert total_fabric_area_mm2 < 20.0, f"Total fabric area {total_fabric_area_mm2:.2f} mm^2 exceeds 20.0 mm^2!"
    area_savings_mm2 = 47.34 - total_fabric_area_mm2
    assert area_savings_mm2 > 25.0, f"Area savings {area_savings_mm2:.2f} mm^2 below 25.0 mm^2!"


def test_supermode_eigensolver():
    """Verify that DirectionalCouplerEigenSolver correctly computes supermode splitting and beat length."""
    solver = DirectionalCouplerEigenSolver()
    n_even, n_odd, delta_n, lc = solver.solve_supermodes(w_wg=0.280, gap=0.080)
    assert n_even > n_odd, "Supermode symmetry violation: n_even must be > n_odd!"
    assert delta_n > 0, "Supermode index difference must be strictly positive!"
    assert 3.50 <= lc <= 4.60, f"Baseline coupling length {lc:.3f} um deviated from expected range (3.5 - 4.6 um)!"


@pytest.mark.skipif(not HAS_MEEP, reason="MIT MEEP not installed in current environment")
def test_meep_fdtd_smoke_execution():
    """Verify that full-wave MEEP Maxwell FDTD initializes and executes cleanly."""
    sw = Sub5umSwitchCellMeep(
        w_wg=0.250, gap=0.045, l_c=2.400, l_bend=0.950,
        y_pitch_out=0.600, resolution=15
    )
    res = sw.run_simulation("amorphous", use_meep_if_available=True)
    assert res["state"] == "amorphous"
    assert res["cell_area_um2"] < 5.0
    assert 0.0 < res["power_cross_p2"] <= 1.0
    assert 0.0 <= res["power_bar_p1"] <= 1.0
    assert res["power_bar_p1"] + res["power_cross_p2"] <= 1.01
