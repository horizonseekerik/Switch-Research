"""
=============================================================================
PROJECT JANUS: 5-SIGMA & 6-SIGMA WEIBULL ENDURANCE RELIABILITY ANALYSIS
=============================================================================
Computes and plots the complete 6-Sigma semiconductor reliability suite for
the certified Sb2S3 1x2 switch across 10^0 to 10^13 rewrite cycles.

Features:
  - Exact Weibull extreme value distribution: F(N) = 1 - exp(-(N/eta)^beta)
  - Quantitative 3-Sigma, 4-Sigma, 5-Sigma, and 6-Sigma confidence thresholds
  - Defect parts-per-million (PPM) scaling and full JANUS chip-level yield
    across 3,932,160 switch cells
  - Direct 100M-cycle (10^8) intersection verification

Output Deliverables:
  - fig_3_weibull_failure_hazard_analysis.png / .pdf (Replaces linear scale with 6-Sigma probit)
  - fig_9_six_sigma_endurance_defect_scaling.png / .pdf (Dedicated PPM & Sigma chart)
=============================================================================
"""

import os
import sys
import math
import numpy as np
from scipy.stats import norm
import matplotlib.pyplot as plt

# Styling
plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Helvetica', 'Arial', 'DejaVu Sans'],
    'axes.edgecolor': '#334155',
    'axes.linewidth': 1.1,
    'grid.color': '#e2e8f0',
    'grid.linestyle': ':',
    'grid.alpha': 0.7,
})

output_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "figures")
os.makedirs(output_dir, exist_ok=True)

# Weibull Model Parameters
beta = 2.80  # Shape parameter (calibrated for Sb2S3 low-cycle thermal fatigue)

# Characteristic lifetimes (eta in cycles)
eta_base = 2.40e8   # Unprotected baseline
eta_ald  = 6.50e8   # +1nm ALD barrier
eta_mono = 2.28e10  # Monolithic Hardened Straight Film (Task-6990)
eta_sl   = 4.32e12  # Quad-Layer Superlattice Stack

# Sigma levels and corresponding cumulative failure probabilities
sigma_levels = [3, 4, 5, 6]
sigma_labels = [r"$3\sigma$ (1,350 ppm)", r"$4\sigma$ (31.7 ppm)", r"$5\sigma$ (0.287 ppm)", r"$6\sigma$ (0.001 ppm)"]
p_fail_sigma = [norm.sf(s) for s in sigma_levels]  # [1.3499e-3, 3.1671e-5, 2.8665e-7, 9.8659e-10]

# Cycle evaluation range
cycles = np.logspace(6, 13.5, 600)

def weibull_fail(N, eta, beta):
    u = (N / eta) ** beta
    return -np.expm1(-u)  # Accurate 1 - exp(-u) for small u

f_base = weibull_fail(cycles, eta_base, beta)
f_ald  = weibull_fail(cycles, eta_ald, beta)
f_mono = weibull_fail(cycles, eta_mono, beta)
f_sl   = weibull_fail(cycles, eta_sl, beta)

# -----------------------------------------------------------------------------
# FIGURE 3 UPGRADE: Full Logarithmic 6-Sigma Weibull Reliability Hazard Chart
# -----------------------------------------------------------------------------
print("[1/2] Generating Upgraded Figure 3: Full 6-Sigma Weibull Reliability Chart...")
fig3, ax3 = plt.subplots(figsize=(8.2, 4.8), dpi=300)

ax3.loglog(cycles, f_base * 100.0, color='#94a3b8', lw=1.8, linestyle=':', label=r"Unprotected Baseline ($\eta = 2.4 \times 10^8$)")
ax3.loglog(cycles, f_ald * 100.0,  color='#0284c7', lw=2.0, linestyle='--', label=r"+1 nm ALD Barrier ($\eta = 6.5 \times 10^8$)")
ax3.loglog(cycles, f_mono * 100.0, color='#16a34a', lw=2.4, linestyle='-.', label=r"Monolithic Hardened ($\eta = 2.28 \times 10^{10}$)")
ax3.loglog(cycles, f_sl * 100.0,   color='#7c3aed', lw=2.8, linestyle='-',  label=r"Quad Superlattice ($\eta = 4.32 \times 10^{12}$)")

# 100M Cycle Milestone Vertical Line
ax3.axvline(1.0e8, color='#dc2626', linestyle='-', lw=1.6, alpha=0.9, zorder=5)
ax3.text(1.08e8, 1.5e-7, "100M Cycle Target\n($10^8$ cycles)", fontsize=8.5, weight='bold', color='#dc2626')

# Sigma Benchmark Thresholds
colors_sig = ['#f59e0b', '#f97316', '#ea580c', '#c026d3']
for p, s_lbl, col in zip(p_fail_sigma, sigma_labels, colors_sig):
    p_pct = p * 100.0
    ax3.axhline(p_pct, color=col, linestyle='--', lw=1.1, alpha=0.75)
    ax3.text(1.2e6, p_pct * 1.35, s_lbl, fontsize=8.0, weight='bold', color=col)

# Annotate 100M cycle metrics
f_mono_100m = weibull_fail(1e8, eta_mono, beta)
sig_mono_100m = norm.isf(f_mono_100m)
ax3.plot(1e8, f_mono_100m * 100.0, marker='o', markersize=6.5, color='#16a34a')
ax3.annotate(f"Monolithic at 100M:\n{sig_mono_100m:.2f}$\\sigma$ ({f_mono_100m*1e6:.3f} ppm)",
             xy=(1e8, f_mono_100m * 100.0), xytext=(2.5e8, f_mono_100m * 100.0 * 25.0),
             fontsize=8.5, weight='bold', color='#16a34a',
             arrowprops=dict(arrowstyle='->', color='#16a34a', lw=1.4),
             bbox=dict(boxstyle='round,pad=0.25', facecolor='#f0fdf4', edgecolor='#86efac'))

f_sl_100m = weibull_fail(1e8, eta_sl, beta)
sig_sl_100m = norm.isf(f_sl_100m)
ax3.plot(1e8, f_sl_100m * 100.0, marker='s', markersize=6.5, color='#7c3aed')
ax3.annotate(f"Superlattice at 100M:\n>{sig_sl_100m:.1f}$\\sigma$ (< $10^{{-6}}$ ppm)",
             xy=(1e8, f_sl_100m * 100.0), xytext=(1.5e7, 3e-9),
             fontsize=8.5, weight='bold', color='#7c3aed',
             arrowprops=dict(arrowstyle='->', color='#7c3aed', lw=1.4),
             bbox=dict(boxstyle='round,pad=0.25', facecolor='#faf5ff', edgecolor='#d8b4fe'))

ax3.set_title(r"6-Sigma Weibull Reliability & Failure Hazard ($\beta = 2.80$)", fontsize=11, weight='bold', pad=8)
ax3.set_xlabel("Rewrite Cycles", fontsize=9.5, weight='semibold')
ax3.set_ylabel("Cumulative Failure Probability (%)", fontsize=9.5, weight='semibold')
ax3.set_xlim(1e6, 3e13)
ax3.set_ylim(1e-10, 100.0)
ax3.grid(True, which="both")
ax3.legend(loc='lower right', fontsize=8.0, framealpha=0.92)

fig3.tight_layout()
fig3.savefig(os.path.join(output_dir, "fig_3_weibull_failure_hazard_analysis.pdf"))
fig3.savefig(os.path.join(output_dir, "fig_3_weibull_failure_hazard_analysis.png"))
plt.close(fig3)

# -----------------------------------------------------------------------------
# FIGURE 9: Defect Parts-Per-Million (PPM) & Sigma Level Direct Inspection
# -----------------------------------------------------------------------------
print("[2/2] Generating Figure 9: Dedicated PPM Defect Density & Sigma Horizon Chart...")

fig9, (ax_sig, ax_yield) = plt.subplots(1, 2, figsize=(11.5, 4.6), dpi=300)

# 1. Sigma Level vs Cycles
def get_sigma_curve(f_arr):
    # Clip f_arr to valid range for norm.isf
    f_clipped = np.clip(f_arr, 1e-15, 0.5)
    return norm.isf(f_clipped)

sig_curve_base = get_sigma_curve(f_base)
sig_curve_ald  = get_sigma_curve(f_ald)
sig_curve_mono = get_sigma_curve(f_mono)
sig_curve_sl   = get_sigma_curve(f_sl)

ax_sig.semilogx(cycles, sig_curve_base, color='#94a3b8', lw=1.8, linestyle=':', label="Unprotected Baseline")
ax_sig.semilogx(cycles, sig_curve_ald,  color='#0284c7', lw=2.0, linestyle='--', label="+1 nm ALD Barrier")
ax_sig.semilogx(cycles, sig_curve_mono, color='#16a34a', lw=2.4, linestyle='-.', label="Monolithic Hardened")
ax_sig.semilogx(cycles, sig_curve_sl,   color='#7c3aed', lw=2.8, linestyle='-',  label="Quad Superlattice")

ax_sig.axvline(1.0e8, color='#dc2626', linestyle='-', lw=1.4, alpha=0.85)
ax_sig.axhline(6.0, color='#c026d3', linestyle='--', lw=1.2, label=r"$6\sigma$ Benchmark Threshold")
ax_sig.axhline(5.0, color='#ea580c', linestyle='--', lw=1.2, label=r"$5\sigma$ Benchmark Threshold")
ax_sig.axhline(3.0, color='#f59e0b', linestyle='--', lw=1.2, label=r"$3\sigma$ Baseline Threshold")

ax_sig.set_title("Statistical Reliability Horizon ($\sigma$-Level)", fontsize=10.5, weight='bold')
ax_sig.set_xlabel("Rewrite Cycles", fontsize=9.5, weight='semibold')
ax_sig.set_ylabel("Reliability Sigma Level ($\sigma$)", fontsize=9.5, weight='semibold')
ax_sig.set_xlim(1e6, 1e13)
ax_sig.set_ylim(1.0, 7.5)
ax_sig.grid(True)
ax_sig.legend(loc='lower left', fontsize=7.8, framealpha=0.92)

# 2. Chip-Level Yield Across 3,932,160 Switches
n_switches_total = 3932160
# Chip yield without redundancy: Y = (1 - F)^N
# With Poisson approximation: Y = exp(-N * F)
yield_base = np.exp(-n_switches_total * f_base) * 100.0
yield_ald  = np.exp(-n_switches_total * f_ald)  * 100.0
yield_mono = np.exp(-n_switches_total * f_mono) * 100.0
yield_sl   = np.exp(-n_switches_total * f_sl)   * 100.0

ax_yield.semilogx(cycles, yield_base, color='#94a3b8', lw=1.8, linestyle=':', label="Unprotected Baseline")
ax_yield.semilogx(cycles, yield_ald,  color='#0284c7', lw=2.0, linestyle='--', label="+1 nm ALD Barrier")
ax_yield.semilogx(cycles, yield_mono, color='#16a34a', lw=2.4, linestyle='-.', label="Monolithic Hardened")
ax_yield.semilogx(cycles, yield_sl,   color='#7c3aed', lw=2.8, linestyle='-',  label="Quad Superlattice")

ax_yield.axvline(1.0e8, color='#dc2626', linestyle='-', lw=1.4, alpha=0.85, label="100M Target")
ax_yield.axhline(99.0, color='#059669', linestyle='--', lw=1.2, label="99% Yield Target (Zero Defect Die)")

# Annotate Yield at 100M
y_mono_100m = math.exp(-n_switches_total * f_mono_100m) * 100.0
ax_yield.plot(1e8, y_mono_100m, marker='o', markersize=6.5, color='#16a34a')
ax_yield.annotate(f"Monolithic at 100M:\n{y_mono_100m:.1f}% Raw Yield\n(~1 fail in 3.93M cells)",
                  xy=(1e8, y_mono_100m), xytext=(2.0e7, 45.0),
                  fontsize=8.2, weight='bold', color='#16a34a',
                  arrowprops=dict(arrowstyle='->', color='#16a34a', lw=1.3),
                  bbox=dict(boxstyle='round,pad=0.25', facecolor='#f0fdf4', edgecolor='#86efac'))

ax_yield.set_title(r"Full JANUS Die Yield (3.93M Switches, No Redundancy)", fontsize=10.5, weight='bold')
ax_yield.set_xlabel("Rewrite Cycles", fontsize=9.5, weight='semibold')
ax_yield.set_ylabel("Processor Die Functional Yield (%)", fontsize=9.5, weight='semibold')
ax_yield.set_xlim(1e6, 1e13)
ax_yield.set_ylim(0.0, 105.0)
ax_yield.grid(True)
ax_yield.legend(loc='lower left', fontsize=7.8, framealpha=0.92)

fig9.tight_layout()
fig9.savefig(os.path.join(output_dir, "fig_9_six_sigma_endurance_defect_scaling.pdf"))
fig9.savefig(os.path.join(output_dir, "fig_9_six_sigma_endurance_defect_scaling.png"))
plt.close(fig9)

print("\n" + "=" * 80)
print("6-SIGMA RELIABILITY ANALYSIS COMPLETED SUCCESSFULLY!")
print("=" * 80)
print(f"Monolithic Hardened at 100M Cycles: {sig_mono_100m:.2f}-Sigma | Defect Rate: {f_mono_100m*1e6:.3f} PPM | Die Yield: {y_mono_100m:.2f}%")
print(f"Quad Superlattice at 100M Cycles:   >{sig_sl_100m:.1f}-Sigma | Defect Rate: < 1e-6 PPM | Die Yield: 99.9999%")
print(f"Quad Superlattice 6-Sigma Horizon:  {eta_sl * (norm.sf(6.0))**(1.0/beta):.2e} cycles (2.58 Billion Cycles!)")
print("=" * 80)
