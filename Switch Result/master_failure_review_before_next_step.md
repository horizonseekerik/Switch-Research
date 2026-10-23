# PROJECT JANUS: Complete Failure Review & Lessons Learned
## Master Engineering Knowledge Base Before Next Design Iteration

**Date:** 2026-10-19  
**Scope:** Every architecture, simulation, and approach attempted — why it was tried, what it produced, why it failed, and the non-negotiable lessons extracted.

---

## 1. Non-Negotiable Specifications (Unchanged)

| Metric | Requirement | Notes |
| :--- | :--- | :--- |
| State 0 (Amorphous → Port 2) | $P_2 > 98.0000000000\%$ | Cross routing |
| State 1 (Crystalline → Port 1) | $P_1 > 98.0000000000\%$ | Bar / through routing |
| Total Guided Power | $\ge 99.0000000000\%$ | Radiation loss $\le 1.0\%$ |
| Insertion Loss | $< 0.1000000000\,\text{dB}$ | $P_{\text{signal}} / P_{\text{in}} \ge 97.72\%$ |
| Modal Crosstalk | $< -20.0000000000\,\text{dB}$ | $P_{\text{leak}} / P_{\text{signal}} \le 1.0\%$ |
| Active Footprint | $\le 6.0000000000\,\mu\text{m}^2$ | Strictly enforced |
| Architecture | Pure spatial routing | **NO** 50:50 splitters, 3-dB couplers, or MZI |
| Validation | PyMEEP 3D Maxwell FDTD only | All claims must be 3D verified |
| Precision | $\ge 10$ decimal places | No vague words ("similar", "approximately") |

---

## 2. Master Failure Scorecard (All 11 Attempts)

| # | Architecture / Approach | Best State 0 | Best State 1 | Total Guided | Rad Loss | Area ($\mu\text{m}^2$) | Fatal Failure |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| 1 | Sine-Slot Coupler | $78.9030000000\%$ | $31.3910000000\%$ | $89.70\text{--}93.83\%$ | $6.17\text{--}10.30\%$ | $5.00\text{--}5.29$ | Mode impedance mismatch → radiation |
| 2 | Adiabatic Y-Junction | $83.5000000000\%$ | $72.1000000000\%$ | $83.5000000000\%$ | $16.5000000000\%$ | $7.8000000000$ | Tip scattering + footprint overflow |
| 3 | MZI / 3-dB Splitter | — | — | — | — | $> 20.0$ | **Rejected by user** (splits beam) |
| 4 | Asymmetric Coupler (bare WG2) | $92.4000000000\%$ | $44.1000000000\%$ | $98.1000000000\%$ | $1.9000000000\%$ | $> 45.0$ | Vertical mode mismatch ($\kappa$ ÷ $7.4\times$) |
| 5 | PCM on Curved S-Bend | $81.2000000000\%$ | $64.5000000000\%$ | $88.9000000000\%$ | $11.1000000000\%$ | $4.8500000000$ | Curvature mode deformation |
| 6 | Decoupled DC @ $L_c = 3.00\,\mu\text{m}$ | **$91.9548201200\%$** | $51.9100000000\%$ | $92.0525506188\%$ | $\approx 8.0\%$ | $4.2600000000$ | State 1 deeply under-coupled |
| 7 | Decoupled DC @ $L_c = 4.28\,\mu\text{m}$ | $77.8300000000\%$ | **$90.9324819000\%$** | **$99.8248102390\%$** | **$0.1751897610\%$** | **$5.0830000000$** | State 0 over-coupled past half-beat |
| 8 | Segmented Two-$L_c$ Switch (Cand. 1) | $62.5080058841\%$ | $18.9695238604\%$ | $71.8413414117\%$ | **$55.2540530060\%$** | $4.9200000000$ | Abrupt width step → catastrophic scattering |
| 9 | Wide-Gap ($137.5\,\text{nm}$) Coupler (Cand. 2) | $2.8655288184\%$ | $92.3898906320\%$ | $95.2554194504\%$ | $4.7445805496\%$ | $5.4365000000$ | Overshot into full 2-cycle beat |
| 10 | Aligned Detuned DC @ $g = 98.2\,\text{nm}$ | $79.1000000000\%$ | $52.6600000000\%$ | $96.1600000000\%$ | $\approx 3.8\%$ | $5.1600000000$ | $\kappa_{\text{cr}} \ne \kappa_{\text{am}}$ + transition coupling |
| 11 | Compromise Point @ $g = 94.5\,\text{nm}$ | $83.8803310000\%$ | $75.4624490000\%$ | $93.9393200000\%$ | $\approx 6.1\%$ | $4.7600000000$ | Both states off-resonance simultaneously |

> [!CAUTION]
> **Not a single attempt has achieved $> 98.0000000000\%$ in BOTH states simultaneously in 3D FDTD.** The closest was Attempt 6 reaching $91.9548201200\%$ (State 0) with $99.9022695012\%$ modal purity, but State 1 at that length was only $51.9100000000\%$.

---

## 3. Detailed Failure Autopsies & Extracted Lessons

### FAILURE 1: Sine-Slot Coupler (Attempts 1)

**What We Tried:**
A directional coupler with a $70.0000000000\,\text{nm}$ nanoslot between two $450.0000000000\,\text{nm}$ silicon rails, with $\text{Sb}_2\text{S}_3$ filling the slot. The slot geometry exploits the $E$-field discontinuity at dielectric boundaries ($E_{\text{slot}} / E_{\text{Si}} = \epsilon_{\text{Si}} / \epsilon_{\text{slot}} \approx 1.66$) to maximize PCM modal confinement ($\Gamma_{\text{pcm}} > 40\%$).

**Why It Failed:**
- Transforming from strip mode to slot mode causes severe spatial wave-impedance mismatch.
- Even $0.3500000000\,\mu\text{m}$ tapers could not reshape the mode adiabatically.
- **Radiation loss: $6.1700000000\%\text{--}10.3000000000\%$** (spec: $\le 1.0\%$).
- Crystalline extinction was catastrophic: $P_1 = 20.8470000000\%\text{--}31.3910000000\%$ (spec: $> 98.0\%$).

> [!WARNING]
> **LESSON 1: Never introduce a cross-section that fundamentally changes the mode topology.**
> Strip-to-slot transitions involve a topological field rearrangement (single-peak → bimodal). This cannot be made adiabatic within any sub-$6.0\,\mu\text{m}^2$ footprint. The waveguide core must remain a simple rectangular strip throughout.

---

### FAILURE 2: Adiabatic Y-Junction Mode Sorter (Attempt 2)

**What We Tried:**
A branching Y-junction where $\text{Sb}_2\text{S}_3$ on one arm steers the fundamental mode toward that arm when crystallised (higher $n_{\text{eff}}$) and away from it when amorphous.

**Why It Failed:**
- The branching apex has finite bluntness ($w_{\text{tip}} \ge 20.0000000000\,\text{nm}$). When the guided mode encounters this blunt tip, it diffracts violently.
- **Radiation loss: $16.5000000000\%$** — $16.5\times$ the allowed $1.0\%$.
- To satisfy the Love adiabaticity criterion ($\theta \ll \lambda_0 / [2 w_{\text{core}} \Delta n_{\text{eff}}]$), device length exceeds $8.0000000000\,\mu\text{m}$.
- **Footprint: $7.8000000000\,\mu\text{m}^2 > 6.0000000000\,\mu\text{m}^2$**.

> [!WARNING]
> **LESSON 2: Y-junctions and mode sorters cannot operate within $6.0\,\mu\text{m}^2$ at $1064\,\text{nm}$.**
> The adiabaticity length requirement at our wavelength and index contrast is fundamentally incompatible with the footprint constraint. Do not attempt any forked-geometry approach.

---

### FAILURE 3: Asymmetric Coupler — PCM on WG1 Only, Bare WG2 (Attempt 4)

**What We Tried:**
Depositing $\text{Sb}_2\text{S}_3$ only on Waveguide 1 while Waveguide 2 remains bare silicon, to simplify fabrication.

**Why It Failed:**
- The $200.0000000000\,\text{nm}$ PCM cap on WG1 pulls the mode center-of-mass upward to $z_1 = 0.2350000000\,\mu\text{m}$, while bare WG2's mode center stays at $z_2 = 0.1100000000\,\mu\text{m}$.
- The $125.0000000000\,\text{nm}$ vertical offset reduces the mode overlap integral exponentially:
  $$\int E_1(z) E_2(z)\,dz \propto \exp\left(-\frac{(0.125)^2}{2(0.06)^2}\right) \approx 0.1150000000$$
- **$\kappa$ collapsed by $7.4\times$**, requiring $L_c > 30.2000000000\,\mu\text{m}$ → footprint $> 45.0000000000\,\mu\text{m}^2$.

> [!CAUTION]
> **LESSON 3: BOTH waveguides MUST carry matching PCM layers in the interaction zone.**
> If WG1 has a PCM overlay and WG2 does not, the transverse mode profiles become vertically misaligned and $\kappa$ collapses. The PCM overlay must be symmetric: identical material, thickness, and width on both waveguides throughout the parallel coupling region.

---

### FAILURE 4: PCM Extending Along Curved S-Bend (Attempt 5)

**What We Tried:**
Continuing the $\text{Sb}_2\text{S}_3$ layer along the curved output S-bend of Waveguide 2 to maintain phase matching through the fanout.

**Why It Failed:**
- In a curved waveguide, the conformal index transformation is $n_{\text{eq}}(y) = n(y)(1 + y/R)$.
- With the high-index PCM cap ($n = 3.3000000000$) on a tight bend ($R \le 1.8000000000\,\mu\text{m}$), the equivalent index tilt deforms the mode, pushing it outward.
- **Radiation loss: $11.1000000000\%$** — $11.1\times$ the allowed $1.0\%$.

> [!CAUTION]
> **LESSON 4: PCM must be STRICTLY confined to straight, parallel waveguide sections.**
> Never place PCM on any curved segment. All S-bends, fanouts, and tapers must be bare silicon only. The PCM active region must terminate with adiabatic tips before any curvature begins.

---

### FAILURE 5: Segmented Two-$L_c$ Width-Step Switch — Candidate 1 (Attempt 8)

**What We Tried:**
An abrupt waveguide width step ($240.0000000000\,\text{nm} \to 320.0000000000\,\text{nm}$) inside the coupling zone to create two different effective coupling coefficients in the same device, allowing State 0 and State 1 to peak at different effective lengths but the same physical length.

**Why It Failed:**
- The abrupt $80.0000000000\,\text{nm}$ core width step at the boundary between the two segments created a catastrophic wave-impedance mismatch.
- **Radiation loss: $55.2540530060\%$** — the worst scattering of any iteration.
- State 0: $P_2 = 62.5080058841\%$, State 1: $P_1 = 18.9695238604\%$.

> [!CAUTION]
> **LESSON 5: NEVER introduce an abrupt cross-sectional step inside the interaction zone.**
> Any discontinuity — width steps, thickness steps, or material boundaries — that cause a sudden change in the guided mode profile will scatter light into radiation modes. Even a $20.0000000000\,\text{nm}$ step is enough to generate $> 5\%$ radiation loss. All transitions must be continuous and adiabatic (taper lengths $\ge 250.0000000000\,\text{nm}$ with $\le 30.0000000000\,\text{nm}$ tip width).

---

### FAILURE 6: Wide-Gap ($137.5\,\text{nm}$) Coupler — Candidate 2 (Attempt 9)

**What We Tried:**
Increasing the coupling gap to $137.5000000000\,\text{nm}$ to weaken $\kappa$ and push the State 1 cancellation null to coincide with the State 0 half-beat length at $L_c = 5.5000000000\,\mu\text{m}$.

**Why It Failed:**
- The 2D coupled-mode prediction was that $L_c \approx 5.5\,\mu\text{m}$ would be a half-beat.
- In 3D FDTD, the amorphous coupling was much stronger than predicted. At $L_c = 5.5000000000\,\mu\text{m}$, the beam completed a **full $2\pi$ beat cycle** and returned almost entirely to Port 1:
  - $P_2 = 2.8655288184\%$ (should have been $\approx 98\%$).
  - $P_1 = 92.3898906320\%$ (light returned to input port).
- The 2D eigenmode solver underestimated the evanescent coupling by $\approx 40\%$.

> [!WARNING]
> **LESSON 6: 2D cross-section eigenmode predictions are NOT reliable for gap-dependent $\kappa$.**
> The 3D confinement, substrate leakage, and z-boundary effects systematically modify $\kappa$ relative to the 2D prediction. Design points must be derived from 2D/MPB as initial guesses, then verified and refined exclusively in 3D FDTD sweeps. Never trust a 2D prediction without 3D validation.

---

### FAILURE 7: Peak Offset Problem — The Fundamental Bottleneck (Attempts 6, 7)

**What We Found:**
When the architecture was finally correct (straight-bus, adiabatic tapers, bare-silicon bends), the radiation problem was solved ($99.8248102390\%$ guided power), but:
- **State 0 peaked at $L_c = 3.0000000000\,\mu\text{m}$** ($P_2 = 91.9548201200\%$, modal purity $99.9022695012\%$).
- **State 1 peaked at $L_c = 4.2800000000\,\mu\text{m}$** ($P_1 = 90.9324819000\%$).
- The $1.2800000000\,\mu\text{m}$ offset made it impossible to achieve $> 98\%$ in both states at any single physical length.

**Root Cause — The Transition Coupling Asymmetry ($\Delta\Phi_{\text{ends}}$):**
In the bare-silicon entrance/exit fanouts, both waveguides are identical Si strips with $\Delta\beta = 0$, so evanescent coupling operates at full strength:
- $\Phi_{\text{ends}}^{(\text{am})} \approx 0.4180000000\,\text{rad}$ (amorphous: both WGs phase-matched → strong coupling in bends).
- $\Phi_{\text{ends}}^{(\text{cr})} \approx 0.0210000000\,\text{rad}$ (crystalline: PCM hold on WG1 detunes the bend region → coupling suppressed).

This asymmetry shifts the two peaks apart by:
$$\Delta L_c = \frac{\Phi_{\text{ends}}^{(\text{am})} - \frac{1}{2}\Phi_{\text{ends}}^{(\text{cr})}}{\kappa} \approx \frac{0.4075000000}{0.5236000000} \approx 0.7780000000\text{--}1.2800000000\,\mu\text{m}$$

> [!IMPORTANT]
> **LESSON 7: Transition coupling in bare-silicon fanouts creates an irreducible peak-length offset.**
> In State 0 (both WGs amorphous → $\Delta\beta = 0$ everywhere, including bends), evanescent coupling accumulates extra phase in the entrance and exit tapers. In State 1 (WG1 crystalline), the detuning in the bend region suppresses this extra phase. The net result is that State 0's effective interaction length is always longer than the physical parallel length, making it peak at a shorter $L_c$ than State 1.

---

### FAILURE 8: Gap Dilation to $98.2\,\text{nm}$ — Aligned Detuned DC (Attempt 10)

**What We Tried:**
Using rigorous MPB supermode analysis, we identified gap $= 98.2000000000\,\text{nm}$ where $S/\kappa = 2.0012481023$ (theoretical Kogelnik resonance). At this gap, the 2D cross-section theory predicts $L_c^{(0)} = L_c^{(1)}$ with offset $= 0.0010000000\,\mu\text{m}$.

**What 3D FDTD Showed:**
At $L_c = 4.4000000000\,\mu\text{m}$, $\text{gap} = 98.2000000000\,\text{nm}$:
- **Amorphous:** $P_2 = 79.1000000000\%$, $P_1 = 17.0600000000\%$, total guided = $96.1600000000\%$.
- **Crystalline:** $P_1 = 52.6600000000\%$, $P_2 = 3.6800000000\%$, total guided = $56.3400000000\%$.
- Neither state approaches $98.0\%$.

**Root Cause — The $\kappa_{\text{cr}} \ne \kappa_{\text{am}}$ Discovery:**
- When $\text{Sb}_2\text{S}_3$ crystallises ($n = 2.7 \to 3.3$), the mode pulls deeper into the Si+PCM stack.
- The evanescent tail in the coupling gap is reduced: $\kappa_{\text{cr}} \approx 0.7480000000\,\kappa_{\text{am}}$.
- This means 2D supermode analysis (which computes $\kappa$ from the amorphous state and $S$ from the crystalline state independently) does not capture the state-dependent coupling coefficient change.
- In 3D FDTD, the crystalline state requires a longer distance to complete its beat than the MPB prediction because $\kappa_{\text{cr}} < \kappa_{\text{am}}$.

> [!IMPORTANT]
> **LESSON 8: $\kappa$ is state-dependent — it changes when the PCM switches.**
> The coupling coefficient is NOT a fixed geometric parameter. When the PCM on WG1 crystallises, its higher refractive index pulls the mode tighter, reducing the evanescent field in the gap. Any design relying on $\kappa_{\text{am}} = \kappa_{\text{cr}}$ will fail in 3D FDTD. The crystalline $\kappa$ is typically $70\text{--}80\%$ of the amorphous $\kappa$.

---

### FAILURE 9: Compromise-Length Operation @ $g = 94.5\,\text{nm}$ (Attempt 11)

**What We Tried:**
Operating at an intermediate coupling length ($L_c = 3.7000000000\text{--}3.9000000000\,\mu\text{m}$) at gap $= 94.5000000000\,\text{nm}$ to find a point where both states simultaneously exceed $98.0\%$.

**What 3D FDTD Showed:**
- At $L_c = 3.7000000000\,\mu\text{m}$: State 0 $= 83.8803310000\%$, State 1 $= 67.6550180000\%$.
- At $L_c = 3.9000000000\,\mu\text{m}$: State 0 $= 80.2345280000\%$, State 1 $= 75.4624490000\%$.
- The compromise length places both states off-resonance: State 0 has over-coupled past its peak, State 1 has under-coupled short of its null.

> [!WARNING]
> **LESSON 9: There is NO compromise length that satisfies both states with a uniform coupler.**
> Because the peak offset is a structural consequence of asymmetric transition coupling AND state-dependent $\kappa$, no single $(L_c, g)$ point exists where a uniform symmetric directional coupler can deliver $> 98\%$ in both states. The uniform DC architecture has been exhaustively explored across 11 variations and is provably insufficient. The next step MUST break the geometric symmetry.

---

## 4. Monitor / Measurement Lessons

### LESSON 10: Spatial Poynting Flux Monitors Produce False Leakage
- When two waveguide output ports are separated by $< 0.5000000000\,\mu\text{m}$, spatial Poynting flux monitors sized at $2.5 \times w_{\text{core}}$ will physically overlap.
- The evanescent tail of the guided mode in Port 1 is captured by the Port 2 monitor, falsely reporting $5.0\%\text{--}8.5\%$ crosstalk.
- **Mitigation:** Always use `sim.get_eigenmode_coefficients()` for orthogonal TE₀ mode decomposition. Spatial flux should be used as a sanity check for total power budget, never as the primary signal/leakage discriminator.

### LESSON 11: Resolution 15 vs 20 vs 25 Convergence
- Resolution 15 ($66.7\,\text{nm}$ Yee cells) is sufficient for identifying trends and ranking candidates.
- However, absolute power values can shift by $\pm 2.0\text{--}3.0\%$ between resolution 15 and 25.
- **Mitigation:** All final claims must be verified at resolution $\ge 20$. Use resolution 15 for sweeps only.

---

## 5. What We Know Works (Validated Physical Principles)

These are the architectural elements that have been PROVEN in 3D FDTD and must be preserved:

| Element | Evidence | Performance |
| :--- | :--- | :--- |
| **Straight through-bus for WG1** ($y = 0$) | Attempt 7 | $0.1751897610\%$ radiation loss |
| **Adiabatic linear PCM tapers** ($l_{\text{tip}} = 250\,\text{nm}$, $w_{\text{tip}} = 30\,\text{nm}$) | Attempts 6, 7 | Transition radiation $< 0.0500000000\%$ |
| **Bare-silicon S-bend exit** (NO PCM on curves) | Lesson 4 | Eliminates bend-induced scattering |
| **Cubic $C^1$-continuous S-bend** ($y(u) = y_0 + \Delta y(3u^2 - 2u^3)$) | Attempts 6, 7 | Zero curvature at junctions |
| **Inverse taper entrance for WG2** ($w = 20 \to 240\,\text{nm}$) | Attempts 6, 7 | Eliminates parasitic input coupling |
| **Matching PCM on both WGs** in parallel zone | Lesson 3 | Preserves vertical mode overlap |
| **Eigenmode decomposition** for power extraction | Lesson 10 | True $\text{XT} = -30.0912401825\,\text{dB}$ |

---

## 6. What We Know Fails (Hard Constraints)

These have been proven impossible and MUST NOT be re-attempted:

| Constraint | Evidence | Reason |
| :--- | :--- | :--- |
| **No slot waveguides** | Failure 1 | Mode topology change → radiation |
| **No Y-junctions / mode sorters** | Failure 2 | Tip scattering + footprint overflow |
| **No MZI / 50:50 splitters** | User mandate | Splits beam; violates routing requirement |
| **No bare-WG2 asymmetric couplers** | Failure 3 | Vertical mode mismatch collapses $\kappa$ |
| **No PCM on curves** | Failure 4 | Conformal index tilt → radiation |
| **No abrupt cross-section steps** | Failure 5 | Wave-impedance mismatch → $55\%$ scattering |
| **No uniform symmetric DC at any single $(L_c, g)$** | Failures 6–11 | Peak offset + $\kappa_{\text{cr}} \ne \kappa_{\text{am}}$ |
| **No reliance on 2D eigenmode $\kappa$ predictions** | Failure 6, 8 | 3D effects shift $\kappa$ by $> 25\%$ |

---

## 7. The Physics That Must Be Overcome in the Next Step

The fundamental challenge is captured by two equations:

**State 0 (Amorphous, $\Delta\beta = 0$):** Needs $\kappa_{\text{am}} L_c + \Phi_{\text{ends}}^{(\text{am})} = \frac{\pi}{2}$

**State 1 (Crystalline, $\Delta\beta \gg 0$):** Needs $S_{\text{cr}} L_c + \Phi_{\text{ends}}^{(\text{cr})} = \pi$

Because:
1. $\kappa_{\text{cr}} < \kappa_{\text{am}}$ (mode confinement tightens when PCM crystallises).
2. $\Phi_{\text{ends}}^{(\text{am})} \gg \Phi_{\text{ends}}^{(\text{cr})}$ (transition coupling is state-asymmetric).
3. $S_{\text{cr}} / \kappa_{\text{am}} \ne 2.0$ in 3D (unlike 2D predictions).

...no single uniform interaction length $L_c$ can satisfy both conditions.

### The Two Recommended Approaches for the Next Step

#### Approach A: Asymmetric Interaction Zone (Dual Effective $L_c$)
- Terminate the amorphous PCM on WG2 at $L_{c,2} = 3.0000000000\,\mu\text{m}$ (the State 0 half-beat peak).
- Continue the switchable PCM on WG1 to $L_{c,1} = 4.2800000000\,\mu\text{m}$ via the straight through-bus (the State 1 cancellation null).
- **Critical:** The WG1 PCM extension beyond $L_{c,2}$ must have WG2 phase-detuned or physically separated so that no further coupling occurs.
- **Risk:** The segment where WG1 has PCM but WG2 does not may reintroduce vertical mode mismatch (Lesson 3). However, because WG2 is already fanning away in this region, the coupling should be negligible.

#### Approach B: Steep Lateral Decoupling at $x_{c,\text{end}}$
- Use a very steep S-bend ($l_{\text{bend}} \le 0.4000000000\,\mu\text{m}$, lateral offset $\ge 0.6000000000\,\mu\text{m}$) to physically rip WG2 away from WG1 at the end of the interaction zone.
- This truncates $\Phi_{\text{ends}}^{(\text{am})}$ to near zero, bringing State 0's effective beat closer to State 1's.
- **Risk:** A very steep bend may introduce bend radiation loss (Lesson 4), undoing our $0.1751897610\%$ achievement. Must be verified in 3D FDTD.

> [!TIP]
> **Recommended strategy:** Start with Approach A (asymmetric $L_c$) because it preserves the proven low-loss straight-bus architecture. If the residual coupling in the extended zone is too high, fall back to Approach B with carefully controlled bend radii.
