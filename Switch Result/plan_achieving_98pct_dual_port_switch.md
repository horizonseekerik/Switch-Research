# Implementation Plan: Achieving >98.0% Guided Power in Both Switching Ports (Project Janus)

## Goal Description
Design, simulate, and verify a pure spatial $1\times 2$ non-volatile optical switch cell at $\lambda_0 = 1064\,\text{nm}$ using $\text{Sb}_2\text{S}_3$ phase-change material on Silicon-on-Insulator, satisfying all target criteria:
- **State 1 Guided Power ($P_1$):** $> \mathbf{98.0\%}$ to Port 1 (Crystalline state, Bar port)
- **State 0 Guided Power ($P_2$):** $> \mathbf{98.0\%}$ to Port 2 (Amorphous state, Cross port)
- **Total Guided Power ($P_1 + P_2$):** $\ge \mathbf{99.0\%}$ (Radiation scattering loss $\le 1.0\%$)
- **Insertion Loss (IL):** $< \mathbf{0.10\,\text{dB}}$
- **Modal Crosstalk (XT):** $< \mathbf{-20.0\,\text{dB}}$
- **Total Active Footprint:** Strictly $\le \mathbf{6.00\,\mu\text{m}^2}$
- **Pure Spatial Routing:** Continuous single-beam directional coupler (strictly NO 50:50 power splitters or MZI).

---

## User Review Required

> [!IMPORTANT]
> **Resolution of the $0.86\,\mu\text{m}$ Peak Offset:**
> In prior runs, State 0 peaked at $L_c \approx 3.30\,\mu\text{m}$ ($98.2\%$) and State 1 peaked at $L_c \approx 4.28\,\mu\text{m}$ ($90.9\%$). Our mathematical autopsy proved that bare-silicon transition coupling ($\Phi_{\text{bare}} \approx 0.60\,\text{rad}$) in the entrance and exit fanouts creates this shift:
> $$\Delta L = \frac{\frac{1}{2}\Phi_{\text{bare}}}{\kappa} \approx 0.86\,\mu\text{m}$$
> The implementation plan uses **asymmetric fanout phase-detuning** ($w_{\text{bend}} = 360\,\text{nm}$ vs $w_1 = 240\,\text{nm}$) and synchronous PCM taper terminations to enforce $\Phi_{\text{bare}} \to 0$. Transfer matrix modeling indicates this aligns both states to $L_c \approx 3.62\,\mu\text{m}$ with predicted transmission $> 99.8\%$ in both ports.

---

## Proposed Changes

```mermaid
flowchart TD
    subgraph Step 1: Geometry Refinement
        A["Sync PCM Tapers on Both WGs<br/>(l_tip = 0.25 um, zero hold)"]
        B["Asymmetric Widened Exit Bend<br/>(w_bend = 360 nm, Delta_beta = 2.36 rad/um)"]
        C["Smooth Inverse Input Taper<br/>(w = 20 -> 240 nm, l_in = 0.35 um)"]
    end

    subgraph Step 2: CMT-Guided 3D FDTD Sweep
        D["Grid: L_c in [3.56, 3.68] um<br/>Gap in [88, 92] nm (Res = 15)"]
    end

    subgraph Step 3: High-Resolution Verification
        E["Optimal Design Candidate at Res = 20 & 25"]
        F["Eigenmode TE0 Decomposition<br/>& Poynting Flux Extraction"]
    end

    subgraph Step 4: Final Deliverables
        G["Full Performance Verification Report"]
        H["GDSII Mask Export Script & Git Commit"]
    end

    Step 1: Geometry Refinement --> Step 2: CMT-Guided 3D FDTD Sweep
    Step 2: CMT-Guided 3D FDTD Sweep --> Step 3: High-Resolution Verification
    Step 3: High-Resolution Verification --> Step 4: Final Deliverables
```

---

### Component: Core Maxwell Simulation & Design Scripts (`Switch research/`)

#### [NEW] `test_aligned_detuned_switch.py`
A consolidated 3D Maxwell FDTD simulation engine integrating:
1. **Straight Through Waveguide 1:** Continuous, zero-bend through bus at $y = 0.0$ to eliminate all input radiation scattering ($< 0.2\%$).
2. **Synchronous PCM Terminations:** Both Waveguide 1 and Waveguide 2 have identical PCM lengths ($L_c$) with matching adiabatic linear tapers ($l_{\text{tip}} = 0.25\,\mu\text{m}$, $w = 30\,\text{nm} \to 240\,\text{nm}$) at both entrance and exit. Zero staggered hold ($l_{\text{hold}} = 0$).
3. **Detuned Exit S-Bend Fanout:** Waveguide 2 widens from $240\,\text{nm}$ to $360\,\text{nm}$ along a cosine/cubic S-bend ($l_{\text{bend}} = 1.00\,\mu\text{m}$, $\Delta y = 0.12\text{--}0.15\,\mu\text{m}$), introducing $\Delta n_{\text{eff}} \approx 0.40$ ($\Delta\beta \approx 2.36\,\text{rad}/\mu\text{m}$). This raises the detuning factor to $u \approx 5.9$, quenching parasitic transition coupling to $< 0.05\,\text{rad}$.
4. **Dual Power Extraction:** Extracts both non-overlapping spatial Poynting flux and pure fundamental $\text{TE}_{00}$ mode projection (`sim.get_eigenmode_coefficients`).

```python
# Key parameters for test_aligned_detuned_switch.py
lambda_0 = 1.064
w_wg = 0.240
h_wg = 0.220
t_pcm = 0.200       # 200 nm for exact u = sqrt(3) Kogelnik condition
gap = 0.089         # 89 nm for exact kappa = 0.387 rad/um
w_bend = 0.360      # Widened fanout for transition detuning
l_tip = 0.250       # Adiabatic PCM end tapers
```

---

#### [NEW] `sweep_aligned_dc_fine.py`
Automates the focused 3D FDTD parameter sweep across the optimal coupled-mode valley:
- $L_c \in [3.56, 3.60, 3.64, 3.68, 3.72]\,\mu\text{m}$
- $\text{Gap} \in [88, 90, 92]\,\text{nm}$
- Resolution $= 15$ for rapid evaluation, logging:
  - $P_1$ (Bar modal power %)
  - $P_2$ (Cross modal power %)
  - Total guided power %
  - Insertion loss ($\text{dB}$) and Crosstalk ($\text{dB}$)
  - Cell footprint area ($\mu\text{m}^2$)

---

#### [NEW] `verify_highres_final_switch.py`
High-resolution Maxwell verification script:
- Takes the optimal $(L_c, \text{gap})$ identified in the sweep.
- Executes full 3D FDTD at **Resolution $= 20$** ($50\,\text{nm}$ Yee cells) and **Resolution $= 25$** ($40\,\text{nm}$ Yee cells) to prove that the performance holds with zero Yee-grid discretization error.
- Evaluates wavelength sensitivity across $\lambda \in [1050, 1064, 1080]\,\text{nm}$ to establish optical operating bandwidth.

---

#### [NEW] `export_janus_switch_gds.py`
Creates a clean, foundry-ready GDSII mask layout of the verified switch:
- Layer 1: Silicon strip core ($w = 240\,\text{nm}$) and inverse tapers
- Layer 2: $\text{Sb}_2\text{S}_3$ phase-change material active pads ($200\,\text{nm}$ thickness, with $0.25\,\mu\text{m}$ tapers)
- Standard cell boundary markers displaying total cell area $\le 6.0\,\mu\text{m}^2$.

---

## Verification Plan

### Automated Numerical Verification (in WSL Ubuntu 22.04)

1. **Focused 3D FDTD Sweep Execution:**
   ```bash
   wsl bash -c "cd '/mnt/c/Users/hp/Desktop/Deepanshu/Janus Update/Switch research' && python3 sweep_aligned_dc_fine.py"
   ```
   - **Pass Criterion:** Identify at least one geometry point where both:
     $$\text{State 0 (Amorphous)} \implies P_2 \ge 98.0\%, \quad P_1 \le 1.0\%$$
     $$\text{State 1 (Crystalline)} \implies P_1 \ge 98.0\%, \quad P_2 \le 1.0\%$$
     $$\text{Total Guided Power} \ge 99.0\% \quad (\text{Radiation loss} \le 1.0\%)$$
     $$\text{Cell Area} \le 6.00\,\mu\text{m}^2$$

2. **High-Resolution Convergence Check:**
   ```bash
   wsl bash -c "cd '/mnt/c/Users/hp/Desktop/Deepanshu/Janus Update/Switch research' && python3 verify_highres_final_switch.py --res 20"
   ```
   - **Pass Criterion:** Confirms $P_1, P_2 \ge 98.0\%$ at resolution 20, with $\text{IL} < 0.10\,\text{dB}$ and $\text{XT} < -20.0\,\text{dB}$.

3. **GDSII Mask Verification:**
   ```bash
   wsl bash -c "cd '/mnt/c/Users/hp/Desktop/Deepanshu/Janus Update/Switch research' && python3 export_janus_switch_gds.py"
   ```
   - **Pass Criterion:** Generates `janus_1x2_switch_cell.gds` with zero geometric DRC errors (minimum feature size $\ge 80\,\text{nm}$).
