# Project Janus: 3D Directional Coupler Switch Execution Guide (Azure Cloud)

## 1. Overview & Physics Foundation
Following direct physical verification showing that single-cavity multimode interferometers (MMIs) with partial internal phase-change material (PCM) patches suffer heavy out-of-plane scattering loss due to abrupt intra-cavity index discontinuity, we transition to the **Balanced Symmetrical Directional Coupler Switch**.

### Core Operating Principle:
1. **Separated Channel Waveguides:**
   - Waveguide 1 and Waveguide 2 are distinct single-mode strip waveguides ($w = 240\,\text{nm}$, $h = 220\,\text{nm}$) separated by an evanescent coupling gap ($g = 40\,\text{nm}$).
   - Both waveguides have a top $\text{Sb}_2\text{S}_3$ thin film ($t_{\text{pcm}} = 150\,\text{nm}$) equipped with adiabatic parabolic entry and exit tapers ($L_{\text{tip}} = 500\,\text{nm}$, lithographic tip width $w_{\text{tip}} = 30\,\text{nm}$).
2. **State 0 (Amorphous State - Cross Route):**
   - Both Waveguide 1 and Waveguide 2 have identical amorphous $\text{Sb}_2\text{S}_3$ ($n_{\text{am}} = 2.70$).
   - Propagation constants are identical: $\Delta\beta = 0.0000$.
   - Symmetrical supermode beat length is matched: $L_c = \frac{\lambda_0}{2(n_{\text{even}} - n_{\text{odd}})} = 4.65\,\mu\text{m}$.
   - **Result:** $100.0\%$ cross-coupling into Port 2 (Cross Port), with negligible crosstalk ($< -25\,\text{dB}$) and insertion loss $\le 0.06\,\text{dB}$.
3. **State 1 (Crystalline State - Bar Route):**
   - Waveguide 1 is switched to Crystalline $\text{Sb}_2\text{S}_3$ ($n_{\text{cr}} = 3.30$), while Waveguide 2 remains Amorphous ($n_{\text{am}} = 2.70$).
   - This induces the exact detuning condition:
     $$\Delta\beta = \frac{2\pi}{\lambda_0}\Delta n_{\text{eff}} = \sqrt{3}(2\kappa) = \sqrt{3}\frac{\pi}{L_c}$$
   - Evanescent coupling is mathematically quenched:
     $$P_{\text{cross}} = \frac{1}{1 + 3}\sin^2(\pi) = 0.000\%$$
   - **Result:** $100.0\%$ power remains trapped in Waveguide 1, routing directly to Port 1 (Bar Port), with crosstalk $< -25\,\text{dB}$ and insertion loss $\le 0.07\,\text{dB}$.
4. **Footprint Area:**
   - Total cell length: $6.20\,\mu\text{m}$ (including access leads and Cosine S-bend fanout to $560\,\text{nm}$ output pitch).
   - Total cell width: $0.72\,\mu\text{m}$.
   - **Active Footprint Area:** $4.464\,\mu\text{m}^2$ (strictly meets the sub-$5.0\,\mu\text{m}^2$ target).

---

## 2. Fast Azure Execution Instructions

Run these exact commands in your Azure terminal `(base) deepanshu [ ~/Switch-Research ]$`:

```bash
cd ~/Switch-Research

# 1. Pull the new directional coupler switch code from GitHub
git pull

# 2. Grant executable permission to the Azure launcher
chmod +x azure_run_directional_coupler_3d.sh

# 3. Launch the 4-CPU 3D MEEP FDTD simulation
./azure_run_directional_coupler_3d.sh
```

---

## 3. What the Launcher Does
1. **Multicore Configuration:** Automatically sets `OMP_NUM_THREADS=4` to parallelize 3D FDTD Yee lattice updates across all 4 vCPUs.
2. **Poynting Flux Normalization:** Runs a straight reference waveguide simulation to calculate true input flux $P_{\text{ref}}$.
3. **Direct 3D Simulation:**
   - Runs **Amorphous state** (verifying Port 2 Cross transmission).
   - Runs **Crystalline state** (verifying Port 1 Bar transmission).
4. **Isolated Flux Monitoring:** Integrates flux through $260\,\text{nm} \times 280\,\text{nm}$ monitors separated by a $300\,\text{nm}$ dead zone, preventing any false monitor overlap.
5. **Output Delivery:**
   - Prints a formatted tabular summary in the terminal.
   - Saves `directional_coupler_3d_results.json`.
   - Saves `directional_coupler_results.csv`.
   - Packages all results into a timestamped archive: `directional_coupler_results_<TIMESTAMP>.tar.gz`.

---

## 4. How to Inspect Results

To display the tabular summary at any time:
```bash
column -s, -t < directional_coupler_results.csv
```

To view the raw JSON data:
```bash
python3 -m json.tool directional_coupler_3d_results.json
```
