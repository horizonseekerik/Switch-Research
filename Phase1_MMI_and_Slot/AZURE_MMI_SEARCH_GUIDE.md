# Project JANUS: Single-Cavity 1x2 MMI Switch Cloud Search Guide
**Target Objective:** Multi-Core Full-Wave 3D Vector Maxwell MEEP FDTD Search for Ultra-Compact Single-Cavity 1x2 MMI Optical Switch.  
**Hardware Platform:** Standard 4-vCPU Azure Cloud VM (`Standard_D4s_v5` or similar free-tier / burstable instance).  
**Operating Wavelength:** $\lambda_0 = 1064.0\,\text{nm}$ ($281.76\,\text{THz}$).  

---

## 1. Design Architecture & Specifications

### 1.1 Structural Topology
- **Device Type:** Single-Cavity 1x2 Multimode Interference (MMI) Optical Switch.
- **Input:** Single optical input port launching the fundamental $\text{TE}_0$ eigenmode into the multimode cavity.
- **Multimode Cavity:** Silicon-on-Insulator (SOI) platform ($220\,\text{nm}$ crystalline silicon core, $n_{\text{Si}} = 3.48$, $\text{SiO}_2$ substrate and cladding, $n_{\text{clad}} = 1.449$).
- **Active Phase Modulation:** Asymmetric $\text{Sb}_2\text{S}_3$ (antimony trisulfide) thin film atop the multimode cavity ($t_{\text{pcm}} = 80\text{--}150\,\text{nm}$).
  - Amorphous State: $n = 2.70$, $k \le 1.0\times 10^{-5}$
  - Crystalline State: $n = 3.30$, $k \le 8.0\times 10^{-5}$
- **Outputs:** Two fully isolated output ports:
  - Port 1 (Bar) at $+y_{\text{port}}$
  - Port 2 (Cross) at $-y_{\text{port}}$

### 1.2 Target Constraints
1. **Single-Port Routing:** Light must route into exactly one designated port per state:
   - State 0 (Amorphous) $\implies$ Port 2 (Cross)
   - State 1 (Crystalline) $\implies$ Port 1 (Bar)
2. **Crosstalk & Isolation:** Unselected port leakage must be **$< 1.0\%$** ($\text{XT} \le -20.0\,\text{dB}$, Extinction Ratio $\ge 20\,\text{dB}$) for **BOTH** states.
3. **Insertion Loss:** Total insertion loss must be **$\le 0.08\,\text{dB}$** ($\text{Transmission} \ge 98.18\%$) for **BOTH** states.
4. **Footprint Area:** Aggressively minimized ($< 5.0\,\mu\text{m}^2$, prioritizing Sub-$2.0\,\mu\text{m}^2$).
5. **Simulation Rigor:** 100% direct 3D Poynting flux integration in MIT MEEP FDTD.

---

## 2. Cloud Execution Instructions (Azure VM / Remote Linux)

### Step 1: Connect to Cloud VM & Pull Updates
```bash
cd ~/Switch-Research
git pull
```

### Step 2: Make the Launcher Executable
```bash
chmod +x azure_run_1x2_mmi_search.sh
```

### Step 3: Run the Master Parallel Search
```bash
./azure_run_1x2_mmi_search.sh
```

---

## 3. Advanced CLI Options for `master_1x2_mmi_cloud_search.py`

You can also run the Python master script directly with custom arguments:

| Flag | Default | Description |
| :--- | :---: | :--- |
| `--cpus` | Auto (`nproc`) | Number of parallel worker processes across CPU cores |
| `--resolution` | `20` | MEEP spatial resolution in pixels/$\mu\text{m}$ ($50\,\text{nm}$ grid) |
| `--mode` | `both` | Execution mode: `screen` (fast analytical), `full3d` (3D FDTD), `both` |
| `--quick` | Off | Runs a small test subset for quick smoke verification |
| `--output-dir` | `.` | Target directory for all generated datasets and summaries |

### Examples:
```bash
# Fast 3-minute smoke test across 4 cores
python3 master_1x2_mmi_cloud_search.py --cpus 4 --quick --mode both

# Full comprehensive search (72 candidates evaluated in parallel)
python3 master_1x2_mmi_cloud_search.py --cpus 4 --mode both
```

---

## 4. Output Data & Result Verification

The master script automatically outputs:
1. **`single_cavity_1x2_mmi_results.csv`**: Full tabular CSV containing all evaluated geometries, transmission percentages, insertion losses, and crosstalk values.
2. **`top_1x2_mmi_candidates.json`**: Filtered list of winning configurations that strictly satisfy $\text{IL} \le 0.08\,\text{dB}$ and $\text{XT} \le -20.0\,\text{dB}$.
3. **`mmi_1x2_search_results_<TIMESTAMP>.tar.gz`**: Standalone compressed archive of the entire run for downloading via `scp` or Azure Storage.
