# Switch Research Repository: Project JANUS

This directory contains research, specifications, simulation models, and Azure 4-CPU search engines for optimizing and shrinking the optical switches utilized in **Project JANUS (Spatial Optical RNS Photonic AI Computing Architecture)**.

---

## 🎯 Primary Research Goals
1. **Aggressive Footprint Shrinkage:** Shrink the unit cell footprint from **$12.04\ \mu\text{m}^2$** to **$< 5.00\ \mu\text{m}^2$** ($> 58\%$ area reduction).
2. **Fabric-Wide Scaling:** Slash total on-chip switch fabric for 3,932,160 switches from **$47.34\text{ mm}^2$ down to $< 18.0\text{ mm}^2$**, freeing $> 29\text{ mm}^2$ on the $100\text{ mm}^2$ die.
3. **Insertion Loss (IL) Parity or Improvement:** Maintain or lower insertion loss ($\text{IL} \le 0.15\text{ dB}$ Amorphous Cross, $\le 0.20\text{ dB}$ Crystalline Bar).
4. **Crosstalk Immunity:** Maintain optical crosstalk $\le -20.0\text{ dB}$ (Extinction Ratio $\ge 20.0\text{ dB}$, $\text{SCR} > 150:1$).

---

## 📑 Directory Contents & Roadmap

### 🔬 1. Sub-5 $\mu\text{m}^2$ Switch Simulation & Azure 4-CPU Search Engine
- **[`sub5um_switch_meep.py`](./sub5um_switch_meep.py)**: Full-wave 2D MEEP FDTD and Coupled-Mode Theory (CMT) simulation class modeling the ultra-compact switch cell with integrated Hermite S-bend tapers and apodized $\text{Sb}_2\text{S}_3$ patches.
- **[`sub5um_switch_3d_meep.py`](./sub5um_switch_3d_meep.py)**: Full-wave 3D Vector Maxwell MEEP FDTD simulation model capturing realistic vertical layer stack ($220\,\text{nm}$ Si core, $20\,\text{nm}$ $\text{Sb}_2\text{S}_3$ film, $\text{SiO}_2$ BOX substrate & cladding).
- **[`mpb_coupler_analyzer.py`](./mpb_coupler_analyzer.py)**: Maxwell eigenmode supermode solver computing $n_{\text{even}}, n_{\text{odd}}, \Delta n, L_c, \kappa$, and the exact $\pi$-null decoupling condition $\Delta\beta L_c = \sqrt{3}\pi$.
- **[`azure_4cpu_grid_search.py`](./azure_4cpu_grid_search.py)**: Multi-process grid-search engine optimized for a 4-vCPU Azure Cloud node. Evaluates 5,040 design points across gaps, widths, lengths, and fanouts.
- **[`azure_run_switch_search.sh`](./azure_run_switch_search.sh)**: Executable bash script for running the search, executing the test suite, and archiving results on Azure.
- **[`AZURE_4CPU_SEARCH_GUIDE.md`](./AZURE_4CPU_SEARCH_GUIDE.md)**: Complete step-by-step instructions for provisioning a 4-vCPU Azure VM (`Standard_D4s_v5` or `Standard_F4s_v2`, ~$0.19/hr) and executing the search.
- **[`test_sub5um_switch.py`](./test_sub5um_switch.py)**: Automated pytest suite verifying geometric footprint ($< 5.0\ \mu\text{m}^2$), 3D cell consistency, insertion loss, crosstalk, energy conservation, supermode eigensolver, and MEEP FDTD smoke execution (**8/8 Tests PASSED**).

### 📊 2. Generated Optimization Datasets & Results
- **[`sub5um_switch_search_results.csv`](./sub5um_switch_search_results.csv)**: Complete data table of evaluated parameter combinations with losses, crosstalk, and areas.
- **[`sub5um_top_candidates.json`](./sub5um_top_candidates.json)**: Top candidate switch configurations achieving $\text{Area} < 5.0\ \mu\text{m}^2$ and sub-$0.15\text{ dB}$ insertion loss.

### 🏛️ 3. Baseline Switch Data & References
- **[`EXISTING_JANUS_SWITCH_DATA.md`](./EXISTING_JANUS_SWITCH_DATA.md)**: Comprehensive master dossier documenting the existing $12.04\ \mu\text{m}^2$ switch.
- **[`switch_spec_data.json`](./switch_spec_data.json)**: Machine-readable JSON specification registry of the baseline switch.
- **[`pcm_materials_benchmark.csv`](./pcm_materials_benchmark.csv)**: Solid-state physics comparison of $\text{Sb}_2\text{S}_3$ vs. GST-225, GST-467, and GSST.

---

## 🏆 Current Leading Candidate Design (Sub-5 $\mu\text{m}^2$ Winner)

| Metric | Baseline JANUS Switch | **Optimized Sub-5 $\mathbf{\mu m^2}$ Switch** | Progress / Margin |
| :--- | :---: | :---: | :---: |
| **Unit Cell Length** | $8.60\ \mu\text{m}$ | **$4.05\ \mu\text{m}$** | **$52.9\%$ shorter** |
| **Unit Cell Width** | $1.40\ \mu\text{m}$ | **$1.00\ \mu\text{m}$** | **$28.6\%$ narrower** |
| **Unit Cell Footprint Area** | $12.04\ \mu\text{m}^2$ | **$4.05\ \mu\text{m}^2$** | **$66.4\%$ smaller ($< 5.0\ \mu\text{m}^2$ Target MET)** |
| **3.93M Switches Fabric Area** | $47.34\text{ mm}^2$ | **$15.93\text{ mm}^2$** | **$+31.41\text{ mm}^2$ surplus chip area saved!** |
| **Amorphous Insertion Loss** | $0.057\text{ dB}$ | **$0.069\text{ dB}$** | **Target Met ($\le 0.15\text{ dB}$)** |
| **Crystalline Insertion Loss** | $0.142\text{ dB}$ | **$0.070\text{ dB}$** | **$51\%$ loss reduction ($0.072\text{ dB}$ lower!)** |
| **Worst-Case Optical Crosstalk**| $-21.86\text{ dB}$ | **$-22.67\text{ dB}$** | **$+0.81\text{ dB}$ higher extinction ratio** |
| **Signal-to-Crosstalk Ratio (SCR)**| $163.6:1$ | **$184.9:1$** | **Zero false-detection margin enhanced** |

---
*Maintained under Project JANUS Core Photonic Architecture.*
