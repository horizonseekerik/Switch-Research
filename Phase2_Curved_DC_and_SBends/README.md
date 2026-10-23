# PHASE 2: CURVED DIRECTIONAL COUPLERS & S-BEND RADIATION AUTOPSY

## Overview
Phase 2 focused on sub-2 um, sub-3 um, and sub-5 um compact directional couplers utilizing Hermite cubic S-bends to separate the waveguides from the narrow coupling gap out to an 800 nm output port pitch.

## Key Outcomes & The Critical Flaw Uncovered
* **The 44% S-Bend Radiation Disaster (`test_sbend_loss.py`):**  
  Full 3D Maxwell FDTD simulations revealed that narrow 240 nm silicon strip waveguides at 1064 nm have low effective index (neff ~ 1.86). Bending these weakly-confined cores across sub-micron radii causes severe out-of-plane radiative leakage:
  - 0.75 um S-bend: **44.46% radiation loss** into cladding.
  - 1.00 um S-bend: **38.04% radiation loss** into cladding.
  - 1.50 um S-bend: **37.15% radiation loss** into cladding.
* **Conclusion:** It is physically impossible for curved S-bends on narrow silicon cores to satisfy the >= 99% total guided power specification. This finding forced the complete elimination of waveguide bends.
