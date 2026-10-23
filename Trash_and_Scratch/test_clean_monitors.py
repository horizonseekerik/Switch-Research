"""
Test with non-overlapping monitors
"""
import sys, math
import numpy as np
import meep as mp
from calibrate_coupler_3d import monitor_coupler_evolution

# Modify calibrate_coupler_3d to use non-overlapping monitors:
# mon_w = 0.24 um, perfectly centered on each waveguide
