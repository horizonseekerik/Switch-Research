import os
import sys
import math
import contextlib

@contextlib.contextmanager
def nostdout():
    with open(os.devnull, 'w') as devnull:
        old_stdout = sys.stdout
        sys.stdout = devnull
        try:
            yield
        finally:
            sys.stdout = old_stdout

from solve_switch_modes import solve_coupler_modes

configs = [
    ('Top 20nm (gap=45nm)', dict(pcm_loc='top', t_pcm=0.020, gap=0.045)),
    ('Top 60nm (gap=45nm)', dict(pcm_loc='top', t_pcm=0.060, gap=0.045)),
    ('Top 80nm (gap=45nm)', dict(pcm_loc='top', t_pcm=0.080, gap=0.045)),
    ('Top 100nm (gap=45nm)', dict(pcm_loc='top', t_pcm=0.100, gap=0.045)),
    ('Slot 45nm gap', dict(pcm_loc='slot', gap=0.045)),
    ('Slot 60nm gap', dict(pcm_loc='slot', gap=0.060)),
    ('Slot 80nm gap', dict(pcm_loc='slot', gap=0.080)),
    ('Slot 100nm gap', dict(pcm_loc='slot', gap=0.100)),
]

print("=" * 85)
print("MPB 3D VECTOR MAXWELL MODE BEAT LENGTH & COUPLING ANALYSIS")
print("=" * 85)
for name, kwargs in configs:
    with nostdout():
        am = solve_coupler_modes(pcm_state='amorphous', **kwargs)
        cr = solve_coupler_modes(pcm_state='crystalline', **kwargs)
    lc_am = am['l_c_um']
    lc_cr = cr['l_c_um']
    dn_am = am['delta_n']
    dn_cr = cr['delta_n']
    ratio = lc_am / max(lc_cr, 1e-4)
    print(f"{name:24s} | Am: Lc={lc_am:5.2f} um, dn={dn_am:6.4f} | Cr: Lc={lc_cr:5.2f} um, dn={dn_cr:6.4f} | Ratio={ratio:4.2f}")
print("=" * 85)
