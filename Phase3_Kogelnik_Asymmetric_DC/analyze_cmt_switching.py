"""
Analyze P2(u) and P1(u) for delta-beta switched directional coupler
"""
import numpy as np

u = np.linspace(0, 5, 501)
arg = 0.5 * np.pi * np.sqrt(1 + u**2)
p2 = (1.0 / (1.0 + u**2)) * np.sin(arg)**2
p1 = 1.0 - p2

print(f"{'u':<8} {'arg/pi':<10} {'P1 (Bar)':<12} {'P2 (Cross)':<12} {'XT (dB)':<12} {'IL (dB)':<12}")
print("-" * 66)
for ui in [0.0, 1.0, 1.5, 1.6, 1.7, np.sqrt(3), 1.8, 1.9, 2.0, 2.5, 3.0, np.sqrt(15), 4.0]:
    a = 0.5 * np.sqrt(1 + ui**2)
    p_cross = (1.0 / (1.0 + ui**2)) * np.sin(np.pi * a)**2
    p_bar = 1.0 - p_cross
    xt = 10 * np.log10(max(p_cross / max(p_bar, 1e-9), 1e-9))
    il = -10 * np.log10(max(p_bar, 1e-9))
    print(f"{ui:<8.3f} {a:<10.3f} {p_bar*100:<11.2f}% {p_cross*100:<11.2f}% {xt:<12.2f} {il:<12.3f}")

# Find all ranges where XT <= -20 dB (P2 <= 0.01)
valid = np.where(p2 <= 0.01)[0]
print("\nRanges of u where Crosstalk <= -20 dB (P1 >= 99.0%, IL <= 0.044 dB):")
# group contiguous ranges
import itertools
from operator import itemgetter
for k, g in itertools.groupby(enumerate(valid), lambda ix: ix[0] - ix[1]):
    group = list(map(itemgetter(1), g))
    u_min = u[group[0]]
    u_max = u[group[-1]]
    print(f"  u in [{u_min:.3f}, {u_max:.3f}] -> XT <= -20 dB, IL <= 0.044 dB")
