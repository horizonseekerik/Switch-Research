import sys, os, math
sys.path.append('/mnt/c/Users/hp/Desktop/Deepanshu/Janus Update/Switch research')
from solve_exact_coupler_parameters import solve_isolated_wg, solve_supermodes

w = 0.240
print("="*80)
print("ISOLATED WG neff AND IDEAL Lc (MPB)")
print("="*80)
for t in [0.180, 0.200, 0.220, 0.240]:
    nam = solve_isolated_wg(w, t, 2.70, resolution=32)
    ncr = solve_isolated_wg(w, t, 3.30, resolution=32)
    dn = ncr - nam
    lc = (math.sqrt(3.0) * 1.064) / (2.0 * dn)
    print(f"t_pcm={t*1000:3.0f}nm: nam={nam:.10f}, ncr={ncr:.10f}, dn={dn:.10f}, ideal_Lc={lc:.10f}um")
