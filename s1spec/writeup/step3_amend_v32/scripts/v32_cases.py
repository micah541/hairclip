"""Step 3 amendment v3.2 (draft): Lemma O constants for the cone family (unchanged), the v3.1 far family and the v3.2 far family.
Read-only import of writeup/lemmaO/lemmaO_v3_interval.case (python-flint arb, 200 bits; outputs rounded outward at 1e-6).
Needed time ratios carry the horizon factor 1/(1 - 4 abar), abar <= 2.2e-7:  cone 50^-4 -> 1.6001e-7,  far 113^-4 -> 6.134e-9."""
import sys, json
sys.dont_write_bytecode = True
sys.path.insert(0, '/home/claude/s1spec/writeup/lemmaO')
from fractions import Fraction as F
import lemmaO_v3_interval as L
KEYS = ('w', 'C1_upper', 'C2_upper', 'C3_upper', 't1_star_lower', 't2_star_lower', 't3_star_lower', 't_star_lower', 'needed',
        'margin_factor_lower', 'deviation_upper_at_needed_time')
out = {}
def run(tag, X, M, K, need):
    c = L.case(tag, F(X), M, K, K, (K - M) / M, need, grid=14)
    out[tag] = {k: c[k] for k in KEYS}
    d = c['deviation_upper_at_needed_time']
    print(f"{tag}: t* >= {c['t_star_lower']:.4g} need {need:.4g} margin x{c['margin_factor_lower']:.3f} | C0,C1,C2 sqrt t = "
          f"{d['C0']:.4g}, {d['C1']:.4g}, {d['C2']:.4g}, sqrt(C3 t) = {d['sqrtC3t']:.4g}", flush=True)
for K in (F(18, 100),):
    run(f'cone X=3 M=0.12 K={float(K)}', 3, F(12, 100), K, 1.6001e-7)
for M, K in ((F(1), F(3, 2)), (F(15, 100), F(3, 10))):
    run(f'far X=2 M={float(M)} K={float(K)}', 2, M, K, 6.134e-9)
json.dump(out, open('/home/claude/s1spec/writeup/step3_amend_v32/scripts/v32_cases.json', 'w'), indent=1)
