"""H6 analysis (review A1 follow-up; analysis only, step 3 unchanged): Lemma O constants for the candidate families, from
writeup/lemmaO/lemmaO_v3_interval.case (read-only import), at the needed time ratios t1/L^4 of step 3.
cone: X = 3, M = 0.12, K = K2 = K3 in (0.18, 0.25, 0.30), need 1.6e-7.  far: X = 2, need 6.13e-9 (L = 113), (M, K) data-driven candidates."""
import sys, json
sys.path.insert(0, '/home/claude/s1spec/writeup/lemmaO')
from fractions import Fraction as F
import lemmaO_v3_interval as L
out = {}
for K in (F(18, 100), F(25, 100)):
    M = F(12, 100); c = L.case(f'cone X=3 M=0.12 K={float(K)}', F(3), M, K, K, (K - M) / M, 1.6e-7, grid=14)
    out[c['name']] = {k: c[k] for k in ('C2_upper', 't_star_lower', 'margin_factor_lower', 'deviation_upper_at_needed_time')}
for M, K in ((F(1, 10), F(15, 100)), (F(1, 10), F(2, 10)), (F(2, 10), F(3, 10))):
    c = L.case(f'far X=2 M={float(M)} K={float(K)}', F(2), M, K, K, (K - M) / M, 6.133187276793765e-09, grid=14)
    out[c['name']] = {k: c[k] for k in ('C2_upper', 't_star_lower', 'margin_factor_lower', 'deviation_upper_at_needed_time')}
json.dump(out, open('/home/claude/s1spec/step4/p2/h6/h6_cases.json', 'w'), indent=1)
for k, v in out.items(): print(k, {a: v[a] for a in ('t_star_lower', 'margin_factor_lower')}, v['deviation_upper_at_needed_time'])
