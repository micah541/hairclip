"""H6 open items (analysis only).
(1) Cone K3 by covering: lemmaO_v2's conclusion int_{-1}^{1} g_xxx^2 <= C3 t per member (g = h - h0, own units). Partition [-3, 3] of cylinder c
    into pieces, each inside one member's |x| <= 1 (cone members r' in [200, 380]: region [3r'/4, 5r'/4], unit r'/4; far members L = 113,
    region [x - 113, x + 113]); a piece in member m contributes (L_c/L_m)^3 C3_m t (int (h_c''')^2 dX_c = L_c^3 int (h''')^2 dx).  Greedy covering
    by the largest available unit; S = sum (L_c/L_m)^3 C3_m t over the members used.  The inner strip r < 150 is left to (2).
(2) Inner strip K3 from step 3's time-uniform core energies E1 <= 1.81e-4, E2 <= 0.248 (omega = 28/y, D = y d/dy, windows [y, 1.5y]):
    28 int y^2 v''^2 dy <= 2 E2 + 2 E1 per window, so L_c^3 int (v'')^2 dy <= L_c^3 (E2 + E1) / (14 y_w^2) summed over the windows meeting [r_c/4, 150)."""
import json, math
C = json.load(open('/home/claude/s1spec/step4/p2/h6/h6_cases.json'))
c3c = C['cone X=3 M=0.12 K=0.25']['deviation_upper_at_needed_time']['sqrtC3t'] ** 2; c3f = C['far X=2 M=0.1 K=0.2']['deviation_upper_at_needed_time']['sqrtC3t'] ** 2
rows = []
for rc in range(200, 381, 10):
    Lc = rc / 4; lo, hi = max(rc / 4, 150.0), 1.75 * rc; x = lo; S = 0.0; used = []
    while x < min(hi, 475.0) - 1e-9:                       # cone members: largest r' with 3r'/4 <= x, i.e. r' = min(4x/3, 380)
        rp = min(4 * x / 3, 380.0); assert 3 * rp / 4 <= x + 1e-9 and rp >= 200
        S += (Lc / (rp / 4)) ** 3 * c3c; used.append(round(rp)); x = 5 * rp / 4
    while x < hi - 1e-9:                                   # far members (L = 113) beyond 475
        S += (Lc / 113) ** 3 * c3f; used.append('far'); x += 226
    E1, E2 = 1.81e-4, 0.248; strip = 0.0; y = 50.0
    while y < 150:
        if 1.5 * y > rc / 4: strip += Lc ** 3 * (E2 + E1) / (14 * y * y)
        y *= 1.5
    rows.append(dict(r_c=rc, members=used, dev_L2=math.sqrt(S), strip_L2=math.sqrt(strip)))
w = max(rows, key=lambda r: r['dev_L2'])
out = dict(args=dict(K3=0.25, C3t_cone=c3c, C3t_far=c3f, E1=1.81e-4, E2=0.248), worst_cover=w, cone_K3_margin=(0.25 - 1.21e-3) / w['dev_L2'],
           strip_max=max(r['strip_L2'] for r in rows), rows=rows)
json.dump(out, open('/home/claude/s1spec/step4/p2/h6/h6_k3.json', 'w'), indent=1)
print('cover worst', w['r_c'], w['members'], 'dev %.4f' % w['dev_L2'], 'K3 margin (excl. strip) x%.3f' % out['cone_K3_margin'], '| inner-strip L2 from E2 <= 0.248: max %.3f' % out['strip_max'])
