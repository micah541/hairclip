"""[cert given step-4's JSON as read] bound of the true-curve frame shift |(x_c, ln lam) - base| <= sum_j |(M^-T)_ij| ||g_j||_rho (1+t) Ybar
(same formula as clam.py's C_alpha, rows 0,1), from the M, f_norms, tail_t stored in step4/p2/out/t9_h1_clam_s180.json.
Compares with the 1e-9 widening used by GS_tau1.py.  Read-only; writes c_frame_shift.json here."""
import json, os, re, math
from flint import arb, arb_mat, ctx
ctx.prec = 128
J = json.load(open('/home/claude/s1spec/step4/p2/out/t9_h1_clam_s180.json'))
def P(s):
    s = s.strip()
    m = re.match(r'\[\+/- ([0-9.e+-]+)\]', s)
    if m: return arb(0, float(m.group(1)) * 1.0000001)
    m = re.match(r'\[([0-9.e+-]+) \+/- ([0-9.e+-]+)\]', s)
    return arb(float(m.group(1)), float(m.group(2)) * 1.0000001 + 1e-12)
print(J['M'])
M = [[P(x) for x in row] for row in J['M']]
Mi = arb_mat([[M[i][j] for i in range(3)] for j in range(3)]).inv()
fn = [arb(math.nextafter(v, math.inf)) for v in J['f_norms']]; t = arb(math.nextafter(J['tail_t'], math.inf))
rows = [sum((abs(Mi[i, j]) * fn[j] for j in range(3)), arb(0)) * (1 + t) for i in range(3)]
up = lambda z: math.nextafter(float(z.upper()), math.inf)
out = dict(C_rows=[up(r) for r in rows], C_alpha_stored=J['C_alpha'], Ybar_box=1.7e-12, shift_at_Ybar_box=[up(r * arb('1.7e-12')) for r in rows], widening_used=1e-9,
           covered=[bool(up(r * arb('1.7e-12')) < 1e-9) for r in rows[:2]], note='M read from 5-digit ball strings (wide); row 2 must reproduce C_alpha within the string rounding')
json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'c_frame_shift.json'), 'w'), indent=1); print(json.dumps(out, indent=1))
