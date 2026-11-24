"""Independent cross-check of the tail alpha-condition along the WHOLE tail (V4_REVIEW R-3), from the production rate tables (read only):
for l in [L_F, L*] (endpoints and 4001 log-spaced points), in arb:  (-LAMT - rho/2) kappa sqrt(l) - beta(kappa sqrt(l), l) > 0,
beta = closure.beta (BETA0, BETA1, BETAF Fn, QBL, STRIP_BETA), rho = -U_s (1 - 1e-9), kappa = a*/sqrt(L*) as in closure.tail.
The rigorous statement is the monotone split in closure.tail (key alpha_condition_whole_tail); this grid check is supporting [num-arb].
usage (env as run_v4b.sh): python3 -B t_tail_alpha.py RUNDIR TAG -> RUNDIR/t_tail_alpha_TAG.json"""
import sys, os, json
sys.dont_write_bytecode = True
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import closure as C
from flint import arb
A = lambda v: arb(repr(v))
rd, tag = sys.argv[1], sys.argv[2]; out = {}
for al in ('0.00015', '0.000125'):
    t = json.load(open(os.path.join(rd, f'rate_table_{tag}_a{al}.json')))['tail']
    Ls = A(t['L_star']); a_s = A(t['a_star']); sq = Ls.sqrt(); kappa = a_s / sq
    rho = A(-t['rate_s_U']) * (1 - arb('1e-9')); yF = A(t['sqrt_LF']); LF = yF * yF
    lhs_c = (-C.LAMT - rho / 2) * kappa
    worst = None; n = 4001
    pts = [LF, Ls] + [LF * (Ls / LF) ** (arb(i) / (n - 1)) for i in range(n)]
    for l in pts:
        m = (lhs_c * l.sqrt() - C.beta(kappa * l.sqrt(), l)) / l.sqrt()
        if worst is None or m.lower() < worst[0]: worst = (float(m.lower()), float(l.mid()))
    out[al] = dict(min_margin_per_sqrt_l=worst[0], at_l=worst[1], L_F=float(LF.mid()), L_star=t['L_star'], ok=worst[0] > 0,
                   closure_whole_tail=t.get('alpha_condition_whole_tail'))
json.dump(out, open(os.path.join(rd, f't_tail_alpha_{tag}.json'), 'w'), indent=1)
print(json.dumps(out, indent=1))
