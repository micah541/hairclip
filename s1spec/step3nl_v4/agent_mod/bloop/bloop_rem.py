"""alpha^2 remainders of E_13, E_23 (|alpha| <= 1e-3), L1 hull integrals with |f3| (hull) :
 |E13 - alpha e13| <= alpha^2 r1,  r1 = int rho|f3| ( |h2 - q h2'| + |q k| (|psi'|+AM|h2'|)(|psi|+AM|h2|)/alpha_g )
 |E23 - alpha e23| <= alpha^2 r2,  r2 = int rho|f3| ( |h2' cos| + |k| (|psi'|+AM|h2'|)(|psi|+AM|h2|)/alpha_g ) (|cos|<=1).
 Also |E33| <= alpha * 2 int rho |chi h2| |f3| (crude check)."""
import sys, json
from flint import arb
sys.argv_saved = sys.argv; sys.argv = sys.argv[:1]
sys.path.insert(0, '/home/claude/s1spec/step3nl/agent_mod/bloop')
import mg_local as MG
sys.argv = sys.argv_saved
HR = MG.HR; import gen_hs3 as G3
BUMP = (0.8324 * 9.0923e-8 + 0.4245 * 2.5707e-7) * 1.01
AM = 1e-3; AG = 0.98305
ab = lambda z: arb(abs(z).upper())
lines = [l for l in open('/home/claude/s1spec/step3/outA5.txt') if l[:2] in ("MB", "MC")]
r1 = arb(0); r2 = arb(0)
for b, line in zip(MG.boxes(), lines):
    v = list(map(float, line.split()[1:])); tag = line[:2]
    h = [MG.A(v[2 + 2 * i], v[3 + 2 * i]) for i in range((len(v) - 2) // 2)]
    ws0 = h[5] if tag == "MB" else h[9]
    f3 = ab(b['chi'] * ws0) + (arb(BUMP) if (tag == 'MC' and b['ta'] >= G3.T0 - G3.LB - 1e-12) else arb(0))
    X = b['X']; xc = float(X.mid()); y = HR.jets_any(HR.F, X, xc)[:2]; hh = HR.jets_any(HR.Hf, X, xc)[:2]
    q, k = b['q'], b['k']
    quad = ab(q) * ab(k) * (ab(y[1]) + AM * ab(hh[1])) * (ab(y[0]) + AM * ab(hh[0])) / AG
    r1 += b['w'] * f3 * (ab(hh[0] - q * hh[1]) + quad)
    r2 += b['w'] * f3 * (ab(hh[1]) + ab(k) * (ab(y[1]) + AM * ab(hh[1])) * (ab(y[0]) + AM * ab(hh[0])) / AG)
print("r1", float(r1.upper()), "r2", float(r2.upper()))
json.dump({"r1": float(r1.upper()), "r2": float(r2.upper()), "alpha_max": AM, "alpha_g_min": AG}, open('rem.json', 'w'))
