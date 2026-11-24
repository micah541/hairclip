"""Error transfer for H1 over W' (T5_DESIGN l.161): with Pi the step-3 projection (IFT linearised at the W' parameter hulls),
  sqrt Lam[w] <= sqrt Lam[w~] + C_L Ybar,   |alpha - alpha~| <= C_a Ybar,
  C_L = 1.01 sqrt(0.08/mu_e)(1 + 1e-12) + K (1 + t),  K = sum_i N_i sum_j |(M^-T)_ij| ||g_j||_rho,  C_a = sum_j |(M^-T)_2j| ||g_j||_rho (1 + t).
|<e, g_j>_rho| <= (||chi32 e||_rho + ||(1 - chi32) e||_rho) ||g_j||_rho (Cauchy-Schwarz on the IBP boxes); ||chi32 e||_rho <= A_e <= Ybar;
on boxes with chi32 != 1 (x <= 34.25) |e| <= 1.25 * 2e3 Ybar / sqrt(mu_e) (T8's certified tip box), so t = 1.6031e4 (sum w_b sup(1 - chi_b)^2)^(1/2) (T7 Lemma 4b).
N_i = ||mode_i||_Lam: the certified mode norms of the 32 W' runs (h1.py, max over runs).  Factor 1.01: the graph map base -> S1 distorts e by 1 + O(|u'|^2 + |k1 u| + |k1 u|^2) <= 1 + 1e-3 (|u'| <= 2e-4, |u| <= 2e-4, k1 <= 1.55)."""
import sys, os, json, math
sys.argv = ['h1.py', 'm']; sys.path.insert(0, '/home/claude/s1spec/step4/p2/t9h')
import h1
from flint import arb, arb_mat
U = h1.U; OUT = '/home/claude/s1spec/step4/p2/out/'
R = [json.load(open(OUT + f't9_h1_p{-1 + (2 * k + 1) / 16}_coh_s180.json')) for k in range(16)]      # slab 180 (review A1)
SM = json.load(open(OUT + 'h1slab/summary.json')); PH = SM['P_hull']                             # verified parameter boxes
hull = lambda key: arb(PH[key][0]).union(arb(PH[key][1]))
xc, lam, al = hull('xc'), hull('lam'), hull('alpha')
def find(o, key):
    if isinstance(o, dict):
        for k, v in o.items():
            if k == key: return v
            r = find(v, key)
            if r is not None: return r
T8C = os.environ.get('T8CONSTS', '/home/claude/s1spec/step4/p2/t8se/constants/constants_fin.json')   # env T8CONSTS: the T8 constants of the chain
mu_e = arb(repr(find(json.load(open(T8C)), 'muE')))
g = h1.basegeo.Geo(180, 64.0, h1.TCK); S = h1.s1chart2.S1(); S.x13 = float((h1.s1chart2.T0 + 5).lower()) - 1e-6
G, M = h1.ift_pass(g, S, xc, lam, al)
Mi = arb_mat([[M[i][j] for i in range(3)] for j in range(3)]).inv()           # M^-T
fn2 = [arb(0)] * 3; tw = arb(0)
for b in h1.IB.boxes():
    for j in range(3):
        v = arb(U(b['g'][j + 1][0])); fn2[j] += b['w'] * v * v
    c_ = arb(U(1 - b['chi'][0])); tw += b['w'] * c_ * c_                    # ||(1 - chi32) e||_rho^2 <= E_tip^2 sum w (1 - chi)^2
fn = [arb(U(v.sqrt())) for v in fn2]
Nn = {n: max(r['Lambda']['norms'][n] for r in R) for n in ('tr', 'sc', 'al')}
Ni = [arb(Nn['tr']) / lam, arb(Nn['sc']), arb(Nn['al'])]
t = arb('1.6031e4') * tw.sqrt()                       # T7 Lemma 4b: sup|e| <= 1.6031e4 Lambda_e^(1/2) on [0, 34.25] (review A3)
K = sum((Ni[i] * sum((abs(Mi[i, j]) * fn[j] for j in range(3)), arb(0)) for i in range(3)), arb(0))
CL = arb('1.01') * (arb('0.08') / mu_e).sqrt() * (1 + arb('1e-12')) + K * (1 + t)
Ca = sum((abs(Mi[2, j]) * fn[j] for j in range(3)), arb(0)) * (1 + t)
base = SM['sup_sqrtLam_base']                          # slab 180: 16 Krawczyk points + 32 derivative balls
out = dict(args=dict(h1=[f't9_h1_p{-1 + (2 * k + 1) / 16}_coh_s180.json' for k in range(16)], summary='h1slab/summary.json', tip_const=1.6031e4,
                     t8consts=T8C, env={k: os.environ[k] for k in ('COHERENT', 'S1PROF', 'S1TT', 'T8CONSTS') if k in os.environ}), mu_e=float(mu_e.mid()), M=[[M[i][j].str(5) for j in range(3)] for i in range(3)], f_norms=[U(v) for v in fn], tail_t=U(t),
           K=U(K), C_Lambda=U(CL), C_alpha=U(Ca), sqrtLam_base_Wp=base, err_at_1e12=U(CL * arb('1e-12')),
           sqrtLam_total_Wp=math.nextafter(base + U(CL * arb('1e-12')), math.inf), Ybar_limit_H1=math.nextafter(float(((arb('5e-8') - base) / CL).lower()), -math.inf))   # a limit: lower end (CL is a wide ball)
json.dump(out, open(OUT + 't9_h1_clam_s180.json', 'w'), indent=1); print(json.dumps(out))
