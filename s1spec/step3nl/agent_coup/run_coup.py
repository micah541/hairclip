"""Rigorous coupling constants Q1-Q6 for one P.   usage: python3 run_coup.py PX [out.json] [--fast]"""
import sys, json, time, math
from functools import lru_cache
from flint import arb, ctx
import coup_common as C
from coup_common import *

PXs = sys.argv[1]; OUT = sys.argv[2] if len(sys.argv) > 2 and not sys.argv[2].startswith('--') else f"coup_P{PXs}.json"
FAST = '--fast' in sys.argv
T0 = time.time(); LOG = []
def log(*a):
    s = " ".join(str(x) for x in a); print(s, flush=True); LOG.append(s)
Wt = Weight(PXs); log(f"P={PXs}: phi_I table built {time.time()-T0:.1f}s")
rel = 5e-3 if FAST else 1e-3
n_strip, n_ramp, n_out = (128, 64, 32) if FAST else (512, 256, 128)
ABf7, ABf5, CF5, CF7 = get_AB(7), get_AB(5), get_coef(5), get_coef(7)
NP = {}
def rec(name, res):
    up_, lo_, n = res; NP[name] = dict(n_pieces=n, lower_witness=lo_, upper=up_)
    return up_

# ------------------------------------------------------------------ strip (sigma in [29,32], variable sigma, weight rho, cut chi)
@lru_cache(maxsize=None)
def strip_cell(a, b):
    SIG, X = Iv(a, b), None
    X = SIG - TIP
    c = cut_chi(SIG); phi, P = rho_data(X); k, m, p, q = kmpq(SIG)
    Ap, Bp = ABf5(*c, *P, k, m, p, q); coefs = CF5(*c, *P, k, m, p, q)
    rho = (-phi).exp()
    e1 = Wt.phiI_ball(X) - phi          # log(rho/W) = phi_I - phi
    return Ap, Bp, coefs, rho, e1, X
STRIP = [(29.0, 32.0)]
def f_sB(a, b):
    Ap, Bp, _, _, e1, _ = strip_cell(a, b); return pos(Bp) * e1.exp()
def f_IA(a, b):
    Ap, _, _, rho, _, _ = strip_cell(a, b); return pos(Ap * rho)
Hd = PhiHardy(Wt, 16.25, 35.0); log(f"Phi table {time.time()-T0:.1f}s")
def f_sA(a, b):
    Ap, _, _, rho, _, X = strip_cell(a, b)
    xl, xh = lof(X), upf(X)
    return pos(Ap * rho) * Hd.upper(xl).union(Hd.upper(xh))
def f_dj(j):
    def f(a, b):
        _, _, coefs, _, e1, _ = strip_cell(a, b); return absb(coefs[j]) * (e1 / 2).exp()
    return f
R = {}
R['sB'] = rec('sB', bnb_sup(f_sB, STRIP, n_strip, rel))
R['IA'] = rec('IA', bnb_int(f_IA, STRIP, n_strip, rel))
R['sA_half'] = rec('sA', bnb_int(f_sA, STRIP, n_strip, rel)); R['sA'] = 2 * R['sA_half']; del R['sA_half']
R['K_rho_16.25'] = upf(phi_ball(arb('16.25')).exp())
R['Phi_34.25'] = upf(Hd.upper(34.25))
log(f"strip Q1 done {time.time()-T0:.1f}s", {k: R[k] for k in ('sB', 'IA', 'sA', 'K_rho_16.25', 'Phi_34.25')})
for j in (3, 2, 1, 0): R[f'd{j}'] = rec(f'd{j}', bnb_sup(f_dj(j), STRIP, n_strip, rel))
log(f"strip Q4 done {time.time()-T0:.1f}s", {f'd{j}': R[f'd{j}'] for j in range(4)})

# ------------------------------------------------------------------ inner ramp (x in [10.25,16.25], weight W, cut zeta, 7km)
@lru_cache(maxsize=None)
def ramp_cell(a, b):
    X = Iv(a, b); SIG = X + TIP
    c = cut_zeta(X); phi = phi_ball(X); P = Wt.Pj(a, b); k, m, p, q = kmpq(SIG)
    Ap, Bp = ABf7(*c, *P, k, m, p, q); coefs = CF7(*c, *P, k, m, p, q)
    e2 = phi - Wt.phiI(a, b)           # log(W/rho)
    return Ap, Bp, coefs, e2, c, q
RAMP = [(10.25, 16.25)]
R['a_f'] = rec('a_f', bnb_sup(lambda a, b: pos(ramp_cell(a, b)[0] * ramp_cell(a, b)[3].exp()), RAMP, n_ramp, rel))
R['b_f'] = rec('b_f', bnb_sup(lambda a, b: pos(ramp_cell(a, b)[1] * ramp_cell(a, b)[3].exp()), RAMP, n_ramp, rel))
def f_drift(a, b):
    _, _, _, e2, c, q = ramp_cell(a, b); return (q / 4) * c[0] * c[1] * e2.exp()
R['drift_f'] = rec('drift', bnb_sup(f_drift, RAMP, n_ramp, rel))
log(f"ramp Q2 done {time.time()-T0:.1f}s", {k: R[k] for k in ('a_f', 'b_f', 'drift_f')})
for j in (3, 2, 1, 0):
    R[f'c{j}'] = rec(f'c{j}', bnb_sup(lambda a, b, j=j: absb(ramp_cell(a, b)[2][j]) * (ramp_cell(a, b)[3] / 2).exp(), RAMP, n_ramp, rel))
log(f"ramp Q3 done {time.time()-T0:.1f}s", {f'c{j}': R[f'c{j}'] for j in range(4)})

# ------------------------------------------------------------------ outer ramp (weight W, cut zeta_out, 7km)
OUTER = {}
for Rr in (60, 40):
    @lru_cache(maxsize=None)
    def out_cell(a, b, Rr=Rr):
        X = Iv(a, b); SIG = X + TIP
        c = cut_zo(X, Rr); P = Wt.Pj(a, b); k, m, p, q = kmpq(SIG)
        Ap, Bp = ABf7(*c, *P, k, m, p, q); coefs = CF7(*c, *P, k, m, p, q)
        phiI = Wt.phiI(a, b); W = (-phiI).exp(); Wh = (-phiI / 2).exp()
        return Ap * W, Bp * W, [cf * Wh for cf in coefs]
    rng = [(Rr - 3.0, float(Rr))] if Rr - 3 >= B2 else [(Rr - 3.0, B2), (B2, float(Rr))]
    o = {}
    o['sup|A_zo|'] = rec(f'Azo{Rr}', bnb_sup(lambda a, b: absb(out_cell(a, b)[0]), rng, n_out, rel))
    o['sup|B_zo|'] = rec(f'Bzo{Rr}', bnb_sup(lambda a, b: absb(out_cell(a, b)[1]), rng, n_out, rel))
    for j in (3, 2, 1, 0):
        o[f'sup|coef{j}|W^1/2'] = rec(f'oc{j}_{Rr}', bnb_sup(lambda a, b, j=j: absb(out_cell(a, b)[2][j]), rng, n_out, rel))
    o['W(R-3)'] = upf(Wt.W_point(Rr - 3.0)); o['W(R)'] = upf(Wt.W_point(float(Rr)))
    o['W(R-3)_ball'] = Wt.W_point(Rr - 3.0).str(8); o['W(R)_ball'] = Wt.W_point(float(Rr)).str(8)
    OUTER[Rr] = o
    log(f"outer R={Rr} done {time.time()-T0:.1f}s", o)

# ------------------------------------------------------------------ Q6 mu window
MU1 = {'2': (0.1552, 0.3357), '1.5': (0.1241, 0.2211)}
def mu_window(Rd, mu1_min, mu1_strip, thE=arb(1) / 3, thI=arb(1) / 3, c1=arb('0.5804')):
    om0 = arb(1) / 32 - arb('0.01607')
    af, bf = arb(Rd['a_f']), arb(Rd['b_f'])
    lam_lo, lam_hi = 0.0, float((om0 / (af * c1)).upper())
    def ok(lam):
        l = arb(lam); A = arb(1) / 9 - l * bf; Cc = om0 - l * af * c1; Bq = l * af
        return bool(A > 0) and bool(Cc >= 0) and bool(Bq * Bq <= 4 * A * Cc)
    for _ in range(200):
        mid = 0.5 * (lam_lo + lam_hi)
        if ok(mid): lam_lo = mid
        else: lam_hi = mid
    mu_max = thE * arb(lam_lo)
    mu_min = arb(Rd['sB']) / (thI * arb(str(mu1_strip))) + arb(Rd['sA']) / (thI * arb(str(mu1_min)))
    return dict(mu_min_upper=upf(mu_min), mu_max_lower=lof(mu_max), ratio_lower=lof(mu_max / mu_min), om0=float(om0.mid()),
                mu1_min=mu1_min, mu1_strip=mu1_strip, thE=1 / 3, thI=1 / 3, c1=0.5804)
MU = mu_window(R, *MU1[PXs]); log("mu window:", MU)
res = dict(P=PXs, tip=[float(TIP.lower()), float(TIP.upper())], quantities=R, outer=OUTER, mu_window=MU, n_pieces=NP,
           runtime_s=time.time() - T0, fast=FAST, log=LOG,
           note="all numbers are rigorous UPPER bounds (mu_max is a rigorous LOWER bound); floats rounded outward")
json.dump(res, open(OUT, 'w'), indent=1)
print("done", time.time() - T0)
