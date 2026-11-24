#!/usr/bin/env python3
"""Stage C1 (fin13v4): certified cutoff / weight sups for the E_I inequality (EI_SPEC 1.2, 3, 4.2-4.4, 4.9).
Copied from writeup/lemmaU/bulk/fix3/audit/g0_ramp_cert.py (definitions identical: zeta inner 1 - S7((16.25-x)/6), outer S7((60-x)/3),
psi' = (1-s) 4x/(1+x^2) + s (x/3)(1+x^2)^(-1/3), s = S7((x-30.25)/8), W = (1+x^2)^-2 on x <= 30.25, rho = exp(-(1+x^2)^(2/3)/4)),
extended by: k3_j, Z2i, KWIR, sup|zeta'| on IR and OR, T0/T2 with the exact b-forms (psi'/2 + |k|/2, psi'/2 + 5|k|/2), TC0 and TC2 as Euclidean
per-piece sups, Q2B (profile only), and the (1 - |k| U)^-1 factor for Q_u, C_u (|k U| from stage C2 is used there; here only the profile part).
All sups: max over a cover of closed pieces, each piece evaluated in arb (200 bits); upper endpoints; float conversion nextafter(+inf).
Output out/c1_cut_sups.json.  Run: PYTHONHASHSEED=0 nice -n 19 python3 c1_cut_sups.py"""
import sys, os, json, time, math
from flint import arb, arb_poly, ctx
ctx.prec = 200
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from v4common import args_block, dump, OUT, fup, fdn
sys.path.insert(0, '/home/claude/s1spec/step3')
import hstail as H
T0_ = time.time()
TIP = arb('-2.2490776515195').union(arb('-2.2490776505754'))
kb, mb, pb, qs = H.profile_bounds(arb(8))
KB, MB, PB, QS = arb(kb.upper()), arb(mb.upper()), arb(pb.upper()), arb(qs.upper())
S7 = arb_poly([0, 0, 0, 0, 35, -84, 70, -20])
S7d = [S7]
for j in range(5): S7d.append(S7d[-1].derivative())
def hull(lo, hi): return arb(lo).union(arb(hi))
def zeta_in(X, j):
    u = (arb('16.25') - X) / 6
    v = -S7d[j](u) * (arb(-1) / 6) ** j
    return (1 + v) if j == 0 else v
def zeta_out(X, j):
    u = (arb(60) - X) / 3
    return S7d[j](u) * (arb(-1) / 3) ** j
def psip(X, a, b):
    Ap = 4 * X / (1 + X * X); Bp = X * (1 + X * X) ** (-arb(1) / 3) / 3
    if b <= 30.25: return Ap
    if a >= 38.25: return Bp
    s = S7((X - arb('30.25')) / 8)
    return (1 - s) * Ap + s * Bp
def psipp(X, a, b):
    Ap2 = 4 * (1 - X * X) / (1 + X * X) ** 2
    Bp2 = (1 + X * X) ** (-arb(1) / 3) / 3 - 2 * X * X * (1 + X * X) ** (-arb(4) / 3) / 9
    if b <= 30.25: return Ap2
    if a >= 38.25: return Bp2
    t = (X - arb('30.25')) / 8; s = S7(t); s1 = S7d[1](t) / 8
    Ap = 4 * X / (1 + X * X); Bp = X * (1 + X * X) ** (-arb(1) / 3) / 3
    return (1 - s) * Ap2 + s * Bp2 + s1 * (Bp - Ap)
def W_over_rho(X): return ((1 + X * X) ** (arb(2) / 3) / 4 - 2 * (1 + X * X).log()).exp()     # x <= 30.25 only
def qball(X): return X + TIP + arb(0, 1) * QS
def cover(a, b, n):
    h = (b - a) / n
    for i in range(n):
        lo, hi = a + i * h, (a + (i + 1) * h if i < n - 1 else b)
        yield lo, hi, hull(lo, hi)
def sup(f, a, b, n):
    m = None
    for lo, hi, X in cover(a, b, n):
        v = f(X, lo, hi); assert v.is_finite()
        v = arb(v.upper())
        m = v if m is None or v > m else m
    return m
def inf(f, a, b, n):
    m = None
    for lo, hi, X in cover(a, b, n):
        v = f(X, lo, hi); assert v.is_finite()
        v = arb(v.lower())
        m = v if m is None or v < m else m
    return m
def regs(a, b):
    pts = [a] + [c for c in (30.25, 38.25) if a < c < b] + [b]
    return list(zip(pts[:-1], pts[1:]))
def amax(lst): return max(lst, key=lambda z: float(z.mid()))
def amin(lst): return min(lst, key=lambda z: float(z.mid()))
def sup_regs(f, a, b, n_per_unit=400): return amax([sup(f, u, v, max(8, int((v - u) * n_per_unit))) for u, v in regs(a, b)])
def inf_regs(f, a, b, n_per_unit=400): return amin([inf(f, u, v, max(8, int((v - u) * n_per_unit))) for u, v in regs(a, b)])
def zeta_any(X, a, b, j):
    """zeta^(j) on a piece lying in one of IR, Z1, OR (or the flat parts)"""
    if b <= 16.25: return zeta_in(X, j)
    if a >= 57.0: return zeta_out(X, j)
    return arb(1) if j == 0 else arb(0)
R = {}
NIR = 6000
# -------- weights: PSI1 = sup psi', PSI2 = sup|psi'^2 - psi''| on [10.25,60]; RWR = sup_IR (W/rho)^1/2; KWIR = 1/W(16.25)
psi1 = sup_regs(lambda X, a, b: psip(X, a, b), 10.25, 60.0)
psimin = inf_regs(lambda X, a, b: psip(X, a, b), 10.25, 60.0)
psi2 = sup_regs(lambda X, a, b: abs(psip(X, a, b) ** 2 - psipp(X, a, b)), 10.25, 60.0)
wr = sup(lambda X, a, b: W_over_rho(X), 10.25, 16.25, NIR)
KWIR = (1 + arb('16.25') ** 2) ** 2          # 1/W(16.25), exact (W = (1+x^2)^-2 on x <= 30.25, decreasing)
R['weight'] = dict(PSI1=fup(psi1), two_PSI1=fup(2 * psi1), PSI1_inf=fdn(psimin), PSI2=fup(psi2), RWR=fup(wr.sqrt()), sup_W_over_rho_IR=fup(wr), KWIR=fup(KWIR))
# -------- KAP_zone = sup(|zeta'| + psi' zeta)
k_in = sup(lambda X, a, b: abs(zeta_in(X, 1)) + psip(X, a, b) * zeta_in(X, 0), 10.25, 16.25, NIR)
k_mid = sup_regs(lambda X, a, b: psip(X, a, b), 16.25, 57.0)
k_out = sup(lambda X, a, b: abs(zeta_out(X, 1)) + psip(X, a, b) * zeta_out(X, 0), 57.0, 60.0, NIR)
R['kappa'] = dict(KAP_IR=fup(k_in), KAP_Z1=fup(k_mid), KAP_OR=fup(k_out), KAP=fup(amax([k_in, k_mid, k_out])))
# -------- cutoff derivative sups
zi = {j: sup(lambda X, a, b, j=j: abs(zeta_in(X, j)), 10.25, 16.25, NIR) for j in (1, 2, 3, 4)}
zo = {j: sup(lambda X, a, b, j=j: abs(zeta_out(X, j)), 57.0, 60.0, NIR) for j in (1, 2, 3, 4)}
R['zeta'] = dict(IR={str(j): fup(v) for j, v in zi.items()}, OR={str(j): fup(v) for j, v in zo.items()},
                 Z1i=fup(zi[1]), Z2i=fup(zi[2]), Z1o=fup(zo[1]),
                 piece_sup_minus_35_96=fup(zi[1] - arb(35) / 96), piece_sup_minus_35_48=fup(zo[1] - arb(35) / 48))
# exact: S7'(t) = 140 t^3 (1-t)^3 (polynomial identity, checked coefficientwise), max on [0,1] at t = 1/2: 140/64;  zeta' = S7'/6 (IR), S7'/3 (OR)
import sympy as _sp
_t = _sp.symbols('t'); _S = _t**4 * (35 - 84*_t + 70*_t**2 - 20*_t**3)
assert _sp.expand(_sp.diff(_S, _t) - 140 * _t**3 * (1 - _t)**3) == 0
R['zeta']['exact_sup_zeta1_IR'] = '35/96 (= 140/64/6)'; R['zeta']['exact_sup_zeta1_OR'] = '35/48 (= 140/64/3)'
# k4_j = sup_IR |binom * zeta^(j)| (W/rho)^1/2 ; k3_j likewise with the K3 binomials (3, 3, 1)
b4 = {1: 4, 2: 6, 3: 4, 4: 1}; b3 = {1: 3, 2: 3, 3: 1}
k4 = [sup(lambda X, a, b, j=j: b4[j] * abs(zeta_in(X, j)) * W_over_rho(X).sqrt(), 10.25, 16.25, NIR) for j in (1, 2, 3, 4)]
k3 = [sup(lambda X, a, b, j=j: b3[j] * abs(zeta_in(X, j)) * W_over_rho(X).sqrt(), 10.25, 16.25, NIR) for j in (1, 2, 3)]
R['k4'] = [fup(v) for v in k4]; R['k3'] = [fup(v) for v in k3]
# -------- transport T0, T2 (EI_SPEC 4.9 T8), exact b-forms; Euclidean norm per piece
qpb = arb(1) + arb(0, 1) * KB * PB              # q' = 1 + k p
def T_euclid(lev):
    ca = (lambda X, a, b: psip(X, a, b) * qball(X) / 2 + qpb / 2 - 1) if lev == 0 else (lambda X, a, b: psip(X, a, b) * qball(X) / 2 + 5 * qpb / 2 - 1)
    cb = (lambda X, a, b: psip(X, a, b) / 2 + KB / 2) if lev == 0 else (lambda X, a, b: psip(X, a, b) / 2 + 5 * KB / 2)
    return sup_regs(lambda X, a, b: (ca(X, a, b) ** 2 + cb(X, a, b) ** 2).sqrt(), 10.25, 60.0)
T0 = T_euclid(0); T2 = T_euclid(2)
# T9 on IR: TC0 = sup |(q, 1)| zeta |zeta'| W/rho ; TC2 = sup |(q, 1)| |zeta'| (W/rho)^1/2
TC0 = sup(lambda X, a, b: (qball(X) ** 2 + 1).sqrt() * zeta_in(X, 0) * abs(zeta_in(X, 1)) * W_over_rho(X), 10.25, 16.25, NIR)
TC2 = sup(lambda X, a, b: (qball(X) ** 2 + 1).sqrt() * abs(zeta_in(X, 1)) * W_over_rho(X).sqrt(), 10.25, 16.25, NIR)
# Q2B = sup |(q'', (cos theta)'')|, q'' = m p - k^2 q, (cos)'' = -m sin - k^2 cos ; |sin|,|cos| <= 1; sup over [10.25,60] (q <= 60 - TIP + QS)
qmax = arb(60) - TIP.lower() + QS
Q2B = ((MB * PB + KB * KB * qmax) ** 2 + (MB + KB * KB) ** 2).sqrt()
R['transport'] = dict(T0=fup(T0), T2=fup(T2), TC0=fup(TC0), TC2=fup(TC2), Q2B=fup(Q2B))
R['profile_sigma8'] = dict(kb=fup(KB), mb=fup(MB), pb=fup(PB), qs=fup(QS))
# -------- Phi(x) = int_16.25^x 1/W on Z1 is computed in C4 (hardy); here only Phi(34.25) cross-check of the typed PHI34 (W = (1+x^2)^-2 on [16.25,30.25], blended beyond)
ins = ['/home/claude/s1spec/step3/hstail.py']
R['args'] = args_block(ins, dict(prec=ctx.prec, pieces_IR=NIR))
R['runtime_s'] = round(time.time() - T0_, 1)
h = dump(os.path.join(OUT, 'c1_cut_sups.json'), R)
print(json.dumps({k: v for k, v in R.items() if k not in ('args',)}, indent=1, sort_keys=True)); print('sha256', h)
