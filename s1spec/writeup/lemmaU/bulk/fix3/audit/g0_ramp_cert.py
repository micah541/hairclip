"""G0 ramp audit: small certified (arb, outward) sups used in AUDIT_RAMP.md.  Read-only w.r.t. the proof tree.
usage: nice -n 19 python3 g0_ramp_cert.py  ->  g0_ramp_cert.json (same directory)
All sups are rigorous UPPER bounds (max over a cover of the interval by closed ball pieces, each piece evaluated in arb);
infs are rigorous LOWER bounds.  Definitions follow the verifier sources:
  zeta inner = 1 - S7((16.25-x)/6) on [10.25,16.25]      (agent_coup/coup_common.py cut_zeta)
  zeta outer = S7((60-x)/3) on [57,60]                    (cut_zo, R=60)
  S7(u) = u^4 (35 - 84u + 70u^2 - 20u^3)
  psi' = (1-s) 4x/(1+x^2) + s (x/3)(1+x^2)^(-1/3),  s = S7((x-30.25)/8) on [30.25,38.25]   (coup_common.Weight, P=2)
  rho = exp(-(1+x^2)^(2/3)/4);  W = (1+x^2)^-2 on x <= 30.25
  q = sigma + [-qs,qs], sigma = x + TIP, |q'-1| <= kb*pb  (step3/hstail.profile_bounds at sigma >= 8, monotone)"""
import sys, os, json, hashlib, time
from flint import arb, arb_poly, ctx
ctx.prec = 200
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, '/home/claude/s1spec/step3')
import hstail as H
T0 = time.time()

TIP = arb('-2.2490776515195').union(arb('-2.2490776505754'))
kb, mb, pb, qs = H.profile_bounds(arb(8))
KB, PB, QS = arb(kb.upper()), arb(pb.upper()), arb(qs.upper())
QP_DEV = KB * PB                      # |q' - 1| = |k p|
S7 = arb_poly([0, 0, 0, 0, 35, -84, 70, -20])
S7d = [S7]
for j in range(5): S7d.append(S7d[-1].derivative())

def hull(lo, hi): return arb(lo).union(arb(hi))
def up(z): return float(arb(z.upper()).mid()) if False else z.upper()
def zeta_in(X, j):
    u = (arb('16.25') - X) / 6
    v = -S7d[j](u) * (arb(-1) / 6) ** j
    return (1 + v) if j == 0 else v
def zeta_out(X, j):
    u = (arb(60) - X) / 3
    return S7d[j](u) * (arb(-1) / 3) ** j
def psip(X, a, b):
    Ap = 4 * X / (1 + X * X)
    Bp = X * (1 + X * X) ** (-arb(1) / 3) / 3
    if b <= 30.25: return Ap
    if a >= 38.25: return Bp
    s = S7((X - arb('30.25')) / 8)
    return (1 - s) * Ap + s * Bp
def W_over_rho(X):     # x <= 30.25 only
    return ((1 + X * X) ** (arb(2) / 3) / 4 - 2 * (1 + X * X).log()).exp()
def qball(X): return X + TIP + arb(0, 1) * QS

def cover(a, b, n):
    h = (b - a) / n
    for i in range(n):
        lo, hi = a + i * h, (a + (i + 1) * h if i < n - 1 else b)
        yield lo, hi, hull(lo, hi)
def sup(f, a, b, n):
    m = None
    for lo, hi, X in cover(a, b, n):
        v = arb(f(X, lo, hi).upper())
        m = v if m is None or v > m else m
    return m
def inf(f, a, b, n):
    m = None
    for lo, hi, X in cover(a, b, n):
        v = arb(f(X, lo, hi).lower())
        m = v if m is None or v < m else m
    return m
def regs(a, b):
    pts = [a] + [c for c in (30.25, 38.25) if a < c < b] + [b]
    return list(zip(pts[:-1], pts[1:]))
def sup_regs(f, a, b, n_per_unit=400):
    return max((sup(f, u, v, max(8, int((v - u) * n_per_unit))) for u, v in regs(a, b)), key=lambda z: float(z.mid()))
def inf_regs(f, a, b, n_per_unit=400):
    return min((inf(f, u, v, max(8, int((v - u) * n_per_unit))) for u, v in regs(a, b)), key=lambda z: float(z.mid()))
F = lambda z: float(arb(z.upper()).mid()) if hasattr(z, 'upper') else float(z)
def fu(z):  # float upper (rounded up)
    import math
    x = float(z.upper().mid()); return math.nextafter(x, math.inf) if x != 0 else 0.0
def fl(z):
    import math
    x = float(z.lower().mid()); return math.nextafter(x, -math.inf) if x != 0 else 0.0

R = {}
# ---------------------------------------------------------------- (1) rounding reserves of the typed constants
raw = json.load(open('/home/claude/s1spec/step3nl_v33/agent_coup/coupling_constants.json'))['P']['2']['rigorous']
typed = dict(c3='0.4421', c2='0.5259', c1='0.4235', c0='0.4250', a_f='6.9555e-2', b_f='3.9359e-2',
             sB='3.5955e-5', sA='1.3372e-4', IA='9.5997e-12', d3='6.841e-3', d2='1.5418e-2', d1='2.3648e-2', d0='4.5627e-2',
             Phi_34_25='9.8547e6')
rawkey = dict(Phi_34_25='Phi_34.25')
res = {}
for k, t in typed.items():
    r = arb(repr(raw[rawkey.get(k, k)])); tt = arb(t)
    rel = (tt - r) / r
    res[k] = dict(raw_certified_upper=raw[rawkey.get(k, k)], typed_in_verifier=t, typed_minus_raw=fu(tt - r) if (tt - r) > 0 else -fu(r - tt),
                  relative=float(rel.mid()), covered=bool(tt >= r))
R['rounding_reserve'] = res

# ---------------------------------------------------------------- (2) level-0 IBP factor sup(|zeta'| + psi' zeta) on supp zeta (verifier 1.67)
k_in = sup(lambda X, a, b: abs(zeta_in(X, 1)) + psip(X, a, b) * zeta_in(X, 0), 10.25, 16.25, 6000)
k_mid = sup_regs(lambda X, a, b: psip(X, a, b), 16.25, 57.0)
k_out = sup(lambda X, a, b: abs(zeta_out(X, 1)) + psip(X, a, b) * zeta_out(X, 0), 57.0, 60.0, 6000)
R['kappa_IBP'] = dict(inner=fu(k_in), middle_sup_psi1=fu(k_mid), outer=fu(k_out), verifier_value=1.67,
                      covered=bool(max(k_in, k_mid, k_out, key=lambda z: float(z.mid())) <= arb('1.67')))
# ---------------------------------------------------------------- (3) transport (a) and (b) coefficients, levels 0 and 2, on [10.25,60]
qp_hi, qp_lo = 1 + QP_DEV, 1 - QP_DEV
qpb = arb(1) + arb(0, 1) * QP_DEV
lev0 = lambda X, a, b: psip(X, a, b) * qball(X) / 2 + qpb / 2 - 1            # coefficient of int W v^2 (a-transport, level 0)
lev2 = lambda X, a, b: psip(X, a, b) * qball(X) / 2 + 5 * qpb / 2 - 1        # coefficient of int W v''^2 (a-transport, level 2)
s0, i0 = sup_regs(lev0, 10.25, 60.0), inf_regs(lev0, 10.25, 60.0)
s2, i2 = sup_regs(lev2, 10.25, 60.0), inf_regs(lev2, 10.25, 60.0)
psi1 = sup_regs(lambda X, a, b: psip(X, a, b), 10.25, 60.0)
psimin = inf_regs(lambda X, a, b: psip(X, a, b), 10.25, 60.0)
cb0 = psi1 / 2 + 2 * KB          # b-transport level 0: (1/2) sup sqrt(psi'^2 + k^2) <= psi1/2 + k
cb2 = psi1 / 2 + 3 * KB          # level 2: (1/2)(psi' cos + k sin) - 2 k sin  ->  <= psi1/2 + 5k/2
ca0 = s0 if s0 >= -i0 else -i0
comb0 = (ca0 * ca0 + cb0 * cb0).sqrt(); comb2 = (s2 * s2 + cb2 * cb2).sqrt()
R['transport'] = dict(q_minus_sigma_ball=fu(QS), qprime_dev=fu(QP_DEV),
                      level0_a_coef_sup=fu(s0), level0_a_coef_inf=fl(i0), level0_b_coef=fu(cb0), level0_euclid=fu(comb0),
                      level0_verifier=39.16, level0_covered=bool(comb0 <= arb('39.16')),
                      level2_a_coef_sup=fu(s2), level2_b_coef=fu(cb2), level2_euclid=fu(comb2), level2_verifier_V2sq=39.16,
                      level2_deficit_V2sq=fu(comb2 - arb('39.16')) if comb2 > arb('39.16') else 0.0,
                      sup_psi_prime=fu(psi1), inf_psi_prime=fl(psimin),
                      hardy_V2_le_c_V3=fu(2 / psimin),
                      level2_deficit_as_V2V3_coef=fu((comb2 - arb('39.16')) * 2 / psimin) if comb2 > arb('39.16') else 0.0,
                      verifier_V2V3_coef_per_ab=fu(arb(62) / arb('7.57')))
# ---------------------------------------------------------------- (4) inner-ramp weights and commutator coefficient sups
wr = sup(lambda X, a, b: W_over_rho(X), 10.25, 16.25, 6000)
qz_wr = sup(lambda X, a, b: qball(X) * zeta_in(X, 1) * W_over_rho(X), 10.25, 16.25, 6000)
qz_swr = sup(lambda X, a, b: qball(X) * zeta_in(X, 1) * W_over_rho(X).sqrt(), 10.25, 16.25, 6000)
z_swr = sup(lambda X, a, b: zeta_in(X, 1) * W_over_rho(X).sqrt(), 10.25, 16.25, 6000)
binom = {1: 4, 2: 6, 3: 4, 4: 1}
ek = {k: sup(lambda X, a, b, k=k: binom[k] * abs(zeta_in(X, k)) * W_over_rho(X).sqrt(), 10.25, 16.25, 6000) for k in (1, 2, 3, 4)}
# planB ramp table, alpha <= 2.5e-4 row, S_ramp <= 7.48e-4 (fin12f4 certified S_ramp 7.473592e-04)
tab = json.load(open('/home/claude/s1spec/step3nl_v33/agent_int/planB/planB_ramp.json'))['table']
c0 = arb(repr(tab['alpha=0.00025,eta=0']['d001_w4'])); c1 = arb(repr(tab['alpha=0.00025,eta=0.01']['d001_w4']))
d001r = c0 + (c1 - c0) / arb('0.01') * arb('7.48e-4')
CF = [arb(x) for x in ('0.4421', '0.5259', '0.4235', '0.4250')]   # multiply D4 (g'''), W3 (g''), B (g'), W1 (g)
# quasilinear commutator d001*[d^4, zeta]g : coefficient of g^(4-k) is C(4,k) zeta^(k); g''' <-> D4 (k=1), g'' <-> W3 (k=2), g' <-> B (k=3), g <-> W1 (k=4)
ql = [d001r * ek[k] for k in (1, 2, 3, 4)]
R['inner_ramp'] = dict(sup_W_over_rho=fu(wr), sqrt=fu(wr.sqrt()), term_a_factor_0_66_covered=bool(wr.sqrt() <= arb('0.66')),
                       sup_q_zeta1_W_over_rho=fu(qz_wr), verifier_CAB25_per_ab=2.5,
                       sup_q_zeta1_sqrtWrho=fu(qz_swr), sup_zeta1_sqrtWrho=fu(z_swr),
                       level2_mod_commutator_W1_coef_per_ab=fu((qz_swr ** 2 + z_swr ** 2).sqrt()),
                       binom_zeta_k_sqrtWrho={k: fu(v) for k, v in ek.items()},
                       d001_ramp_bound=fu(d001r),
                       quasilinear_commutator_coefs_D4_W3_B_W1=[fu(v) for v in ql],
                       quasilinear_rel_to_CF=[fu(v / c) for v, c in zip(ql, CF)])
# reserve needed vs available (relative to CF_j), for the two inner-ramp uncovered level-2 terms at the certified sup|(a,b)|
abar = arb('2.19103e-7')
modc = (qz_swr ** 2 + z_swr ** 2).sqrt() * abar
R['kill_check'] = dict(abar=2.19103e-7, mod_commutator_W1_coef=fu(modc), mod_rel_CF3=fu(modc / CF[3]),
                       CF_reserve_rel={k: res[k]['relative'] for k in ('c3', 'c2', 'c1', 'c0')},
                       AF_reserve_rel=res['a_f']['relative'], BF_reserve_rel=res['b_f']['relative'])
# ---------------------------------------------------------------- (5) D2v = (W v'')''/W = v'''' - 2 psi' v''' + (psi'^2 - psi'') v''  (verifier: V4 + 2.61 V3 + 1.7 V2)
def psipp(X, a, b):
    Ap2 = 4 * (1 - X * X) / (1 + X * X) ** 2
    Bp2 = (1 + X * X) ** (-arb(1) / 3) / 3 - 2 * X * X * (1 + X * X) ** (-arb(4) / 3) / 9
    if b <= 30.25: return Ap2
    if a >= 38.25: return Bp2
    t = (X - arb('30.25')) / 8; s = S7(t); s1 = S7d[1](t) / 8
    Ap = 4 * X / (1 + X * X); Bp = X * (1 + X * X) ** (-arb(1) / 3) / 3
    return (1 - s) * Ap2 + s * Bp2 + s1 * (Bp - Ap)
c2v = sup_regs(lambda X, a, b: abs(psip(X, a, b) ** 2 - psipp(X, a, b)), 10.25, 60.0)
R['D2v'] = dict(two_psi1=fu(2 * psi1), verifier_2_61=bool(2 * psi1 <= arb('2.61')), sup_abs_psi1sq_minus_psi2=fu(c2v), verifier_1_7=bool(c2v <= arb('1.7')))
# ---------------------------------------------------------------- args block
ins = ['/home/claude/s1spec/step3nl_v33/agent_coup/coupling_constants.json', '/home/claude/s1spec/step3nl_v33/agent_int/planB/planB_ramp.json',
       '/home/claude/s1spec/step3/hstail.py', os.path.abspath(__file__)]
R['args'] = dict(argv=sys.argv, env={k: os.environ.get(k) for k in ('OMP_NUM_THREADS',)}, prec=ctx.prec,
                 inputs={p: hashlib.sha256(open(p, 'rb').read()).hexdigest() for p in ins},
                 profile_bounds_sigma8=dict(kb=fu(KB), pb=fu(PB), qs=fu(QS)), runtime_s=round(time.time() - T0, 1))
json.dump(R, open(os.path.join(HERE, 'g0_ramp_cert.json'), 'w'), indent=1, sort_keys=True)
print(json.dumps({k: v for k, v in R.items() if k != 'args'}, indent=1, sort_keys=True))
