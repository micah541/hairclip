"""fix11 review, part D (2): the E0 cubic form.
(a) hypotheses of cubic_E0fine_0.001_0.01_0.01.json (alpha <= 1e-3, |w|,|w'| <= 0.01 on [0,31.25]) against what closure.box_assert asserts and against the run.
(b) the strip sigma in [29,32]: certified sups of the SECOND derivatives of the weighted divergence-form coefficients C_l(U; X, X1, X2) (tc_weighted.pkl, X = chi^2),
    i.e. of the cubic remainder  R3^X = int_0^1 (1-t) sum_{i,j,l} d^2 C_l/dU_i dU_j (u0 + t w) w^(i) w^(j) w^(l) dt  on the strip, over the box
    |alpha| <= 2.5e-4, |w|,|w'|,|w''| <= 0.02 (asserted: U3_strip_S0..S2).  This is the term the verifier covers by the uncertified scaling ADD_STRIP_CUBIC = 2.6e-3.
    Output: monomial coefficients c_(a,b,c) with |R3^X/rho| <= sum c_mono |w|^a |w'|^b |w''|^c, and a size estimate at the production peak [est].
Read-only on the proof tree.  Run: PYTHONHASHSEED=0 nice -n 19 python3 -B d_cubic_strip.py"""
import sys, os, json, hashlib, pickle, time, math
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); os.chdir(HERE)
from flint import arb, ctx
import sympy as sp
ctx.prec = 120
t0 = time.time()
V = '/home/claude/s1spec/step3nl_v4b'
out = {}
# ---------------------------------------------------------------- (a)
CJ = json.load(open(V + '/cubic_E0fine_0.001_0.01_0.01.json'))
rows = {}
for key, f in (('1.5e-4', V + '/runs/P/rate_table_fin13v4b_a0.00015.json'), ('1.25e-4', V + '/runs/P/rate_table_fin13v4b_a0.000125.json')):
    d = json.load(open(f)); st = d['steps']
    rows[key] = dict(far_w_max=max(s['far']['far_w'] for s in st), sI1_max=max(s['far']['sI1'] for s in st), far_w1_max=max(s['far']['far_w1'] for s in st), S_ramp_max=max(s['S_ramp_sup'] for s in st),
                     alpha_max=max(s['alpha_up'] for s in st), sI2_max=max(s['far']['sI2'] for s in st), peak_sqrtL=d['sqrtL_max'])
out['cubic_box'] = dict(certificate=dict(alpha_max=CJ['alpha_max'], eta0=CJ['eta0'], eta1=CJ['eta1'], regions=list(CJ['sup'].keys())),
                        asserted=dict(tip_0_16p25=1e-3, far_w_16p25_34p25=0.02, far_w1=0.02), run=rows,
                        holds_in_run=all(r['far_w_max'] <= 0.01 and r['far_w1_max'] <= 0.01 and r['S_ramp_max'] <= 0.01 and r['alpha_max'] <= 1e-3 for r in rows.values()),
                        asserted_implies_certificate_box=False)
print(json.dumps(out['cubic_box'], indent=1), flush=True)
# ---------------------------------------------------------------- (b)
PK = pickle.load(open(V + '/v4bcert/out/tc_weighted.pkl', 'rb'))
C = PK['C']
U = sp.symbols('U0 U1 U2')
ARGS = "phi1, phi2, k, m, p, X, X1, X2, U0, U1, U2"
allowed = set(sp.symbols(ARGS.replace(' ', '')))
class _Pr(sp.printing.str.StrPrinter):
    def _print_Rational(self, e): return f"(arb({e.p})/arb({e.q}))"
    def _print_Integer(self, e): return f"arb({e.p})"
    def _print_Pow(self, e, rational=False):
        b, ex = e.as_base_exp(); bb = f"({self._print(b)})"
        if ex.is_Integer:
            n = int(ex); prod = "*".join([bb] * abs(n)); return f"({prod})" if n > 0 else f"(1/({prod}))"
        if ex.is_Rational and ex.q == 2:
            n = int(ex.p); base = f"({bb}.sqrt())"; prod = "*".join([base] * abs(n)); return f"({prod})" if n > 0 else f"(1/({prod}))"
        raise ValueError(f"power {ex}")
D2 = {}
for l in range(3):
    for i in range(3):
        for j in range(i, 3):
            e = sp.diff(C[l], U[i], U[j])
            if e == 0: continue
            assert e.free_symbols <= allowed, (l, i, j, e.free_symbols - allowed)
            D2[(l, i, j)] = eval(f"lambda {ARGS}: " + _Pr().doprint(e), {"arb": arb})
print(f"[{time.time()-t0:.0f}s] {len(D2)} nonzero second derivatives", flush=True)
# also: is C_l independent of q (the curvature part) and at most quadratic in U2?
out['structure'] = dict(n_nonzero=len(D2), keys=[list(k_) for k_ in D2])
sys.path.insert(0, '/home/claude/s1spec/step3nl'); sys.path.insert(0, '/home/claude/s1spec/step3')
import h2res as HR
F, Hf = HR.F, HR.Hf
AMAX = arb('2.5e-4'); Ab = arb(0, AMAX); SB = arb(0, arb('0.02'))
Su = sp.symbols('u'); Sm = Su ** 4 * (35 - 84 * Su + 70 * Su ** 2 - 20 * Su ** 3)
dSm = [eval("lambda u: " + _Pr().doprint(sp.expand(sp.diff(Sm, Su, j))), {"arb": arb}) for j in range(3)]
def chi_jets(sa, sb):
    tl = (32 - sb) / 3; th = (32 - sa) / 3; T = arb((tl + th) / 2, 0).union(tl).union(th)
    return dSm[0](T), -dSm[1](T) / 3, dSm[2](T) / 9
def phis(Xv):
    E = (1 + Xv * Xv) ** (arb(-1) / 3); return Xv * E / 3, E / 3 - 2 * Xv * Xv * E * E * E * E / 9
mono = {}; NS = 200; IR = arb(0); IR2 = arb(0); RWs = 0.0
for ci in range(29, 32):
    for i in range(NS):
        sa = arb(ci) + arb(i) / NS; sb = arb(ci) + arb(i + 1) / NS
        kb, mb, pb, qs = HR.HS.profile_bounds(sa)
        xa = float((sa - HR.ST).lower()); xb = float((sb - HR.ST).upper()); xc = (xa + xb) / 2
        Xv = arb(xc, (xb - xa) / 2); ph1, ph2 = phis(Xv)
        y = HR.jets_any(F, Xv, xc)[:3]; h = HR.jets_any(Hf, Xv, xc)[:3]
        kk, mm, pp = arb(0, kb.upper()), arb(0, mb.upper()), arb(0, pb.upper())
        ch, ch1, ch2 = chi_jets(sa, sb); Xw = ch * ch; X1w = 2 * ch * ch1; X2w = 2 * ch1 * ch1 + 2 * ch * ch2
        Ub = [Ab * y[j] + Ab * Ab * h[j] + SB for j in range(3)]            # u0^(j) + t w^(j), |w^(j)| <= 0.02, t in [0,1]
        for (l, i_, j_), f in D2.items():
            v = abs(f(ph1, ph2, kk, mm, pp, Xw, X1w, X2w, *Ub)); vu = float(v.upper()); assert vu == vu and math.isfinite(vu)
            mult = 1.0 if i_ == j_ else 2.0                                   # sum over ordered (i,j)
            cnt = [0, 0, 0]
            for o in (l, i_, j_): cnt[o] += 1
            key = (tuple(cnt), (l, i_, j_))
            mono[key] = max(mono.get(key, 0.0), 0.5 * mult * vu)             # (1/2) from int_0^1 (1-t) dt
        rho_hi = float(arb(HR.rho_of(xa)).upper()); L = arb(xb) - arb(xa)
        IR += L * arb(rho_hi); IR2 += L * arb(rho_hi) ** 2; RWs = max(RWs, rho_hi)
tot = {}
for (cnt, lij), v in mono.items(): tot[cnt] = tot.get(cnt, 0.0) + v        # sum of sups >= sup of the sum
C4 = json.load(open(V + '/v4cert/out/constants_v4.json'))['constants']
KW = arb(repr(C4['KW']['value'])); RW = arb(RWs) * KW; IRW = (IR2 * KW).sqrt()
out['strip_cubic_coefficients'] = {str(k_): float(arb(repr(v)).upper()) * (1 + 1e-15) for k_, v in sorted(tot.items())}
out['strip_integrals'] = dict(IR=float(IR.upper()), IRW=float(IRW.upper()), RW=float(RW.upper()))
# size estimate at the production peak (alpha1 = 1.5e-4) [est: worst case per factor, not a certified term of the inequality]
r = rows['1.5e-4']; s = r['peak_sqrtL']; MU, EPSI = 0.08, 3e-5; C1W = C4['C1W']['value']
V0 = s / math.sqrt(MU); V2 = s / math.sqrt(MU * EPSI); V1 = math.sqrt(V0 * (V2 + C1W * V0))
Ssup = [r['far_w_max'], r['sI1_max'], r['sI2_max']]
def pair(o1, o2):          # int_strip rho |w^(o1)||w^(o2)|
    Vn = {1: V0, 2: V1}
    if o1 == 0 and o2 == 0: return float(IR.upper()) * Ssup[0] ** 2
    if o2 == 0: return float(IRW.upper()) * Ssup[0] * Vn[o1]
    return float(RW.upper()) * Vn[o1] * Vn[o2]
est = 0.0; parts = {}
for cnt, c in tot.items():
    orders = sorted([0] * cnt[0] + [1] * cnt[1] + [2] * cnt[2], reverse=True)
    # sup on the LOWEST-order factor (best available sup), pair = the two highest
    val = c * pair(orders[0], orders[1]) * Ssup[orders[2]]
    parts[str(cnt)] = val; est += val
Lam = s * s
MU1 = C4['MU1']['value']
out['size_estimate_at_peak'] = dict(strip_cubic_upper_est=est, over_Lambda=est / Lam, parts=parts, V0=V0, V1=V1, V2=V2, sups=Ssup,
                                    note='[est] every norm at its own maximum for Lyap = peak^2 (mutually exclusive maxima): an over-estimate of the rate contribution; the E_I dissipation at the V2-saturated shape is mu(1-theta)/2 V2^2/Lambda = %.3g per unit Lambda' % (0.5 * 0.95 / EPSI))
out['args'] = dict(argv=sys.argv, env={k_: os.environ.get(k_) for k_ in ('PYTHONHASHSEED', 'PYTHONDONTWRITEBYTECODE')},
                   inputs={p_: hashlib.sha256(open(p_, 'rb').read()).hexdigest() for p_ in (V + '/v4bcert/out/tc_weighted.pkl', V + '/cubic_E0fine_0.001_0.01_0.01.json', V + '/runs/P/rate_table_fin13v4b_a0.00015.json',
                                                                                           V + '/runs/P/rate_table_fin13v4b_a0.000125.json', V + '/v4cert/out/constants_v4.json')},
                   script_sha256=hashlib.sha256(open(os.path.abspath(__file__), 'rb').read()).hexdigest(), runtime_s=round(time.time() - t0, 1))
json.dump(out, open(os.path.join(HERE, 'd_cubic_strip.json'), 'w'), indent=1, sort_keys=True)
print(json.dumps({k_: out[k_] for k_ in ('strip_cubic_coefficients', 'strip_integrals', 'size_estimate_at_peak')}, indent=1))
