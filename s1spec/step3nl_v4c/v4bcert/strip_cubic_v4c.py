"""fin13v4c certification stage (REVIEW_v4b_v36_proofs.md D-3, section 4.3): certified coefficients of the STRIP part of the E0 cubic form,
replacing the uncertified scaling ADD_STRIP_CUBIC = 2.6e-3 of fin13v4b.  Promoted from fix11_review/d_cubic_strip.py part (b), same mathematics:
  R3^X = int_0^1 (1-t) sum_{l,i,j} d^2 C_l/dU_i dU_j (u0 + t w) w^(i) w^(j) w^(l) dt     (weighted divergence form of tc_weighted.pkl, X = chi^2)
on the strip sigma in [29,32] (x in [31.25,34.25]), over the box |alpha| <= 2.5e-4, |w|,|w'|,|w''| <= 0.02 (all asserted in closure.box_assert:
E0_alpha, U3_strip_S0..S2).  Output: monomial coefficients c_(a,b,c) with |R3^X / rho| <= sum c_(a,b,c) |w|^a |w'|^b |w''|^c on the strip.
Changes against the review script: every sup and every sum is kept in arb (no float sums); the verifier value is
math.nextafter(float(upper), +inf) of the arb total (raw arb string stored); inputs pinned by sha256; part (a) of the review script
(run-box comparison) is not part of the chain (it is the closure D-2a assert now).
Labels: [cert] arb with outward rounding on certified inputs taken as given (profile hulls, jets of psi3^/h2 from step3nl).
Run (chain): PYTHONHASHSEED=0 nice -n 19 python3 -B strip_cubic_v4c.py  -> out/strip_cubic_v4c.json"""
import sys, os, json, pickle, time, math
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import v4bcommon as VC
from flint import arb, ctx
import sympy as sp
ctx.prec = 120
t0 = time.time()
PKL = os.path.join(VC.OUT, 'tc_weighted.pkl')
E0C = os.path.join(VC.OUT, 'e0v4b_constants.json')
PK = pickle.load(open(PKL, 'rb'))
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
STEP3NL = '/home/claude/s1spec/step3nl'
sys.path.insert(0, STEP3NL); sys.path.insert(0, '/home/claude/s1spec/step3')
import h2res as HR
F, Hf = HR.F, HR.Hf
AMAX = arb('2.5e-4'); Ab = arb(0, AMAX); SBOX = arb(0, arb('0.02'))
Su = sp.symbols('u'); Sm = Su ** 4 * (35 - 84 * Su + 70 * Su ** 2 - 20 * Su ** 3)
dSm = [eval("lambda u: " + _Pr().doprint(sp.expand(sp.diff(Sm, Su, j))), {"arb": arb}) for j in range(3)]
def chi_jets(sa, sb):
    tl = (32 - sb) / 3; th = (32 - sa) / 3; T = arb((tl + th) / 2, 0).union(tl).union(th)
    return dSm[0](T), -dSm[1](T) / 3, dSm[2](T) / 9
def phis(Xv):
    E = (1 + Xv * Xv) ** (arb(-1) / 3); return Xv * E / 3, E / 3 - 2 * Xv * Xv * E * E * E * E / 9
def amax(a, b): return b if (a is None or b.upper() > a.upper()) else a
mono = {}; NS = 200; IR = arb(0); IR2 = arb(0); RWs = None; ncell = 0
for ci in range(29, 32):
    for i in range(NS):
        sa = arb(ci) + arb(i) / NS; sb = arb(ci) + arb(i + 1) / NS
        kb, mb, pb, qs = HR.HS.profile_bounds(sa)
        xa = float((sa - HR.ST).lower()); xb = float((sb - HR.ST).upper()); xc = (xa + xb) / 2
        Xv = arb(xc, (xb - xa) / 2); ph1, ph2 = phis(Xv)
        y = HR.jets_any(F, Xv, xc)[:3]; h = HR.jets_any(Hf, Xv, xc)[:3]
        kk, mm, pp = arb(0, kb.upper()), arb(0, mb.upper()), arb(0, pb.upper())
        ch, ch1, ch2 = chi_jets(sa, sb); Xw = ch * ch; X1w = 2 * ch * ch1; X2w = 2 * ch1 * ch1 + 2 * ch * ch2
        Ub = [Ab * y[j] + Ab * Ab * h[j] + SBOX for j in range(3)]          # u0^(j) + t w^(j), |w^(j)| <= 0.02, t in [0,1]
        for (l, i_, j_), f in D2.items():
            v = abs(f(ph1, ph2, kk, mm, pp, Xw, X1w, X2w, *Ub)); assert v.is_finite()
            mult = 1 if i_ == j_ else 2                                          # sum over ordered (i,j)
            cnt = [0, 0, 0]
            for o in (l, i_, j_): cnt[o] += 1
            key = (tuple(cnt), (l, i_, j_))
            mono[key] = amax(mono.get(key), arb(abs(v).upper()) * mult / 2)  # (1/2) from int_0^1 (1-t) dt
        rho_hi = arb(HR.rho_of(xa)).upper(); L = arb(xb) - arb(xa)
        IR += L * arb(rho_hi); IR2 += L * arb(rho_hi) ** 2; RWs = amax(RWs, arb(rho_hi)); ncell += 1
assert ncell == 3 * NS
tot = {}
for (cnt, lij), v in mono.items(): tot[cnt] = tot.get(cnt, arb(0)) + v       # sum of sups >= sup of the sum
NAME = lambda cnt: 'C3S_%d%d%d' % cnt
CONST = {}
for cnt in sorted(tot):
    r = tot[cnt]; val = VC.fup(r)
    assert val >= float(r.upper()) and arb(val) >= r.upper()
    CONST[NAME(cnt)] = dict(raw=r.str(radius=True, more=True), raw_upper=float(r.upper()), value=val, dir='up',
                            source=f'strip_cubic_v4c: sum over (l,i,j) of sup_strip (1/2) mult |d^2 C_l/dU_i dU_j|, monomial |w|^{cnt[0]} |w\'|^{cnt[1]} |w\'\'|^{cnt[2]}')
expected = {(0, 1, 2), (0, 2, 1), (0, 3, 0), (1, 0, 2), (1, 1, 1), (1, 2, 0), (2, 0, 1), (2, 1, 0), (3, 0, 0)}
assert set(tot) <= expected, set(tot) - expected
for cnt in expected - set(tot):            # absent monomial: coefficient 0 (record it explicitly so the verifier never falls back to a literal)
    CONST[NAME(cnt)] = dict(raw='0', raw_upper=0.0, value=0.0, dir='exact', source='no nonzero second derivative with this monomial')
# cross-check: the strip integrals of this script are covered by the ones the verifier uses (e0v4b_constants.json IR, IRW, RW from e0_strip_drift)
C4 = json.load(open(os.path.join(VC.V4B, 'v4cert', 'out', 'constants_v4.json')))['constants']
KW = arb(repr(C4['KW']['value'])); RW = RWs * KW; IRW = (IR2 * KW).sqrt()
E0 = json.load(open(E0C))['constants']
xchk = {n: dict(here=float(v.upper()), verifier=E0[n]['value'], covered=bool(arb(E0[n]['value']) >= v.upper()))
        for n, v in (('IR', IR), ('IRW', IRW), ('RW', RW))}
out = dict(constants=CONST, strip_integrals_crosscheck=xchk,
           box=dict(alpha_max=2.5e-4, w_jets_strip=0.02, strip_sigma=[29, 32], strip_x=[31.25, 34.25], cells=ncell),
           term_H8='RW [(C3S_030 sI1 + C3S_120 S0s) V0^2 + (C3S_021 sI1 + C3S_111 S0s) V0 V1 + (C3S_012 sI1 + C3S_102 S0s) V1^2] + IRW S0s^2 (C3S_201 V1 + C3S_210 V0) + C3S_300 IR S0s^3 '
                   '(state form, sup on the lowest-order factor; on the strip w\' = v: V0 = ||v||_W, V1 = ||v\'||_W; S0s = sup|w|, sI1 = sup|w\'| on [16.25,34.25])',
           review='fix11_review/REVIEW_v4b_v36_proofs.md section 4.3 (D-3); review values in fix11_review/d_cubic_strip.json')
REV = '/home/claude/s1spec/writeup/lemmaU/bulk/fix11_review/d_cubic_strip.json'
rv = json.load(open(REV))['strip_cubic_coefficients']
out['review_comparison'] = {NAME(eval(k_)): dict(review=v_, here=CONST[NAME(eval(k_))]['value'], rel=(CONST[NAME(eval(k_))]['value'] / v_ - 1) if v_ else None) for k_, v_ in rv.items()}
out['args'] = VC.args_block([PKL, E0C, REV, os.path.join(VC.V4B, 'v4cert', 'out', 'constants_v4.json'), STEP3NL + '/h2res.py', STEP3NL + '/h2fun.py', STEP3NL + '/psi3fun.py'] + VC.EXT,
                            extra=dict(runtime_s=round(time.time() - t0, 1)))
VC.dump(os.path.join(VC.OUT, 'strip_cubic_v4c.json'), out)
print(json.dumps({k_: v_['value'] for k_, v_ in CONST.items()}, indent=1))
print(json.dumps(xchk, indent=1)); print(json.dumps(out['review_comparison'], indent=1))
assert all(x['covered'] for x in xchk.values()), xchk
