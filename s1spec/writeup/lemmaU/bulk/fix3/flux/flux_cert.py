#!/usr/bin/env python3
"""F1 flux lemma, step 2 (certified, python-flint arb, outward rounding): coefficient tables of the reduced forms of flux_reduce.py
for the ACTUAL verifier cutoff zeta = S7((60-x)/3), the verifier weight W = exp(-psi), psi' = phi' = x/(3(1+x^2)^{1/3}) on x >= 38.25,
and the linear operator L1 of the verifier (profile k, m, p, q from the step-1 tail theorem, step3/hstail.profile_bounds).
No jet box: the commutator is that of the LINEAR operator L1 (quasilinear parts are in the verifier's N-terms, multiplied by zeta).

For every reduced pair (a,b): Gam_ab = sup_{x in [57,60]} |c_ab(x)|/W(x) (a != b) or sup max(c_aa/W, 0) (a == b, one-sided: only an
upper bound of K is needed).  kappa_l = sum of diagonal Gam_ll + half of every off-diagonal Gam_ab with l in {a,b} (as v3.1).
Then  K <= sum_l kappa_l N_l^2,  N_l^2 = int_57^60 W (g^(l))^2 dx  (the same N_l as lemmaS_cflux_final.py / assembly_arb.f4b).
Transfer to the certified v35 flux pipeline (asm_design_N2ref4.json, beta_A = 2, with d and chart terms): the pipeline value is a
monotone, degree-1 homogeneous function of the vector (kappa_0, kappa_1, kappa_2), so K_new <= t * K_old with t = max_l kappa_new_l/kappa_old_l.
Writes flux_cert.json (arb strings, doubles rounded outward, labels, args block)."""
import sys, os, json, time, pickle, hashlib, platform, math
import sympy as sp
from flint import arb, ctx
ctx.prec = 200
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import flux_reduce as FR
sys.path.insert(0, '/home/claude/s1spec/step3')
import hstail
x = FR.x
NCELL = int(os.environ.get('NCELL', '3000'))
NODRIFT = os.environ.get('NODRIFT', '0') == '1'      # diagnostic (not a bound): drop C = 7km + q/4
OUTNAME = 'flux_cert_nodrift_DIAG.json' if NODRIFT else 'flux_cert.json'
OLD_CFLUX = '/home/claude/s1spec/writeup/lemmaS/lemmaS_cflux.json'
ASM = '/home/claude/s1spec/writeup/lemmaU/bulk/fix2/n12/t2_A_v35/out/asm_design_N2ref4.json'
VCONST = '/home/claude/s1spec/step3nl_v33/verifier/run_v33.sh'

class _Pr(sp.printing.str.StrPrinter):
    def _print_Rational(self, e): return "R(%d,%d)" % (e.p, e.q)
    def _print_Integer(self, e): return "arb(%d)" % e.p
    def _print_Pow(self, e, rational=False):
        if e.exp.is_Integer: return "IP(%s, %d)" % (self._print(e.base), int(e.exp))     # integer powers by repeated products (balls containing 0)
        raise ValueError("non-integer power in compiled expression: %s" % e)
def R(p, q): return arb(p)/arb(q)
def IP(b, n):
    """b**n for integer n by repeated multiplication (arb ** int is evaluated as exp(n log b) in python-flint: NaN for balls containing 0)"""
    r = arb(1)
    for _ in range(abs(n)): r = r*b
    return r if n >= 0 else 1/r
def compile_expr(e, names):
    src = _Pr().doprint(e)
    code = compile(src, '<expr>', 'eval')
    return lambda env: eval(code, {'arb': arb, 'R': R, 'IP': IP}, env)

# ---- W ratios r_n = W^(n)/W in terms of x and S = (1+x^2)^(1/3) (exact) ----
psi = sp.Function('psi')(x)
S = sp.Symbol('S', positive=True)
phi = (1 + x**2)**sp.Rational(2, 3)/4
def wratio(n):
    e = sp.diff(sp.exp(-psi), x, n)/sp.exp(-psi)
    e = sp.simplify(e)
    for j in range(n, 0, -1):
        e = e.subs(sp.Derivative(psi, (x, j)), sp.diff(phi, x, j))
    e = sp.expand(e)
    e = e.replace(lambda t: t.is_Pow and t.base == 1 + x**2 and t.exp.is_Rational and (t.exp*3).is_Integer, lambda t: S**int(t.exp*3))
    assert not e.has(sp.Pow(1 + x**2, sp.Rational(1, 3))), e
    return e

# ---- profile ODE derivatives of C = 7km + q/4 and K = k^2 ----
kS, mS, pS, qS = sp.symbols('k m p q')
def dode(f): return sp.expand(sp.diff(f, kS)*mS + sp.diff(f, mS)*pS/4 + sp.diff(f, pS)*(-kS*qS) + sp.diff(f, qS)*(1 + kS*pS))
def ode_derivs(f, n):
    out = [f]
    for _ in range(n): out.append(dode(out[-1]))
    return out

uS = sp.Symbol('u')
S7 = uS**4*(35 - 84*uS + 70*uS**2 - 20*uS**3)

def ball(lo, hi): return arb((lo + hi)/2, (hi - lo)/2)
def up(a): return float(arb(a).upper()) if not isinstance(a, float) else a
def fup(a):
    v = float(a.upper()); return math.nextafter(v, math.inf) if v != 0 else (0.0 if a.upper() <= 0 else math.nextafter(0.0, math.inf))

def sha(p): return hashlib.sha256(open(p, 'rb').read()).hexdigest()

def load_arb(j):
    m, e, r = j['arb']
    return arb(int(m))*arb(2)**int(e) + arb(0, r) if r else arb(int(m))*arb(2)**int(e)

def main():
    t0 = time.time()
    forms = pickle.load(open(os.path.join(HERE, 'flux_forms.pkl'), 'rb'))
    Wf, Zf, Cf, Kf = FR.Wf, FR.Zf, FR.Cf, FR.Kf
    # symbols for evaluation
    nW = 6; nZ = 8; nC = 6
    rsym = sp.symbols('r0:%d' % nW); zsym = sp.symbols('z0:%d' % nZ); csym = sp.symbols('c0:%d' % nC); ksym = sp.symbols('kk0:%d' % nC)
    def to_eval(e):
        rep = {}
        for d in e.atoms(sp.Derivative):
            n = d.derivative_count
            rep[d] = {Wf: rsym, Zf: zsym, Cf: csym, Kf: ksym}[d.expr][n]
        e = e.xreplace(rep).xreplace({Wf: rsym[0], Zf: zsym[0], Cf: csym[0], Kf: ksym[0]})
        return sp.expand(e)
    # W ratios, C^(j), K^(j)
    WR = [compile_expr(wratio(n), None) for n in range(nW)]
    Cd = [compile_expr(e, None) for e in ode_derivs(7*kS*mS + qS/4, nC - 1)]
    Kd = [compile_expr(e, None) for e in ode_derivs(kS**2, nC - 1)]
    S7d = [compile_expr(sp.horner(sp.expand(sp.diff(S7, uS, i)), uS), None) for i in range(nZ)]   # S7^(nZ) = 0 for nZ = 8
    stip = hstail.STIP
    s1 = arb(57) + stip; s1lo = arb(float(s1.lower())); s1lo = arb(math.nextafter(float(s1lo.lower()), -math.inf))
    kb, mb, pb, qs = hstail.profile_bounds(s1lo)
    kb, mb, pb, qs = [arb(0, a.upper()) for a in (kb, mb, pb, qs)]
    prof = dict(s1=float(s1lo.mid()), kb=up(abs(kb).upper()), mb=up(abs(mb).upper()), pb=up(abs(pb).upper()), qs=up(abs(qs).upper()))
    print('profile bounds on sigma >= %.6f:' % prof['s1'], prof, flush=True)
    # compile pair coefficients
    comp = {}
    for lev in (0, 2):
        for key, srep in forms[lev].items():
            e = to_eval(eval(srep, {**sp.__dict__, 'x': x}))
            comp[(lev, key)] = compile_expr(e, None)
    # one extra diagnostic: the obstruction coefficient zeta zeta''''  (W-free) and sup|zeta''''|
    obs_c = compile_expr(zsym[0]*zsym[4], None); z4c = compile_expr(zsym[4], None)
    sups = {k: arb(0) for k in comp}; sup_obs = arb(0); sup_z4 = arb(0); sup_c = [arb(0)]*nC
    a0, a1 = arb(57), arb(60)
    for i in range(NCELL):
        lo = a0 + (a1 - a0)*i/NCELL; hi = a0 + (a1 - a0)*(i + 1)/NCELL
        X = ball(lo, hi)
        Sx = (1 + X*X)**(arb(1)/3)
        env = {'x': X, 'S': Sx}
        rv = [f(env) for f in WR]
        u = (arb(60) - X)/3
        ulo = max(0.0, float(u.lower())); uhi = min(1.0, float(u.upper()))
        uL = arb(math.nextafter(ulo, -math.inf)) if ulo > 0 else arb(0)      # outward, clipped to [0, 1]
        uH = arb(math.nextafter(uhi, math.inf)) if uhi < 1 else arb(1)
        ub = uL.union(uH)
        # (ub may exceed [0,1] by a radius-rounding ulp; it only has to CONTAIN the true u-range)
        assert ub.lower() <= max(u.lower(), arb(0)).lower() and ub.upper() >= min(u.upper(), arb(1)).upper()   # ub contains u clipped to [0, 1]
        uc = arb(ub.mid())                                            # mean-value form on the cell (S7 is a polynomial; S7^(8) = 0)
        zv = []
        for j in range(nZ):
            sj = S7d[j]({'u': uc}) + (S7d[j + 1]({'u': ub})*(ub - uc) if j + 1 < len(S7d) else arb(0))
            zv.append((arb(-1)/3)**j*sj)
        qb = X + stip + qs           # sigma = x + sigma_tip (STIP < 0), q in sigma +- qs
        penv = {'k': kb, 'm': mb, 'p': pb, 'q': qb}
        cv = [f(penv) for f in Cd]; kv = [f(penv) for f in Kd]
        if NODRIFT: cv = [arb(0)]*nC            # diagnostic only: drift commutator switched off
        env2 = {}
        for j in range(nW): env2['r%d' % j] = rv[j]
        for j in range(nZ): env2['z%d' % j] = zv[j]
        for j in range(nC): env2['c%d' % j] = cv[j]; env2['kk%d' % j] = kv[j]
        for (lev, key), f in comp.items():
            val = f(env2)
            a, b = map(int, key.split(','))
            if a == b:
                ubv = val.upper()
                cand = arb(ubv) if ubv > 0 else arb(0)
            else:
                cand = arb(abs(val).upper())
            if cand.upper() > sups[(lev, key)].upper(): sups[(lev, key)] = cand
        o = arb(abs(obs_c(env2)).upper()); sup_obs = o if o.upper() > sup_obs.upper() else sup_obs
        o = arb(abs(z4c(env2)).upper()); sup_z4 = o if o.upper() > sup_z4.upper() else sup_z4
        for j in range(nC):
            o = arb(abs(cv[j]).upper()); sup_c[j] = o if o.upper() > sup_c[j].upper() else sup_c[j]
        if i % 500 == 0: print('  cell', i, '/', NCELL, '%.0fs' % (time.time() - t0), flush=True)
    # sign facts used: C >= 0 on the ramp (drift), zeta zeta' <= 0
    # kappa
    def kap(lev, maxl, exclude=()):
        k = [arb(0)]*(maxl + 1)
        for (lv, key), s in sups.items():
            if lv != lev or key in exclude: continue
            a, b = map(int, key.split(','))
            if max(a, b) > maxl: continue
            if a == b: k[a] += s
            else: k[a] += s/2; k[b] += s/2
        return k
    OB_KEYS = ('4,0', '3,1', '3,0')
    old = json.load(open(OLD_CFLUX))
    ko0 = [arb(math.nextafter(v, -math.inf)) for v in old['kappa_I0']]; ko2 = [arb(math.nextafter(v, -math.inf)) for v in old['kappa_I2']]
    def kap_w(lev, lam):
        """weighted AM-GM: |c g_a g_b| <= (|c|/2)(lam_ab g_a^2 + g_b^2/lam_ab), a > b; lam exact rationals (any lam > 0 is valid)"""
        k = [arb(0)]*3
        for (lv, key), s_ in sups.items():
            if lv != lev: continue
            a, b = map(int, key.split(','))
            if max(a, b) > 2: continue
            if a == b: k[a] += s_
            else:
                L = lam.get(key, arb(1)); k[a] += s_*L/2; k[b] += s_/(2*L)
        return k
    def choose_lam(lev, kold):
        """float grid search for the weights minimising max_l kappa_l/kold_l (choice only; the bound is recomputed in arb)"""
        from fractions import Fraction
        G_ = {k: float(v.upper()) for (lv, k), v in sups.items() if lv == lev}
        ko = [float(z.mid()) for z in kold]
        grid = [Fraction(n, 8) for n in range(1, 8*30 + 1)]
        grid = [g_ for g_ in grid if g_ <= 4 or g_.denominator <= 2]
        best = None
        keys = [k for k in ('1,0', '2,0', '2,1') if k in G_]
        import itertools
        coarse = [Fraction(n, 4) for n in (1, 2, 3, 4, 6, 8, 10, 12, 16, 20, 24, 32, 40, 48, 64, 80)] + [Fraction(1, 2**j) for j in (3, 4, 5)]
        for combo in itertools.product(coarse, repeat=len(keys)):
            lam = dict(zip(keys, combo)); k = [G_.get('%d,%d' % (l, l), 0.0) for l in range(3)]
            for key in keys:
                a, b = map(int, key.split(',')); L = float(lam[key]); k[a] += G_[key]*L/2; k[b] += G_[key]/(2*L)
            sc = max(k[l]/ko[l] for l in range(3))
            if best is None or sc < best[0]: best = (sc, lam)
        return {k: arb(v.numerator)/arb(v.denominator) for k, v in best[1].items()}, {k: str(v) for k, v in best[1].items()}
    lam0, lam0s = choose_lam(0, ko0); lam2, lam2s = choose_lam(2, ko2)
    k0 = kap_w(0, lam0)
    k2 = kap_w(2, lam2)
    k2_3 = sups[(2, '3,3')]
    ob = {key: sups[(2, key)] for key in OB_KEYS}
    t0r = max((k0[l]/ko0[l] for l in range(3)), key=lambda z: z.upper())
    t2r = max((k2[l]/ko2[l] for l in range(3)), key=lambda z: z.upper())
    A = json.load(open(ASM))['res']['tight']['flux']
    K0old = load_arb(A['beta']['2']['K0']); K2old = load_arb(A['beta']['2']['K2'])
    FLUX0 = arb('1.64e-20'); FLUX2 = arb('8.0e-20')
    K0new = t0r*K0old; K2pw = t2r*K2old
    res = dict(
        profile=prof, ncell=NCELL,
        Gamma={'%d:%s' % k: str(v) for k, v in sups.items()},
        Gamma_up={'%d:%s' % k: fup(v) for k, v in sups.items()},
        kappa0_new=[fup(z) for z in k0], amgm_weights_level0=lam0s, amgm_weights_level2=lam2s, kappa2_new_le2=[fup(z) for z in k2], kappa2_level3_33=fup(k2_3),
        obstruction_pairs_sup={k: fup(v) for k, v in ob.items()},
        sup_zeta_zeta4=fup(sup_obs), sup_abs_zeta4=fup(sup_z4), sup_abs_Cderiv=[fup(z) for z in sup_c],
        kappa0_old=old['kappa_I0'], kappa2_old=old['kappa_I2'],
        ratio_t0=fup(t0r), ratio_t2=fup(t2r),
        K0_old_beta2=fup(K0old), K2_old_beta2=fup(K2old),
        K0_new_beta2=fup(K0new), K2_pointwise_new_beta2=fup(K2pw),
        margin_FLUX0=float((FLUX0/K0new).lower()), margin_FLUX2_pointwise_part=float((FLUX2/K2pw).lower()),
        labels=dict(Gamma='cert', kappa='cert', K_new='cert (transfer from the certified v35 pipeline, monotone/homogeneous)', margins='cert (lower bounds)'))
    ins = [os.path.abspath(__file__), os.path.join(HERE, 'flux_reduce.py'), os.path.join(HERE, 'flux_forms.pkl'), hstail.__file__,
           '/home/claude/s1spec/lin/lintail2.py', OLD_CFLUX, ASM, VCONST]
    args = dict(argv=sys.argv, cwd=os.getcwd(), python=sys.version, platform=platform.platform(), flint=__import__('flint').__version__, sympy=sp.__version__,
                env={k: os.environ.get(k) for k in ('PYTHONHASHSEED', 'NCELL', 'NODRIFT')}, inputs={p: sha(p) for p in ins}, prec=ctx.prec,
                runtime_s=round(time.time() - t0, 1))
    res['args'] = args
    json.dump(res, open(os.path.join(HERE, OUTNAME), 'w'), indent=1)
    for k in ('kappa0_new', 'kappa0_old', 'kappa2_new_le2', 'kappa2_old', 'kappa2_level3_33', 'obstruction_pairs_sup', 'sup_zeta_zeta4', 'sup_abs_zeta4',
              'ratio_t0', 'ratio_t2', 'K0_new_beta2', 'K2_pointwise_new_beta2', 'margin_FLUX0', 'margin_FLUX2_pointwise_part'):
        print(k, res[k])
    print('done %.1fs' % (time.time() - t0))

if __name__ == '__main__':
    main()
