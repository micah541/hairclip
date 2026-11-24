#!/usr/bin/env python3
"""F1: Lemma S(A) in divergence form (plan STEP3_FIXPLAN F1).  Core y in [28,60], inner ramp [14,28] (s-length ln 2), outer cutoff far away.
Machine-verified identities (sympy) + arb sups.  Ramp commutators are linear in eta-derivatives and are bounded by Cauchy-Schwarz against J_k (k<=j+1);
NO Young against the shell dissipation on the ramp.  Young on the core only (same weight Theta as the dissipation)."""
import pickle, json, math, sys, itertools
import sympy as sp
from flint import arb, ctx
ctx.prec = 100
sys.path.insert(0, '.')
import lemmaS_prop1_sym as P
import lemmaS_prop1_bounds as PB
s = P.s
res = pickle.load(open('prop1_forms.pkl', 'rb'))
X0 = 2.2490776515
STIR = {0: {0: 1}, 1: {1: 1}, 2: {1: 1, 2: 1}, 3: {1: 1, 2: 3, 3: 1}}
BELL = {k: sum(v.values()) for k, v in STIR.items()}

# ---------------- A: divergence form (one IBP of the cross terms with z_{j+2}) ----------------
def div_form(j):
    Q = dict(res[j][0]); n = j + 2; out = {}
    def add(k, v): k = tuple(sorted(k)); out[k] = out.get(k, 0) + v
    for (a, b), c in Q.items():
        if b < n or (a, b) == (n, n): add((a, b), c)
        else:                                       # b == n, a < n
            if a == n - 1: add((a, a), -sp.Rational(1, 2)*P.D(c))
            else: add((a + 1, n - 1), -c); add((a, n - 1), -P.D(c))
    out = {k: sp.expand(v) for k, v in out.items() if v != 0}
    ok = P.is_total_derivative(sp.expand(P.to_u(P.form_expr(Q)) - P.to_u(P.form_expr(out))))
    return out, ok

# ---------------- B: numeric helpers ----------------
tt = sp.Symbol('t')
Sx = tt**3*(10 - 15*tt + 6*tt**2)
def sup_grid(f, n=1500):
    """rigorous upper bound of f(t), t in [0,1], via n ball sub-intervals (f: arb -> arb)"""
    best = 0.0
    for i in range(n):
        lo = arb(i)/n; hi = arb(i + 1)/n
        v = f(arb((lo + hi)/2, (hi - lo)/2)); u = float(abs(v).upper())
        if u != u: u = 1e300
        best = max(best, u)
    return best

class Ramp:
    """y = y0 * (y1/y0)^t ; eta = S(t) (increasing=True) or 1-S(t); omega = (y/28)^-gam; Wf = closure weight"""
    def __init__(self, y0, y1, gam, increasing=True, yref=28.0):
        self.y0, self.y1, self.gam, self.inc = y0, y1, gam, increasing
        self.L = math.log(y1/y0); self.yref = yref
        e = Sx**2 if increasing else (1 - Sx)**2
        self.eta2 = [sp.lambdify(tt, sp.diff(e, tt, l)/self.L**l, modules=[{}]) for l in range(6)]
    def y(self, t): return arb(self.y0)*(t*arb(self.L)).exp()
    def om(self, t): return (-(arb(self.gam))*((self.y(t)/arb(self.yref)).log())).exp()
    def eps(self, t): return self.y(t)**-4
    def W(self, t): return (1 + (self.y(t) + arb(str(X0)))**2)**-2
    def e2(self, t, l): return arb(self.eta2[l](t))
    def rem(self, t, m):          # |omega_E^(m) - (-gam)^m omega_E| / 1   (omega_E = eta^2 omega)
        return abs(self.om(t)*sum(math.comb(m, l)*self.e2(t, l)*(-self.gam)**(m - l) for l in range(1, m + 1)))
    def rho(self, t):             # (1/2)|Theta'' - kappa^2 Theta| with Theta = c0 e^{-kappa s} eta^2 ; here returned divided by eps*omega
        kap = self.gam + 4
        return abs((self.e2(t, 2) - 2*kap*self.e2(t, 1))*self.om(t)*self.eps(t))/2

# ---------------- C: assemble constants ----------------
def build_tables():
    T = {}
    for j in (0, 1, 2):
        Qd, ok = div_form(j); assert ok, j
        T[j] = PB.table(Qd)
    return T

def A_bound_set(mode, box):
    if mode == 'lin':
        Ab, Fa, Fb = PB.A_bounds_s(box, rmax=4)
        return {k: (v if k[1] == 0 else 0.0) for k, v in Ab.items()}
    Ab, Fa, Fb = PB.A_bounds_s(box, rmax=4)
    return Ab

def constants(gam, box, mode, y0=14.0, y1=28.0, tables=None):
    R = Ramp(y0, y1, gam)
    T = tables or build_tables()
    Ab = A_bound_set(mode, box)
    amin = (1 + box['M0']**2)**-2
    kap = gam + 4
    eps_max = y0**-4
    out = dict(amin=amin, eps_max=eps_max, mode=mode)
    G = {}
    for j in (0, 1, 2):
        n = j + 2; tab = T[j]
        YC = {k: 0.0 for k in range(n + 1)}
        pairsR = []        # (a,b, list of (m,i,r,coef))
        for (a, b), mons in tab.items():
            if (a, b) == (n, n): continue
            gC = sum(co*gam**m*Ab[(i, r)] for (m, i, r), co in mons.items())
            if a == b: YC[a] += gC
            else: YC[a] += gC/2; YC[b] += gC/2
            pairsR.append((a, b, mons))
        # ramp remainder coefficient functions, per k (<= j+1):  eps(t) * sum_pairs sum_mons coef * rem_m(t) * A_r   (+ interpolation remainder)
        Bc = 1.5*amin
        Ghat = {k: YC[k] for k in range(j + 1)}
        Ghat[j] += YC[j + 1]*(YC[j + 1]/(4*Bc) + kap**2/2)
        # ramp coefficient sups against data weights: need, for each k<=j+1 and m<=k:  sup eps*Yrem_k(t) * y^{2m-1}/W
        ramp_k = {}
        for k in range(j + 2):
            def Yrem(t, k=k):
                tot = arb(0)
                for (a, b, mons) in pairsR:
                    w = 0.0
                    if a == k: w += 1.0 if a == b else 0.5
                    elif b == k: w += 0.5
                    if w == 0.0: continue
                    for (m, i, r), co in mons.items():
                        if m == 0: continue
                        tot += arb(w*co*Ab[(i, r)])*R.rem(t, m)
                tot = tot*R.eps(t)
                if k == j: tot += arb(YC[j + 1])*R.rho(t)         # interpolation remainder (z_j^2)
                if k <= j: tot += arb(Ghat[k])*R.eps(t)*R.e2(t, 0)*R.om(t)   # ramp share of the core rate (eps(14) instead of eps(28))
                return tot
            ramp_k[k] = {m: sup_grid(lambda t, m=m: Yrem(t)*R.y(t)**(2*m - 1)/R.W(t)) for m in STIR[k]}
        # inflow  c * int (eta^2)'_s omega z_j^2
        infl = {m: sup_grid(lambda t, m=m: (R.e2(t, 1)*R.om(t)).__abs__()*R.y(t)**(2*m - 1)/R.W(t)) for m in STIR[j]}   # per ds: (eta^2)_s omega z^2 ds
        G[j] = dict(Ghat=Ghat, ramp=ramp_k, infl=infl, YC=YC)
        out[f'Ghat_{j}'] = {str(k): v for k, v in Ghat.items()}
        out[f'YC_{j}'] = {str(k): v for k, v in YC.items()}
        out[f'ramp_{j}'] = {str(k): {str(m): v for m, v in d.items()} for k, d in ramp_k.items()}
        out[f'infl_{j}'] = {str(m): v for m, v in infl.items()}
    return out, G

def assemble(G, J, sqLam, gam, alpha1, handoff=(0, 0, 0), y0=14.0, y1=28.0, YO=60.0, Rout=350.0, extra_outer=(0, 0, 0), slopes=(0, 0, 0, 0), abar_c=None, beta_c=None):
    # abar_c, beta_c: certified verifier sups of |(a,b)| and |beta| (Nov 10); if given they replace the legacy formulas
    J0, J1, J2, J3 = J
    Jl = [J0, J1, J2, J3]
    Abar = sqLam; abar = 7.57*Abar if abar_c is None else abar_c; bbar = abar
    cbar = 0.25 + abar + bbar/y0; cmin = 0.25 - abar - bbar/y0
    eps_max = y0**-4; eps_core = y1**-4; w_a = (28.0/y0)**gam
    E = []; parts = []
    p = 1.0246
    for j in range(3):
        g = G[j]
        gp = gam*cmin - (abs(1 - 2*j) + 1)*bbar/y0 - eps_core*g['Ghat'][j]
        Je = [Jl[m] + slopes[m]/gp for m in range(4)]      # variation of constants: int e^{-gp(tau-s)} dI <= J + slope/gp (uniform in tau)
        infl = cbar*BELL[j]*sum(STIR[j][m]*g['infl'][m]*Je[m] for m in STIR[j])
        rampc = sum(BELL[k]*sum(STIR[k][m]*g['ramp'][k][m]*Je[m] for m in STIR[k]) for k in range(j + 2))
        coup = sum(g['Ghat'].get(k, 0.0)*eps_core*E[k] for k in range(j))
        beta = (0.00244 + 63.5*alpha1)*Abar if beta_c is None else beta_c
        r0 = 1.4*(beta + 1.02*alpha1*abar)
        e = 2*p - gam
        Iint = 28.0**gam*((Rout**e - y0**e)/e if abs(e) > 1e-9 else math.log(Rout/y0))
        F2 = (r0*p**j)**2*Iint
        Ej = handoff[j] + infl + rampc + 4*F2/gp**2 + 2*coup/gp + extra_outer[j]
        E.append(Ej); parts.append(dict(gamma_prime=gp, inflow=infl, ramp_comm=rampc, forcing=4*F2/gp**2, coupling=2*coup/gp, handoff=handoff[j]))
    conv = (YO/28.0)**gam; ell = math.log(YO/28.0)
    Eu = [conv*e_ for e_ in E]
    S0 = math.sqrt(Eu[0]/ell + 2*math.sqrt(Eu[0]*Eu[1]))
    S1 = math.sqrt(Eu[1]/ell + 2*math.sqrt(Eu[1]*(Eu[1] + Eu[2])))
    return dict(E=E, Eu=Eu, S0=S0, S1=S1, parts=parts)

if __name__ == '__main__':
    for j in (0, 1, 2):
        Qd, ok = div_form(j)
        print('level', j, 'divergence form: pairs', sorted(Qd.keys()), 'total-derivative check', ok, flush=True)
