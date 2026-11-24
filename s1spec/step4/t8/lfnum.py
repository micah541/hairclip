"""Lemma F constants (float) for one base sample, LEMMA_F.md section 3.  Identity: lfsym.py (lf_forms.pkl).
Per level j:  F_j' <= -lam_j F_j - th_j D_j + 2 sqrt(F_j) g_j + h_j  with
  g_j = sum_{i<j} (c_ji/2) sqrt(F_i) + |x| gam_j + Ybar a_j            (core couplings, modulation, e-terms)
  h_j = sum_{i<j} d_ji F_i + Ybar^2 e_j + [ramp data] + residual terms
  ramp data (non-nested level): k_jj V_j^2 + sum_{i<j} r_ji V_j V_i + sum_i k_ji V_i^2 + re_j V_j Ybar   (W-weighted E_I norms)
  ramp data (nested level, ramp_j inside the core of level j-2):  kD_j D_{j-2} + sum_i sD_ji sqrt(D_{j-2} F_i) + ...
spec: ramps {j: (ya, yb, qr)}, nested {j: j-2}, om (p, Yh, Yh_tau) and ynear (localized hand-off energies)."""
import pickle, numpy as np, sympy as sp
from math import comb
import lfsym as L
import weights as Wt
from fb import X0

PAR = dict(dP=0.25, d6=0.1, aS3=0.01, nu=1e-3)          # Young parameters (LEMMA_F.md 3.2)
LS = ((14.0, 28.0, 0.0),) * 5                  # Lemma S ramp for all levels (T5 design)

def _lamb(expr):
    subs = {}
    for m in range(6, 0, -1):
        subs[sp.Derivative(L.k, (L.x, m))] = sp.Symbol('k%d' % m)
        subs[sp.Derivative(L.T, (L.x, m))] = sp.Symbol('T%d' % m)
    e = sp.expand(expr).subs(subs).subs({L.k: sp.Symbol('k0'), L.T: sp.Symbol('T0'), L.Q: sp.Symbol('Qs')})
    args = [sp.Symbol('T%d' % m) for m in range(5)] + [sp.Symbol('k%d' % m) for m in range(6)] + [sp.Symbol('Qs'), L.c, L.x, L.X0]
    return sp.lambdify(args, e, 'numpy')

_F = None
def forms():
    global _F
    if _F is None:
        st = pickle.load(open('/home/claude/s1spec/step4/t8/lf_forms.pkl', 'rb')); _F = {}
        for j, pt in st.items():
            _F[j] = dict(P=_lamb(pt['P']), q=_lamb(pt['q']), P2=_lamb(L.Dn(pt['P'], 2)),
                         pairs={key: _lamb(cf) for key, cf in pt['pairs'].items()})
    return _F

def g_jets(s, mmax=5):
    K = s['K']; P = [s['p']]; Qq = [s['q']]; N = [s['n']]; T1 = [s['t1']]
    for m in range(mmax):
        P.append(-sum(comb(m, l) * K[l] * Qq[m - l] for l in range(m + 1)))
        Qq.append((1.0 if m == 0 else 0.0) + sum(comb(m, l) * K[l] * P[m - l] for l in range(m + 1)))
        N.append(-sum(comb(m, l) * K[l] * T1[m - l] for l in range(m + 1)))
        T1.append(sum(comb(m, l) * K[l] * N[m - l] for l in range(m + 1)))
    return P, N

def _mx(a, m): return float(np.max(a[m])) if np.any(m) else 0.0

def constants(s, ramps=LS, nested=None, om=None, ynear=None, jmax=4, Ke=1.4e4, muE=0.02, dP=None):
    nested = nested or {}
    dPj = (lambda j: PAR['dP']) if dP is None else (lambda j: dP[j] if isinstance(dP, (list, tuple, dict)) else dP)
    F = forms()
    X = s['X'].ravel(); Yy = X - X0; Wq = s['Wq'].ravel()
    yfar0 = min(r[0] for r in ramps[:jmax + 1])
    far = Yy >= yfar0 - 1e-12
    x = X[far]; y = Yy[far]; wq = Wq[far]
    K = [s['K'][m].ravel()[far] for m in range(6)]; Qs = s['Q'].ravel()[far]
    c = s['c']; ly = s['ell'] - X0; ellt = s['ellt']
    Wf = (1 + x**2) ** -2.0
    near = np.ones_like(y, bool) if ynear is None else (y <= ynear)
    Pd, Nd = g_jets(s)
    Pd = [a.ravel()[far] for a in Pd]; Nd = [a.ravel()[far] for a in Nd]
    gs = [[Pd[m] / s['np_'] for m in range(6)], [Nd[m] / s['nn_'] for m in range(6)],
          [(c * Pd[m] + s['b'] * Nd[m]) / s['nV'] for m in range(6)]]
    It = np.sqrt(np.cumsum(wq / Wf))                     # (int_{x*}^x W^-1)^(1/2)
    Tw = {j: Wt.weight(j, y, ly, ellt, ramps[j], 0.5, om) for j in range(jmax + 1)}
    # e-bound: |e| <= eps0 + H sqrt(F_0) on {eta_0 = 1}; |e| <= (Ke + It/sqrt(muE)) Ybar on the level-0 ramp
    yb0 = ramps[0][1]; C0 = y >= yb0
    T0 = Tw[0][0][0]
    Hc = np.sqrt(np.cumsum(np.where(C0, wq / np.where(T0 > 0, T0, np.inf), 0.0)))
    I1 = It[np.argmax(C0)] if np.any(C0) else It[-1]
    eps0 = Ke + I1 / np.sqrt(muE)
    er = Ke + It / np.sqrt(muE)
    out = dict(k=s['k'], tau=s['tau'], c=c, ell=s['ell'], lev={})
    for j in range(jmax + 1):
        T, dT = Tw[j]; Tj = T[0]
        ya, yb = ramps[j][0], ramps[j][1]
        R = (y > ya) & (y < yb); C = (y >= yb)
        args = [T[m] for m in range(5)] + K + [Qs, c, x, X0]
        ones = np.ones_like(x)
        q = F[j]['q'](*args) * ones; P = F[j]['P'](*args) * ones; P2 = F[j]['P2'](*args) * ones
        pairs = {key: f(*args) * ones for key, f in F[j]['pairs'].items()}
        Tpos = np.where(Tj > 0, Tj, np.inf)
        nS3 = sum(1 for (a, b) in pairs if a == j + 1) + 1          # +1: the residual/rough-modulation Young (LEMMA_F 3.3, S6)
        qh = dT + q + 0.5 * P2 + P**2 / (4 * dPj(j) * Tpos) + (PAR['aS3'] / 4) * nS3 * np.maximum(T[2], 0) + PAR['nu'] * np.abs(T[2])
        rate = qh / Tpos
        sig = c + _mx(rate, C) + nS3 * PAR['aS3'] / 8
        lam = c - sig
        rho = np.where(R, np.maximum(qh + lam * Tj, 0.0), 0.0)
        src = nested.get(j)
        lev = dict(sigma=sig, lam=lam, nested=src, k={}, kD=0.0, cji={}, cji_far={}, dji={}, rji={}, sDji={}, A1=0.0, ramp_e=0.0,
                   e_quad=0.0, Dcoef=2 - dPj(j) - PAR['aS3'] / 2 * nS3 - (PAR['d6'] if (j + 2, -1) in pairs else 0.0) - 0.5,
                   rate_core_max=_mx(rate, C))
        if src is None:
            lev['k'][j] = _mx(rho / Wf, R)
            if np.any(R): lev['k_at_y'] = float(y[R][np.argmax((rho / Wf)[R])])
        else:
            Ts = Tw[src][0][0]
            assert np.all(Ts[R] >= Ts.max() * 0 + Ts[R] * 0) and ramps[src][1] <= ya + 1e-9, 'nested ramp must lie in the core of its source'
            lev['kD'] = _mx(rho / np.where(Ts > 0, Ts, np.inf), R)
        Cn = C & near; Cf = C & ~near
        for (a, b), p in pairs.items():
            ap = np.abs(p)
            if a == j and b >= 0:
                Ti = Tw[b][0][0]
                lev['cji'][b] = lev['cji'].get(b, 0.0) + _mx(ap / np.sqrt(Tj * Ti + 1e-300), Cn)
                lev['cji_far'][b] = lev['cji_far'].get(b, 0.0) + _mx(ap / np.sqrt(Tj * Ti + 1e-300), Cf)
                if src is None: lev['rji'][b] = lev['rji'].get(b, 0.0) + _mx(ap / Wf, R)
                else:
                    Ts = Tw[src][0][0]
                    lev['sDji'][b] = lev['sDji'].get(b, 0.0) + _mx(ap / np.sqrt(Ts * Ti + 1e-300), R)
            elif a == j and b == -1:
                A1 = np.sqrt(np.sum(wq[C] * p[C]**2 / Tj[C])); A2 = np.sqrt(np.sum(wq[C] * p[C]**2 / Tj[C] * Hc[C]**2))
                lev['A1'] += A1 * eps0
                lev['cji'][0] = lev['cji'].get(0, 0.0) + 2 * A2
                if np.any(R):
                    if src is None:
                        lev['ramp_e'] += np.sqrt(_mx(ap / Wf, R) * np.sum(wq[R] * ap[R] * er[R]**2))
                    else:   # ramp_j in the level-0 core: e <= eps0 + H sqrt(F0): via D_src
                        Ts = Tw[src][0][0]
                        lev['ramp_e'] += np.sqrt(np.sum(wq[R] * p[R]**2 / Ts[R] * (eps0 + Hc[R])**2))   # crude (F0 <= 1 per unit)
            elif a == j + 1:
                if b == 0:
                    Tz = Tw[0][0][0]
                    lev['dji'][0] = lev['dji'].get(0, 0.0) + _mx(p**2 / (Tpos * np.where(Tz > 0, Tz, np.inf)), C) / (2 * PAR['aS3'])
                    lev['k'][0] = lev['k'].get(0, 0.0) + _mx(p**2 / (Tpos * Wf), R) / (2 * PAR['aS3'])
                else:
                    lev['e_quad'] += (np.sum(wq[C] * p[C]**2 / Tj[C]) * 2 * eps0**2 + np.sum(wq[R] * p[R]**2 / Tpos[R] * er[R]**2)) / (2 * PAR['aS3'])
                    lev['dji'][0] = lev['dji'].get(0, 0.0) + 2 * np.sum(wq[C] * p[C]**2 / Tj[C] * Hc[C]**2) / (2 * PAR['aS3'])
            elif a == j + 2:
                lev['e_quad'] += (np.sum(wq[C] * p[C]**2 / Tj[C]) * 2 * eps0**2 + np.sum(wq[R] * p[R]**2 / Tpos[R] * er[R]**2)) / (4 * PAR['d6'])
                lev['dji'][0] = lev['dji'].get(0, 0.0) + 2 * np.sum(wq[C] * p[C]**2 / Tj[C] * Hc[C]**2) / (4 * PAR['d6'])
        lev['gam'] = float(np.sqrt(sum(np.sum(wq * Tj * gs[i][j + 1]**2) for i in range(3))))
        lev['gam_cone'] = float(np.sqrt(sum(np.sum((wq * Tj * gs[i][j + 1]**2)[y <= 600.0]) for i in range(3))))
        grf = (-K[1] if j == 0 else -K[j + 1]) / s['nV']
        lev['gam_rough'] = float(np.sqrt(np.sum(wq * Tj * grf**2)))
        out['lev'][j] = lev
    out['gamma0_far'] = max(np.sqrt(np.sum(wq * Wf * gs[i][1]**2)) for i in range(3))
    out['gamma2_far'] = max(np.sqrt(np.sum(wq * Wf * gs[i][3]**2)) for i in range(3))
    return out

if __name__ == '__main__':
    import fb, sys
    i = int(sys.argv[1]) if len(sys.argv) > 1 else 150
    s = fb.sample(i, 0.5)
    for name, kw in [('T5 [14,28]', {}),
                     ('narrow [14,20]', dict(ramps=((14.0, 20.0, 0.0),) * 5)),
                     ('power q=2j+5 on [14,28]', dict(ramps=tuple((14.0, 28.0, 2.0 * j + 5) for j in range(5)))),
                     ('nested 3,4 <- D1,D2', dict(ramps=((14.0, 20.0, 0.0),) * 3 + ((20.0, 28.0, 0.0),) * 2, nested={3: 1, 4: 2}))]:
        o = constants(s, **kw)
        print(name, f"k={o['k']:.4g}")
        for j, lv in o['lev'].items():
            print(f"   j={j} lam={lv['lam']:.4f} k={ {m: f'{v:.3e}' for m, v in lv['k'].items()} } kD={lv['kD']:.3e} cji={ {a: f'{v:.1e}' for a, v in lv['cji'].items()} }"
                  f" sD={ {a: f'{v:.1e}' for a, v in lv['sDji'].items()} } gam={lv['gam']:.2e} gam_cone={lv['gam_cone']:.2e} Dcoef={lv['Dcoef']:.2f}")
