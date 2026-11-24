"""Independent FLOAT reference: design scripts (sB, sA, a_f, b_f, mu window, W(40), W(60)) + own sympy/numpy evaluation of the rest
(k=m=p=0, q = x + tip_mid, explicit sympy differentiation of the smoothsteps; no shared code with the rigorous scripts)."""
import sys, json, os, numpy as np, sympy as sp
os.chdir('/home/claude/s1spec/step3nl/design'); sys.path.insert(0, '.')
import coupling_design as CD, coupling_design2 as CD2
out = {}
MU1 = {2: (0.1552, 0.3357), 1.5: (0.1241, 0.2211)}
X = CD.X; STIP = CD.STIP
def S(u): return u**4 * (35 - 84*u + 70*u**2 - 20*u**3)
def float_extra(P):
    x, h = CD2.weight(P, 30.25, 38.25, xmax=70, h=1e-3)[0:2]
    xs, W, dphi = CD2.weight(P, 30.25, 38.25, xmax=70, h=5e-4)
    Wf = lambda xv: np.interp(xv, xs, W)
    rho = lambda xv: np.exp(-(1 + xv**2)**(2/3)/4)
    r = {}
    # strip: chi(sigma)
    sg = sp.symbols('sg'); chi = S((32 - sg)/3)
    cd = [sp.lambdify(sg, sp.diff(chi, sg, j), 'numpy') for j in range(5)]
    sig = np.linspace(29, 32, 6001); xx = sig - STIP
    c = [f(sig) + 0*sig for f in cd]; q = sig
    drift = q/4
    coef = [c[4] + drift*c[1], 4*c[3], 6*c[2], 4*c[1]]
    rw = np.sqrt(rho(xx)/Wf(xx))
    for j in (3, 2, 1, 0): r[f'd{j}'] = float(np.max(np.abs(coef[j])*rw))
    Ap = CD.ibp_coeffs(CD.rho, CD.smooth_down(X + STIP, 29, 32), CD.q)[0]
    Af = sp.lambdify(X, Ap, 'numpy')(xx)
    r['IA'] = float(np.trapezoid(np.maximum(Af, 0), sig)); r['K_rho_16.25'] = float(1/rho(16.25))
    xi = np.arange(16.25, 34.25 + 1e-9, 5e-4); r['Phi_34.25'] = float(np.trapezoid(1/Wf(xi), xi))
    # inner ramp: zeta
    xr = np.linspace(10.25, 16.25, 6001)
    Wp = lambda xv: (1 + xv**2)**(-P)
    zs = 1 - S((16.25 - X)/6); zc = [sp.lambdify(X, sp.diff(zs, X, j), 'numpy') for j in range(5)]
    z = [f(xr) + 0*xr for f in zc]; qr = xr + STIP
    cf = [z[4] + qr/4*z[1], 4*z[3], 6*z[2], 4*z[1]]
    for j in (3, 2, 1, 0): r[f'c{j}'] = float(np.max(np.abs(cf[j])*np.sqrt(Wp(xr)/rho(xr))))
    r['drift_f'] = float(np.max(qr/4*z[0]*z[1]*Wp(xr)/rho(xr)))
    # outer
    Wsym = sp.Symbol('Wsym')
    def AB_float(dphi, cut):
        Dx = lambda e: sp.diff(e, X) + sp.diff(e, Wsym) * (-dphi * Wsym)
        def Dn(e, n):
            for _ in range(n): e = Dx(e)
            return e
        cdv = [sp.diff(cut, X, j) for j in range(5)]
        g1, g2, g3, g4 = [Wsym * cut * cdv[j] for j in (1, 2, 3, 4)]
        A = -2*Dn(g1, 3) + 3*Dn(g2, 2) - 2*Dn(g3, 1) + g4 + (X + STIP)/4*g1
        B = 6*Dn(g1, 1) - 6*g2
        return sp.lambdify((X, Wsym), A, 'numpy'), sp.lambdify((X, Wsym), B, 'numpy')
    fr = X*(1 + X**2)**sp.Rational(-1, 3)/3; fp = 2*sp.nsimplify(P)*X/(1 + X**2); sb = S((X - sp.Rational(121, 4))/8)
    for R in (60, 40):
        xo = np.linspace(R - 3, R, 6001); zo = S((R - X)/3)
        zc = [sp.lambdify(X, sp.diff(zo, X, j), 'numpy') for j in range(5)]
        z = [f(xo) + 0*xo for f in zc]; qo = xo + STIP
        Aa, Ba = AB_float((1 - sb)*fp + sb*fr, zo); Ar, Br = AB_float(fr, zo)
        msk = xo <= 38.25; Wv = Wf(xo)
        Av = np.where(msk, Aa(xo, Wv), Ar(xo, Wv)); Bv = np.where(msk, Ba(xo, Wv), Br(xo, Wv))
        r[f'Azo_{R}'] = float(np.max(np.abs(Av))); r[f'Bzo_{R}'] = float(np.max(np.abs(Bv)))
        r[f'W{R-3}'] = float(Wf(R - 3)); r[f'W{R}'] = float(Wf(R))
        cf = [z[4] + qo/4*z[1], 4*z[3], 6*z[2], 4*z[1]]
        for j in (3, 2, 1, 0): r[f'oc{j}_{R}'] = float(np.max(np.abs(cf[j])*np.sqrt(Wv)))
    return r
for P in (2, 1.5):
    m1, m1s = MU1[P]
    d = CD2.evaluate(P, 10.25, 16.25, m1, m1s)
    keep = {k: float(v) for k, v in d.items()}
    keep.update(float_extra(P)); out[str(P)] = keep
    print(P, keep, flush=True)
json.dump(out, open('/home/claude/s1spec/step3nl/agent_coup/float_ref.json', 'w'), indent=1)
