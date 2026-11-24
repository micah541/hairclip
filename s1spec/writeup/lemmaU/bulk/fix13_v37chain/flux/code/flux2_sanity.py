#!/usr/bin/env python3
"""Stage X6 (FLUX2_DESIGN.md sec. 10): FLOAT sanity checks, NOT part of the proof.  mpmath, 30 digits.
(1) For 3 test triples (g, j, C_eff) on [57, 60] with |j - 1| <= CBAR_OR, |j'| <= C1_OR, C_eff = (sigma)/4 (k = 0), W = exp(-(phi - phi(57))):
    (a) the defining integral  int W v (j K4g + j' K3g + C_eff zeta' g)  equals the reduced level-0 form of X1 (same j, j' = dj/dx);
    (b) it is <= sum_l kappa0_l N_l^2 (X2 kappa0);
    (c) |int W D2v X0| <= frak D frak X with the X2 coefficients (linear part only: frak D = sum frak a_i N_i + N_4, frak X from the j, j', C_eff rows).
(2) Reflected integral (sec. 7.2): a random phi >= 0 (spikes) with window integrals <= C' + S' T on a grid, s > S';
    Z = P - running min P satisfies 0 <= Z <= 2C' and Z' >= 2 (phi - s) (checked on the grid).
Output out/flux2_sanity.json."""
import sys, os, json, pickle, random
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import f2common as FC
import mpmath as mp
import sympy as sp
mp.mp.dps = 30
T_ = FC.Timer()
X2 = json.load(open(os.path.join(FC.OUT, 'flux2_coef.json')))['res']
forms = pickle.load(open(os.path.join(FC.OUT, 'flux2_forms.pkl'), 'rb'))
CB = X2['CBAR_OR']['up']; C1 = X2['C1_OR']['up']; kap = [r['up'] for r in X2['kappa0']]
aa = {int(k): v['up'] for k, v in X2['frak_a'].items()}; z = {int(k): v['up'] for k, v in X2['z'].items()}
STIP = mp.mpf('-2.249073')
x = sp.Symbol('x', real=True)
Wf, Zf, Cf, Kf, Jf = [sp.Function(n)(x) for n in ('W', 'Z', 'C', 'K', 'J')]
G = sp.symbols('G0:8')
ns = sp.__dict__.copy(); ns['x'] = x
L0 = {tuple(map(int, k.split(','))): eval(s, ns) for k, s in forms['level0'].items()}
u_ = sp.Symbol('u'); S7 = u_**4 * (35 - 84 * u_ + 70 * u_**2 - 20 * u_**3)
zeta = sp.lambdify(x, S7.subs(u_, (60 - x) / 3), 'mpmath')
zd = [sp.lambdify(x, sp.diff(S7.subs(u_, (60 - x) / 3), x, i), 'mpmath') for i in range(5)]
phi = lambda t: (1 + t * t) ** (mp.mpf(2) / 3) / 4
W = lambda t: mp.e ** (-(phi(t) - phi(57)))
psi1 = lambda t: t / (3 * (1 + t * t) ** (mp.mpf(1) / 3))
psi2 = lambda t: mp.diff(lambda s: mp.e ** (-phi(s)), t, 2) / mp.e ** (-phi(t))
TESTS = [dict(g=lambda t: mp.sin(2 * t) * 1e-3 + 1e-4 * t, j=lambda t: 1 + CB * 0.9 * mp.cos(t / 7), name='sin'),
         dict(g=lambda t: mp.e ** (-(t - 58.5) ** 2) * 1e-3, j=lambda t: 1 - CB * 0.5 * mp.sin(t), name='gauss'),
         dict(g=lambda t: (t - 57) ** 3 * 1e-5 - 2e-4, j=lambda t: 1 + CB * 0.99 * mp.sin(0.3 * t), name='cubic')]
wsy = sp.symbols('w0:6'); zsy = sp.symbols('z0:5'); jsy, jpsy, csy = sp.symbols('jsy jpsy csy')
def _lam(c):
    rep = {}
    for d in c.atoms(sp.Derivative):
        n = d.derivative_count
        rep[d] = {Wf: wsy, Zf: zsy}[d.expr][n] if d.expr in (Wf, Zf) else jpsy
    c2 = c.xreplace(rep).xreplace({Wf: wsy[0], Zf: zsy[0], Cf: csy, Kf: 0, Jf: jsy})
    return sp.lambdify([*wsy, *zsy, jsy, jpsy, csy], c2, 'mpmath')
LF = {k: _lam(c) for k, c in L0.items()}
res = []
for T in TESTS:
    g = T['g']; j = T['j']
    gd = lambda t, i: mp.diff(g, t, i) if i else g(t)
    jp = lambda t: mp.diff(j, t)
    assert max(abs(jp(t)) for t in mp.linspace(57, 60, 31)) <= C1, 'test j\' exceeds C1_OR'
    Ce = lambda t: (t + STIP) / 4
    def X0(t):
        k4 = 4 * zd[1](t) * gd(t, 3) + 6 * zd[2](t) * gd(t, 2) + 4 * zd[3](t) * gd(t, 1) + zd[4](t) * g(t)
        k3 = 3 * zd[1](t) * gd(t, 2) + 3 * zd[2](t) * gd(t, 1) + zd[3](t) * g(t)
        return j(t) * k4 + jp(t) * k3 + Ce(t) * zd[1](t) * g(t)
    K0def = mp.quad(lambda t: W(t) * zeta(t) * g(t) * X0(t), [57, 58, 59, 60])
    nodes = mp.linspace(57, 60, 61)
    def red_integrand(t):
        Wt = mp.taylor(W, t, 5); Wd = [Wt[n] * mp.factorial(n) for n in range(6)]
        env = [*Wd, *[zd[i](t) for i in range(5)], j(t), jp(t), Ce(t)]
        return sum(LF[(a, b)](*env) * gd(t, a) * gd(t, b) for (a, b) in L0)
    K0red = mp.quad(red_integrand, [57, 58, 59, 60], maxdegree=6)
    Nl = [mp.sqrt(mp.quad(lambda t, l=l: W(t) * gd(t, l) ** 2, [57, 60])) for l in range(5)]
    bound0 = sum(kap[l] * Nl[l] ** 2 for l in range(3))
    def D2full(t):
        p1 = psi1(t); p2 = psi2(t)
        return (zeta(t) * gd(t, 4) + (4 * zd[1](t) - 2 * p1 * zeta(t)) * gd(t, 3) + (6 * zd[2](t) - 6 * p1 * zd[1](t) + p2 * zeta(t)) * gd(t, 2)
                + (4 * zd[3](t) - 6 * p1 * zd[2](t) + 2 * p2 * zd[1](t)) * gd(t, 1) + (zd[4](t) - 2 * p1 * zd[3](t) + p2 * zd[2](t)) * g(t))
    K2lin = mp.quad(lambda t: W(t) * D2full(t) * X0(t), [57, 58, 59, 60], maxdegree=6)
    jb = 1 + CB
    fD = aa[2] * Nl[2] + aa[1] * Nl[1] + aa[0] * Nl[0] + Nl[4] + aa[3] * Nl[3]
    fX = jb * 4 * z[1] * Nl[3] + jb * (6 * z[2] * Nl[2] + 4 * z[3] * Nl[1]) + X2['X0lin']['up'] * Nl[0] + C1 * (3 * z[1] * Nl[2] + 3 * z[2] * Nl[1] + z[3] * Nl[0])
    r = dict(test=T['name'], K0_defining=float(K0def), K0_reduced=float(K0red), rel_diff=float(abs(K0def - K0red) / max(abs(K0def), mp.mpf('1e-300'))),
             K0_bound=float(bound0), K0_le_bound=bool(K0def <= bound0), K2lin=float(K2lin), K2lin_bound=float(fD * fX), K2lin_le_bound=bool(abs(K2lin) <= fD * fX))
    res.append(r); print(r, flush=True)
# (2) reflected integral
random.seed(12345)
Cp, Sp, s = 1.0, 0.01, 0.05
n = 4000; dt = 0.01
phi_ = [0.0] * n
for _ in range(40):
    i0 = random.randrange(n); phi_[i0] += Cp * 0.9 / dt * random.random() / 40 * 10
def window_ok():
    pre = [0.0]
    for v in phi_: pre.append(pre[-1] + v * dt)
    for a in range(0, n, 7):
        for b in range(a + 1, n + 1, 13):
            if pre[b] - pre[a] > Cp + Sp * (b - a) * dt + 1e-12: return False
    return True
while not window_ok():
    phi_ = [v * 0.9 for v in phi_]
P = [0.0]; Zs = [0.0]; m = 0.0
for i in range(n):
    P.append(P[-1] + 2 * (phi_[i] - s) * dt); m = min(m, P[-1]); Zs.append(P[-1] - m)
zmax = max(Zs); slope_ok = all(Zs[i + 1] - Zs[i] >= 2 * (phi_[i] - s) * dt - 1e-12 for i in range(n))
refl = dict(Cp=Cp, Sp=Sp, s=s, Zmax=zmax, Z_le_2Cp=bool(zmax <= 2 * Cp + 1e-12), Zmin=min(Zs), slope_ok=bool(slope_ok))
print('reflected:', refl)
out = dict(stage='X6 flux2_sanity', label='float diagnostic, not part of the proof', tests=res, reflected=refl,
           all_ok=bool(all(r['rel_diff'] < 1e-15 and r['K0_le_bound'] and r['K2lin_le_bound'] for r in res) and refl['Z_le_2Cp'] and refl['slope_ok'] and refl['Zmin'] >= 0),
           args=FC.args_block([os.path.abspath(__file__), os.path.join(FC.CODE, 'f2common.py'), os.path.join(FC.OUT, 'flux2_coef.json'), os.path.join(FC.OUT, 'flux2_forms.pkl')]))
out['timing'] = dict(runtime_s=T_.s())
h = FC.dump(os.path.join(FC.OUT, 'flux2_sanity.json'), out)
print('X6 done all_ok', out['all_ok'], 'sha', h, '(%.1fs)' % T_.s())
