"""F3b: tail count of sym_rho(L_Gamma + transport) on the macroscopic-phase curves.
usage: python3 f3b.py [--h 0.05] [--sigma 35] statefiles...   (state_t*.npz) ; frame: tip, lam = 1.5398/kappa_tip.
Variants: A = frame terms (q/4, 1/4) + transport tau (the F3 operator);  B = A without tau;  C = frame-free physical operator (no q/4, no 1/4, tau from -kappa'' only)."""
import sys, json, argparse, numpy as np, scipy.linalg as sla
sys.dont_write_bytecode = True
sys.path.insert(0, '/home/claude/s1spec/step4/feas2')
from f3geom import *; from f3fem import *
ap = argparse.ArgumentParser(); ap.add_argument('--h', type=float, default=0.05); ap.add_argument('--sigma', type=float, default=35.0)
ap.add_argument('--out', default='logs/f3b.json'); ap.add_argument('files', nargs='+'); a = ap.parse_args()
CS = (0.05, 0.1, 0.25, 0.5, 1.0)

def coeffsV(mesh, sp, var, g=None):
    x = mesh.xq; k, k1, k2 = sp['kap'](x), sp['kap1'](x), sp['kap2'](x)
    s = sp['var'][var]; mu = s['mu']
    return dict(b2=-k**2, b1=-5*k*k1 - mu*sp['q'](x) + s['wT'](x) + (s['tau'](x) if s['use_tau'] else 0.0), b0=mu - 3*k1**2 - 4*k*k2)
def build_var(g, mu, wc, use_tau):
    """frame data: dilation rate mu, frame-centre velocity wc (complex, lambda^3 * dxc/dt):  V = -kappa'' + mu p - wc.N ;  tau = int kappa V"""
    N_ = 1j * g['T']; wN = (wc * np.conj(N_)).real; wT = (wc * np.conj(g['T'])).real
    V = -g['kap2'] + mu * g['p'] - wN
    return dict(mu=mu, use_tau=use_tau, tau=CubicSpline(g['x'], g['kap'] * V).antiderivative(), wT=CubicSpline(g['x'], wT))
from scipy.interpolate import CubicSpline
def frame_rates(f, T=4.5988639266):
    """mu = lam^4 dln(kappa_tip)/dt,  wc = lam^3 d(x_c)/dt  (x_c = tip + D_S1 lam), by a short forward run"""
    sys.path.insert(0, '/home/claude/s1spec/step4/feas2')
    from run2 import load_state
    import copy
    S = load_state(f); t0 = S.t; dt = min(1e-3, 0.005 * (T - t0)); 
    def obs(S):
        al, th, Z, _ = S.curve(1); N = len(al); kt = abs((2 * np.pi / S.L) * (S.om + S.d(S.th_hat, 1))[N // 2]); lam = K0S1 / kt
        return kt, Z[N // 2] + H.D_S1 * lam, lam
    k0, x0, l0 = obs(S); S.run_to(tend=t0 + dt); k1, x1, l1 = obs(S)
    mu = l0**4 * np.log(k1 / k0) / dt; wc = l0**3 * (x1 - x0) / dt
    return float(mu), complex(wc), float(l0)

def minM(U, S, c):
    SU = U.T @ S @ U; n = SU.shape[0]
    top = lambda M: np.linalg.eigvalsh(SU[M:, M:])[-1]
    lo, hi = 0, min(n - 1, 40)
    if top(0) <= -c: return 0
    if top(hi) > -c: return -1
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if top(mid) <= -c: hi = mid
        else: lo = mid
    return hi

mesh = Mesh(a.sigma, a.h); M = mesh.mass()
spS1 = s1_splines(Smax=a.sigma + 5, dx=0.01)
K1 = mesh.form(coeffsV(mesh, dict(spS1, var={'B': dict(mu=0.25, use_tau=False, tau=None, wT=lambda x: 0*x)}), 'B')); wS1, U1 = eig_sym(K1, M)
res = {'S1': dict(top=[float(x) for x in wS1[:8]], n={str(c): int(np.sum(wS1 > -c)) for c in CS})}
print('S1 top', np.round(wS1[:7], 4), 'n', res['S1']['n'], flush=True)
rho_q = mesh.rho * mesh.wq
for f in a.files:
    z = np.load(f); t = float(z['t'])
    g = numeric_geometry(f, Smax=a.sigma + 5, dsig=0.01)
    # frame-free transport: tauC = int kappa (-kappa'')
    from scipy.interpolate import CubicSpline
    sp = as_splines(g)
    mu_t, wc_t, _ = frame_rates(f)
    sp['var'] = {'A': build_var(g, 0.25, 0j, True), 'B': build_var(g, 0.25, 0j, False), 'E': build_var(g, mu_t, wc_t, True), 'C': build_var(g, 0.0, 0j, True)}
    Y = g['Y']; sig = g['x']
    ax = np.nonzero(np.sign(Y.imag[1:-1]) * np.sign(Y.imag[2:]) < 0)[0] + 1
    cr = [round(float(sig[i] + (sig[i+1]-sig[i]) * Y.imag[i] / (Y.imag[i] - Y.imag[i+1])), 2) for i in ax]
    out = dict(t=t, lam=float(g['lam']), kt=float(g['kt']), axis_crossing_sigma=cr[:12], mu=mu_t, wc=[wc_t.real, wc_t.imag])
    rhoint = lambda lo, hi: float(np.sum(np.where((mesh.xq >= lo) & (mesh.xq < hi), rho_q, 0)))
    tot = rhoint(0, 1e9); out['rho_mass_frac'] = {f'{lo}-{hi}': rhoint(lo, hi) / tot for lo, hi in ((0, 2), (2, 5), (5, 10), (10, 20), (20, 40))}
    for var in 'ABEC':
        K = mesh.form(coeffsV(mesh, sp, var)); S = 0.5 * (K + K.T)
        w, V = eig_sym(K, M)
        o = dict(top=[round(float(x), 4) for x in w[:10]], n={str(c): int(np.sum(w > -c)) for c in CS}, M={str(c): minM(U1, S, c) for c in CS})
        if var == 'A':   # localisation of the modes above -0.25 : median / 90% sigma of rho u^2, and S1-mode overlaps
            xs_ = np.r_[0.0, np.arange(1, 0) ]  # placeholder
            loc = []
            xn = np.linspace(0, a.sigma, mesh.ne + 1)       # node positions (values dofs at even index after dropping slope dof 1)
            for m in range(min(8, V.shape[1])):
                v = V[:, m]; val = np.zeros(mesh.ndof); val[mesh.free] = v
                u = val[0::2]; r = np.exp(-0.25 * (1 + xn**2) ** (2 / 3)); e = r * u**2; cs = np.cumsum(e) / e.sum()
                loc.append(dict(eig=round(float(w[m]), 4), s50=float(xn[np.searchsorted(cs, .5)]), s90=float(xn[np.searchsorted(cs, .9)]),
                                S1overlap=[round(float(abs(U1[:, i] @ M @ v) / np.sqrt((v @ M @ v) * (U1[:, i] @ M @ U1[:, i]))), 3) for i in range(4)]))
            o['modes'] = loc
        out[var] = o
        print(f"t={t:.4f} lam={g['lam']:.4f} {var}: top {o['top'][:8]} n(.05/.1/.25/.5/1)={[o['n'][str(c)] for c in CS]} M={[o['M'][str(c)] for c in CS]}", flush=True)
    print('   axis crossings at sigma =', cr[:8], f'mu={mu_t:.4f} wc={wc_t:.3f}', ' rho mass:', {k: round(v, 3) for k, v in out['rho_mass_frac'].items()}, flush=True)
    for m in out['A']['modes'][:6]: print('   mode', m)
    res[f't={t:g}'] = out
json.dump(res, open(a.out, 'w'), indent=1)
