"""T7 float sanity check (not certified): T5 closure rate rho_c with the certified T7 constants (worst over each certificate interval),
using closure.py's rate() with (E2) in AB form (G02 + G_x) A B instead of G_m A^2, omega_eff = c_I, theta_0 = 0.025.
C_gamma0/2 and G_x use the per-component bound with C_x,i at tau_c (certificate output) when the c2r3 batch is pending."""
import os, sys, json, math
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
from scipy.optimize import minimize
import closure as CL

EPS, EPSI, C1, PHI1, PHI2 = 1e-5, 3e-5, 0.5804, 1.083, 1.3


def rate(c, muE, nstart=30, seed=1):
    rng = np.random.default_rng(seed)
    def f(y):
        A, B, D4, V0, V2, V4 = np.exp(y)
        W1 = np.sqrt(A * B + C1 * A * A); W3 = np.sqrt(B * D4 + C1 * B * B)
        V1 = np.sqrt(V0 * V2 + c['c1W'] * V0 * V0); V3 = np.sqrt(V2 * V4 + c['c1W'] * V2 * V2)
        E0 = -c['om'] * A * A - c['th0'] * B * B + c['strip'] * V0 * V0 + c['sA0'] * (A * A + 2 * A * W1)
        E2 = EPS * (-0.5 * D4 * D4 + c['G22'] * B * B + c['G12'] * W1 * B + (c['G02'] + c['Gx']) * A * B
                    + (c['d3'] * V2 + c['d2'] * V1 + c['d1'] * V0) * (D4 + 2 * PHI1 * W3 + PHI2 * B))
        I0 = -c['mu1'] * V0 * V0 - 0.5 * V2 * V2 + c['a_f'] * W1 * W1 + c['b_f'] * B * B + c['Cg0'] * A * V0
        I2 = EPSI * (-0.5 * V4 * V4 - c['mu2'] * V2 * V2 + (c['c3'] * D4 + c['c2'] * W3 + c['c1'] * B + c['c0'] * W1 + c['Cg2'] * A)
                     * (V4 + 2 * c['Psi1'] * V3 + 2 * c['c1W'] * V2))
        Ds = c['th0'] / 4 * B * B + EPS / 8 * D4 * D4 + muE / 8 * V2 * V2 + muE * EPSI / 8 * V4 * V4
        Q = -(E0 + E2 + muE * (I0 + I2)) - Ds
        L = A * A + EPS * B * B + muE * (V0 * V0 + EPSI * V2 * V2)
        return Q / L
    best = np.inf
    for _ in range(nstart):
        r = minimize(f, rng.normal(0, 3, 6), method='L-BFGS-B', bounds=[(-60, 60)] * 6); best = min(best, r.fun)
    return best


def main():
    C = json.load(open(os.path.join(HERE, 'constants.json'))); W = C['meta']['weight_consts']; out = []
    for cert in C['certificates']:
        lo, hi = cert['slabs']; S = [C['slabs'][str(s)] for s in range(lo, hi + 1)]
        if not all('E_I' in e and 'E2' in e for e in S): continue
        mx = lambda f: max(f(e) for e in S); mn = lambda f: min(f(e) for e in S)
        Cxi = cert['C3_tau_c'] and [cert['C3_tau_c'].get(f'Cx{i}') for i in (1, 2, 3)]
        def mod(key):
            vals = [e['modulation'][key] if 'modulation' in e else float('nan') for e in S]
            return max(vals) if vals else float('nan')
        c = dict(om=cert['c_I'], th0=cert['theta'], strip=mx(lambda e: e['strip']['sA'] + e['strip']['sB']), sA0=mx(lambda e: e['strip']['sA0']),
                 G22=mx(lambda e: e['E2']['G22']), G12=mx(lambda e: e['E2']['G12']), G02=mx(lambda e: e['E2']['G02']),
                 d3=mx(lambda e: e['strip']['d3']), d2=mx(lambda e: e['strip']['d2']), d1=mx(lambda e: e['strip']['d1']),
                 mu1=mn(lambda e: e['E_I']['mu1']), mu2=mn(lambda e: e['E_I']['mu2']), a_f=mx(lambda e: e['ramp']['a_f']), b_f=mx(lambda e: e['ramp']['b_f']),
                 c3=mx(lambda e: e['ramp']['c3']), c2=mx(lambda e: e['ramp']['c2']), c1=mx(lambda e: e['ramp']['c1']), c0=mx(lambda e: e['ramp']['c0']),
                 Psi1=W['Psi1'], c1W=W['c1W'], Cg0=mod('Cgamma0'), Cg2=mod('Cgamma2'), Gx=mod('Gx'))
        if any(math.isnan(v) for v in c.values()): continue
        grid = np.exp(np.linspace(np.log(1e-3), np.log(0.2), 16)); vals = [rate(c, m, 12) for m in grid]
        j = int(np.argmax(vals)); r = rate(c, grid[j], 40)
        out.append(dict(id=cert['id'], k_c=cert['k_c'], rho_c=r, mu_e=float(grid[j]), consts=c))
        print(f"cert {cert['id']:2d} k_c {cert['k_c']:8.1f}: rho_c ~ {r:+.4f} at mu_e {grid[j]:.3g}  (om {c['om']:.4f}, Cg0 {c['Cg0']:.3g}, Cg2 {c['Cg2']:.3g}, Gx {c['Gx']:.3g})", flush=True)
    json.dump(out, open(os.path.join(HERE, 'closure7.json'), 'w'), indent=1)


if __name__ == '__main__':
    main()
