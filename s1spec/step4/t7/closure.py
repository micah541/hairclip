"""T7 float: closure rate rho_c of T5 (E0)+(E2)+mu_e[(I0)+(I2)] in linear mode, per weight and k; Xi, Ybar budget (T5 (Y), surplus
reserve kappa_D = 1/2 => Xi^2 doubled); J-cascade and T8 amplification with the per-weight constants (t8/amp.py, design opt4n).
Inputs: eicon.json (mu1, mu2, ramp c's, Gram of zeta g'), cx.json (modulation map), c1float_clamped_0.1.json (tip coercivity).
Assumed tip input: omega_eff = omega_c + m_c - s_c = 0.045, theta0 = 0.025 (float c(0.03) = 0.069-0.072, certificate loss ~0.005, s_c <= m_c/2).
usage: python3 closure.py  -> closure.json, prints the weight table"""
import sys, json, numpy as np
from scipy.optimize import minimize
sys.path.insert(0, '/home/claude/s1spec/step4/t8')
EPS, EPSI, C1, PHI1, PHI2 = 1e-5, 3e-5, 0.5804, 1.083, 1.3
G22, G12, G02, GM = 14.34, 16.3, 10.7, 43.33          # S1 Gaarding constants of the eps level (T7 C4 to redo along the base)
SA0 = 5.8e-7                                          # 2 IA K_rho(16.25) (strip A^2 + 2 A W1 part)
STRIP = {2: 1.68e-4, 1.5: 6.03e-6, 1: 2.24e-7}        # sA + sB (design formula, chi_32 strip), per P
D2STRIP = (6.841e-3, 1.5418e-2, 2.3648e-2)            # d3, d2, d1 at P = 2;  d_j(P) = d_j(2) 0.1753^(2(2-P))
OMEFF, TH0 = 0.045, 0.025

def lam_max(TT, G):
    return float(np.sqrt(max(np.linalg.eigvals(np.array(TT) @ np.array(G)).real)))

def rate(c, muE, nstart=40, seed=1):
    """min over shapes of Q(z)/Lambda(z) (float, multistart); c = dict of constants"""
    rng = np.random.default_rng(seed); P = c['P']
    dS = [d * 0.1753 ** (2 * (2 - P)) for d in D2STRIP]; c1W = c['c1W']
    def f(y):
        A, B, D4, V0, V2, V4 = np.exp(y)
        W1 = np.sqrt(A * B + C1 * A * A); W3 = np.sqrt(B * D4 + C1 * B * B)
        V1 = np.sqrt(V0 * V2 + c1W * V0 * V0); V3 = np.sqrt(V2 * V4 + c1W * V2 * V2)
        E0 = -OMEFF * A * A - TH0 * B * B + c['strip'] * V0 * V0 + SA0 * (A * A + 2 * A * W1)
        E2 = EPS * (-0.5 * D4 * D4 + G22 * B * B + G12 * W1 * B + G02 * A * B + GM * A * A + (dS[0] * V2 + dS[1] * V1 + dS[2] * V0) * (D4 + 2 * PHI1 * W3 + PHI2 * B))
        I0 = -c['mu1'] * V0 * V0 - 0.5 * V2 * V2 + c['a_f'] * W1 * W1 + c['b_f'] * B * B + c['Cg0'] * A * V0
        I2 = EPSI * (-0.5 * V4 * V4 - c['mu2'] * V2 * V2 + (c['c3'] * D4 + c['c2'] * W3 + c['c1'] * B + c['c0'] * W1 + c['Cg2'] * A) * (V4 + 2 * c['Psi1'] * V3 + 2 * c1W * V2))
        Ds = TH0 / 4 * B * B + EPS / 8 * D4 * D4 + muE / 8 * V2 * V2 + muE * EPSI / 8 * V4 * V4
        Q = -(E0 + E2 + muE * (I0 + I2)) - Ds
        L = A * A + EPS * B * B + muE * (V0 * V0 + EPSI * V2 * V2)
        return Q / L
    best = np.inf
    for _ in range(nstart):
        y0 = rng.normal(0, 3, 6)
        r = minimize(f, y0, method='L-BFGS-B', bounds=[(-60, 60)] * 6)
        best = min(best, r.fun)
    return best

def best_mu(c):
    lo, hi = np.log(1e-9), np.log(1.0); grid = np.linspace(lo, hi, 28)
    vals = [rate(c, np.exp(t), nstart=12) for t in grid]
    j = int(np.argmax(vals)); a, b = grid[max(j - 1, 0)], grid[min(j + 1, len(grid) - 1)]
    for _ in range(14):
        m1, m2 = a + (b - a) / 3, b - (b - a) / 3
        if rate(c, np.exp(m1), 12) < rate(c, np.exp(m2), 12): a = m1
        else: b = m2
    mu = np.exp(0.5 * (a + b)); return mu, rate(c, mu, 60)

def window(c, frac=0.5):
    """mu_e interval where rho_c >= frac * max"""
    mu, r = best_mu(c); out = [mu, mu]
    for s, i in ((1 / 1.25, 0), (1.25, 1)):
        m = mu
        while 1e-10 < m < 10 and rate(c, m * s, 12) >= frac * r: m *= s
        out[i] = m
    return mu, r, out

if __name__ == '__main__':
    tag = sys.argv[1] if len(sys.argv) > 1 else ''
    E = json.load(open(f'eicon{tag}.json')); CX = json.load(open('cx.json'))
    names = []
    for e in E:
        if e['weight'] not in names: names.append(e['weight'])
    res = {}
    for nm in names:
        rows = [e for e in E if e['weight'] == nm]; P = float(nm.split('|')[0][2:])
        per_k = []
        for e, cx in zip(rows, CX):
            c = dict(P=P, mu1=e['mu1'], mu2=e['mu2'], a_f=e['a_f'], b_f=e['b_f'], c0=e['c0'], c1=e['c1'], c2=e['c2'], c3=e['c3'],
                     Psi1=e['Psi1'], c1W=0.5 * e['Psi2'], strip=STRIP[P], Cg0=lam_max(cx['TTt'], e['Gam0']), Cg2=lam_max(cx['TTt'], e['Gam2']))
            mu, r = best_mu(c); c.update(k=e['k'], muE=mu, rho=r); per_k.append(c)
            print(f"{nm:12s} k={e['k']:8.1f}: mu1={c['mu1']:.4f} Cg0={c['Cg0']:.3e} Cg2={c['Cg2']:.3e}  best mu_e={mu:.3e}  rho_c={r:+.4f}", flush=True)
        # one mu_e for the whole run: maximise the minimum over k
        grid = np.exp(np.linspace(np.log(min(p['muE'] for p in per_k)) - 1.5, np.log(max(p['muE'] for p in per_k)) + 1.5, 25))
        mins = [min(rate(p, m, 12) for p in per_k) for m in grid]; j = int(np.argmax(mins)); muE = float(grid[j])
        rhos = [rate(p, muE, 60) for p in per_k]
        res[nm] = dict(P=P, per_k=per_k, muE=muE, rho_k=rhos)
        print(f"   {nm}: common mu_e = {muE:.3e}, rho_c(k) = " + ', '.join(f'{v:+.4f}' for v in rhos), flush=True)
    json.dump(res, open(f'closure{tag}.json', 'w'), indent=1)
