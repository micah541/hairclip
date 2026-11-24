"""T7 float: T8's Ybar requirement per E_I weight, with the J-cascade fed by the per-weight T7 constants (not 'J fixed').
For each weight: mu_e window where rho_c >= 0.95 max (closure.rate, min over the 5 k), then t8/amp.py (design opt4n, all 188 slabs)
with  k_jm -> k_jm W_2/W  at the ramp sup (x = kat + X0),  jcas params = (mu_e, a_f, b_f, mu1, mu2, c_j, c1W, Psi1, Cg0, Cg2, rho_c, th0).
Ybar budget (T5 (Y), Xi^2 doubled for kappa_D = 1/2) from T1's residual table; Y(tau0) <= 2e-12 (T4 budget).  -> t8req.json"""
import sys, json, numpy as np, warnings
warnings.filterwarnings('ignore')
sys.path.insert(0, '/home/claude/s1spec/step4/t8')
import amp, jcas
from fb import X0
import closure as CL
CXF = 7.5                     # |x| <= C_x Ybar for Lemma F (normalised directions; cx.json: 6.5-7.5); T8 used 4
tag = sys.argv[1] if len(sys.argv) > 1 else ''
R = json.load(open(f'closure{tag}.json'))
per, tau1 = amp.load('opt4n')
# T1 residuals (k >= 24 clean slabs): R2, Rs (rho zone, W-independent); R1, RIs for W = (1+x^2)^-2; RF0, RF2 give the slow-weight proxies
T1 = dict(R2=5.4e-16, Rs=4.2e-14, R1_2=6.3e-15, RIs_2=2.9e-18, RF0=9.5e-12, RF2=1.5e-10)
def resid(P):
    if P == 2: return T1['R1_2'], T1['RIs_2']
    # (1+x^2)^-P <= (1+x^2)^-1 and (1+x^2)^-1 <= x^-2 on y >= 14; [10.25,16.25] part bounded by the P=2 value x (1+16.25^2)^(2-P)
    R1 = T1['RF0'] / np.sqrt(28) + T1['R1_2'] * (1 + 16.25**2) ** ((2 - P) / 2)
    RIs = T1['RF2'] / (np.sqrt(28) * 16.25**2) + T1['RIs_2'] * (1 + 16.25**2) ** ((2 - P) / 2)
    return R1, RIs
out = {}
ONLY = sys.argv[2].split(';') if len(sys.argv) > 2 else None
for nm, r in R.items():
    if ONLY and nm not in ONLY: continue
    P = r['P']; pk = r['per_k']
    grid = np.exp(np.linspace(np.log(1e-6), np.log(0.5), 31))
    rho_min = np.array([min(CL.rate(p, m, 10) for p in pk) for m in grid])
    rmax = rho_min.max(); ok = grid[rho_min >= 0.95 * rmax]
    if rmax <= 0 or len(ok) == 0: print(f' => {nm}: no closing mu_e'); continue
    worst = dict(mu1=min(p['mu1'] for p in pk), mu2=min(p['mu2'] for p in pk), af=max(p['a_f'] for p in pk), bf=max(p['b_f'] for p in pk),
                 cI=tuple(max(p[f'c{j}'] for p in pk) for j in range(4)), c1W=max(p['c1W'] for p in pk), Psi1=max(p['Psi1'] for p in pk),
                 Cg0=max(p['Cg0'] for p in pk), Cg2=max(p['Cg2'] for p in pk))
    Ws = (lambda P_: (lambda j, m, kat: (1.0 + (kat + X0)**2) ** (P_ - 2)))(P)
    rows = []
    for mu in np.exp(np.linspace(np.log(ok.min()), np.log(ok.max()), 7)):
        rho = min(CL.rate(p, mu, 30) for p in pk)
        prm = dict(muE=mu, th0=CL.TH0, rho=rho, af=worst['af'], bf=worst['bf'], mu1=worst['mu1'], mu2=worst['mu2'], cI=worst['cI'],
                   c1W=worst['c1W'], Psi1=worst['Psi1'], gam0=worst['Cg0'] / CXF, gam2=worst['Cg2'] / CXF, kD=0.5)
        Rg = amp.integrate(per, tau1, P=prm, Cx=CXF, Wscale=None if P == 2 else Ws)
        h = amp.handoff(Rg, per, tau1); L = h['lim']
        # Ybar budget
        R1, RIs = resid(P); R2, Rs = T1['R2'], T1['Rs']
        Xi = np.sqrt(2 * ((1 + CL.PHI1)**2 * R2**2 / CL.TH0 + 2 * CL.EPS * Rs**2 + 2 * mu * R1**2 + 2 * mu * CL.EPSI * RIs**2))
        ac = (CL.PHI2 + CL.PHI1 * (2 * np.sqrt(CL.C1) + 1) + 5.4 * 3.0) * R2 + np.sqrt(mu) * (1.0 + worst['Cg0']) * R1
        Yres = ac / rho + Xi / np.sqrt(rho)
        req = min(L['item3_E3'], L['item4_K4_goodtime'], min(L['item2_E0'], L['item2_E1'], L['item2_E2']))
        J = Rg['diag'][3]['J'][-1]
        rows.append(dict(muE=float(mu), rho=float(rho), A=h['A'], item3=L['item3_E3'], K4good=L['item4_K4_goodtime'], K4F4=L['item4_K4_viaF4'],
                         req=float(req), Xi=float(Xi), Yres=float(Yres), J=[float(J[m]) for m in range(5)], R1=float(R1), RIs=float(RIs)))
        print(f"{nm:12s} mu_e={mu:.2e} rho_c={rho:.4f}  J0..4=" + ','.join(f'{J[m]:.2e}' for m in range(5)) +
              f"  A3={h['A'][3]:.2e} A4={h['A'][4]:.2e}  Ybar req: item3 {L['item3_E3']:.2e} K4good {L['item4_K4_goodtime']:.2e} -> {req:.2e}"
              f" | Xi={Xi:.1e} Yres={Yres:.1e}", flush=True)
    b = max(rows, key=lambda q: q['req'])
    out[nm] = dict(P=P, window=[float(ok.min()), float(ok.max())], rho_max=float(rmax), worst=worst, rows=rows, best=b)
    print(f" => {nm}: mu_e window [{ok.min():.1e}, {ok.max():.1e}], best req {b['req']:.2e} at mu_e {b['muE']:.1e} (rho_c {b['rho']:.4f})", flush=True)
json.dump(out, open(f't8req{tag}.json', 'w'), indent=1)
