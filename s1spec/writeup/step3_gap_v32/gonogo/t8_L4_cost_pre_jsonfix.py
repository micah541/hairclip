"""CLEAN-GO task 3 (FLOAT estimates): cost of the general-weight L^4 terms (derivation: L4_general_weight.md; identities: l4_identity_check.py)
relative to the power-weight tables.  Weight of Lemma U:  Theta = W(sigma) D(s),  sigma = s - s_rho(tau),  s = ln y,
    W = 1/(e^{-mu(sigma+a)} + 1 + e^{nu(sigma - ln1.5 - a)}),   D = 1 + g e^{-4s},   theta(s) = g e^{-4s}/D,   theta_a = theta(ln 28) = g/(28^4 + g),
on the geometry s >= ln 28 (foot of the inner ramp [28,56]) and s_rho >= ln y_L (crude family: y_L = 56; fine family [150,250]: y_L = 150).
Constants (all sups over that geometry; grid + the exact sigma -> +-inf limits):
    Gamma_m = sup|Theta^(m)/Theta|, m = 1..3 (the F1 tables use weight-derivative orders m <= 3);  Gamma_eff = max_m Gamma_m^(1/m)  [tables: |coef| Gamma_eff^m >= |coef| Gamma_m]
    kappa   = sup|L|,  L = (ln eps Theta)' = -4 + W'/W - 4 theta                       [L^4 class bounds, (I1), (I2)]
    kappa2  = sup(L^2 + L')_+ = sup((eps Theta)''/(eps Theta))_+                       [interpolation step (I5); the power tables use kappa^2]
Rates: Ghat_jj from F.constants(Gamma_eff) (core YC) + l4_extra_YC(kap = kappa, s33 = 'ibp') + ghat_from_YC(kap = sqrt(kappa2)), as t2_rates.py / lemmaS_Bbox.py.
Damping: g >= Ghat_22/(1 - Ghat_22/28^4)  <=>  theta_a = Ghat_22/28^4 (level 2 binds; used for all levels: conservative) -> fixed point theta_a <-> Ghat_22.
References: power weight of Lemma S(A) (gamma = 1, kappa = 5); the previous agent's general-weight shortcut (Gamma = 4, kappa = nu + 4, no damping-factor derivatives).
Inner-ramp L^4 foot terms (est.): (A)'s certified foot weights (rigor/weights_exact: footw_2 = 2862.49, footw_1 = 1.5707) scaled by the ratio of the L^4 extras
(general / power), times the space-time cutoff-free ramp data F_j = eps(28) Gw IW 2 E^A_j (IW = 4 int W dsigma); the sup-free S33 term likewise (formula of lemmaS.tex l.331).
Run from anywhere (reads writeup/lemmaS read-only; no bytecode written).  Output: t8_L4_cost.json (read by t1_crude_shell_model.py)."""
import os, sys, json, math
import numpy as np, sympy as sp
HERE = os.path.dirname(os.path.abspath(__file__)); LS = '/home/claude/s1spec/writeup/lemmaS'
sys.dont_write_bytecode = True; sys.path.insert(0, LS); os.chdir(LS)
import lemmaS_F1 as F, lemmaS_L4terms as L4, run_F1 as R
os.chdir(HERE)
LN28 = math.log(28.)
xs = sp.Symbol('x')

def wfuncs(mu, nu, a):
    W = 1/(sp.exp(-mu*(xs + a)) + 1 + sp.exp(nu*(xs - sp.log(sp.Rational(3, 2)) - a)))
    return [sp.lambdify(xs, sp.diff(W, xs, k)/W, 'numpy') for k in (1, 2, 3)], sp.lambdify(xs, W, 'numpy')
_WF = {}
def wfuncs_c(mu, nu, a):
    if (mu, nu, a) not in _WF: _WF[(mu, nu, a)] = wfuncs(mu, nu, a)
    return _WF[(mu, nu, a)]

def weight_consts(mu, nu, a, theta_a, yL, ns=401, nsig=1601):
    g = theta_a*28.**4/(1 - theta_a)
    wf, _ = wfuncs_c(mu, nu, a)
    s = np.linspace(LN28, LN28 + 16, ns)[None, :]; sig = np.linspace(-16, 16, nsig)[:, None]
    th = g*np.exp(-4*s)/(1 + g*np.exp(-4*s))
    mask = sig <= s - math.log(yL)                       # s_rho = s - sigma >= ln yL
    w = [np.ones_like(sig)] + [f(sig)*np.ones_like(sig) for f in wf]
    d = [np.ones_like(th)] + [(-4.)**k*th for k in (1, 2, 3)]
    G = []
    for m in (1, 2, 3):
        Tm = sum(math.comb(m, k)*w[k]*d[m - k] for k in range(m + 1))
        G.append(float(np.max(np.abs(np.where(mask, Tm, 0.)))))
    L = -4 + w[1] - 4*th; Lp = (w[2] - w[1]**2) + 16*th*(1 - th)
    kap = float(np.max(np.abs(np.where(mask, L, 0.)))); kap2 = float(np.max(np.where(mask, L*L + Lp, 0.)))
    # exact limits sigma -> -inf (w_k = mu^k; any s) and sigma -> +inf (w_k = (-nu)^k; then s -> inf, theta -> 0)
    for wk, thv in (([mu**k for k in (1, 2, 3)], th), ([(-nu)**k for k in (1, 2, 3)], np.zeros_like(th))):
        ww = [1.] + wk; dd = [np.ones_like(thv)] + [(-4.)**k*thv for k in (1, 2, 3)]
        for i, m in enumerate((1, 2, 3)): G[i] = max(G[i], float(np.max(np.abs(sum(math.comb(m, k)*ww[k]*dd[m - k] for k in range(m + 1))))))
        Ll = -4 + ww[1] - 4*thv; kap = max(kap, float(np.max(np.abs(Ll)))); kap2 = max(kap2, float(np.max(Ll*Ll + 16*thv*(1 - thv))))
    return dict(g=g, theta_a=theta_a, Gamma_m=G, Gamma_eff=max(G[m]**(1/(m + 1)) for m in range(3)), kappa=kap, kappa2=kap2)

T = F.build_tables()
for j in T:
    for mons in T[j].values(): assert all(float(co) >= 0 for co in mons.values())     # tables store |coef|: gam^m -> Gamma_eff^m is an upper bound
_AB = {}
def core_YC(gam, M0, M1):
    """the core part of F.constants (lines 86-95 of lemmaS_F1.py), without the (unused, slow) ramp tables"""
    if (M0, M1) not in _AB: _AB[(M0, M1)] = F.A_bound_set('quasi', R.box(0, 0, 0, M0=M0, M1=M1))
    Ab = _AB[(M0, M1)]; G = {}
    for j in (0, 1, 2):
        n = j + 2; YC = {k: 0.0 for k in range(n + 1)}
        for (a_, b_), mons in T[j].items():
            if (a_, b_) == (n, n): continue
            gC = sum(co*gam**m*Ab[(i, r)] for (m, i, r), co in mons.items())
            if a_ == b_: YC[a_] += gC
            else: YC[a_] += gC/2; YC[b_] += gC/2
        G[j] = {k: float(v) for k, v in YC.items()}
    return G
_cc, _G = F.constants(4.0, R.box(0, 0, 0, M0=0.11, M1=0.52), 'quasi', tables=T)
_chk = core_YC(4.0, 0.11, 0.52)
assert all(abs(float(_G[j]['YC'][k]) - _chk[j][int(k)]) <= 1e-12*abs(_chk[j][int(k)]) for j in range(3) for k in _G[j]['YC']), 'core_YC != F.constants'
print('core_YC reproduces F.constants YC (gam = 4, crude box)', flush=True)
_GH = {}
def ghat(gam, kap, kap2, M0=0.11, M1=0.52):
    key = (round(gam, 9), round(kap, 9), round(kap2, 7), M0, M1)
    if key not in _GH: _GH[key] = _ghat(gam, kap, kap2, M0, M1)
    return _GH[key]
def _ghat(gam, kap, kap2, M0=0.11, M1=0.52):
    amin = (1 + M0**2)**-2; G = core_YC(gam, M0, M1); raw, l4, ok = [], [], True
    for j in range(3):
        YC0 = dict(G[j])
        raw.append(L4.ghat_from_YC(YC0, j, amin, kap=math.sqrt(kap2))[j])
        add, log, dis, spare = L4.l4_extra_YC(T, j, M0, M1, gam, kap=kap, s33='ibp'); ok = ok and dis <= spare*(1 + 1e-9)
        YC = {k: YC0[k] + add.get(k, 0) for k in YC0}; l4.append(L4.ghat_from_YC(YC, j, amin, kap=math.sqrt(kap2))[j])
    return raw, l4, bool(ok)

out = dict(refs={}, cases={})
for name, gam, kap in (('power weight, Lemma S(A) core (gamma=1, kappa=5)', 1.0, 5.0), ('previous T2 v1 (Gamma=4, kappa=5)', 4.0, 5.0),
                       ('previous shortcut (Gamma=4, kappa=nu+4=8)', 4.0, 8.0), ('previous shortcut (Gamma=4, kappa=nu+4=6)', 4.0, 6.0)):
    raw, l4, ok = ghat(gam, kap, kap*kap); out['refs'][name] = dict(Ghat_no_L4=raw, Ghat_L4=l4, L4_factor_22=l4[2]/raw[2], budget_ok=ok)
    print(f'{name}: Ghat_jj no L4 {["%.4g" % x for x in raw]}, with L4 {["%.4g" % x for x in l4]}, L4 factor (22) x{l4[2]/raw[2]:.2f}', flush=True)
for (mu, nu, a) in ((4, 4, 0.3), (4, 2, 0.5), (4, 4, 0.0)):
    for fam, yL, M0, M1 in (('crude', 56., 0.11, 0.52), ('fine', 150., 0.21, 0.21)):
        th = 0.39; hist = []
        for it in range(12):
            wc = weight_consts(mu, nu, a, th, yL)
            raw, l4, ok = ghat(wc['Gamma_eff'], wc['kappa'], wc['kappa2'], M0, M1)
            new = l4[2]/28.**4; hist.append((th, l4[2]))
            if new >= 0.999: th = None; break
            if abs(new - th) < 1e-4: th = new; break
            th = new
        key = f'{fam}|mu={mu},nu={nu},a={a}'
        if th is None:
            out['cases'][key] = dict(closed=False, history=hist); print(key, 'NO FIXED POINT (Ghat_22 >= 28^4)', hist, flush=True); continue
        wc = weight_consts(mu, nu, a, th, yL); raw, l4, ok = ghat(wc['Gamma_eff'], wc['kappa'], wc['kappa2'], M0, M1)
        rawp, l4p, _ = ghat(1.0, 5.0, 25.0, M0, M1)                      # power reference on the same box
        extra_ratio = (l4[2] - raw[2])/(l4p[2] - rawp[2])
        out['cases'][key] = dict(closed=True, weight=wc, Ghat_no_L4=raw, Ghat_L4=l4, budget_ok=ok, L4_factor_22=l4[2]/raw[2],
                                 power_ref=dict(Ghat_no_L4=rawp, Ghat_L4=l4p, L4_factor_22=l4p[2]/rawp[2]), L4_extras_ratio_general_over_power=extra_ratio,
                                 rate_margin_28=28.**4/l4[2], rate_margin_40=40.**4/l4[2], Gw=1/(1 - l4[2]/28.**4), history=hist)
        json.dump(out, open(os.path.join(HERE, 't8_L4_cost.json'), 'w'), indent=1)
        print(f"{key}: theta_a = {th:.4f}, g = {wc['g']:.4g}; Gamma_m = {[round(x, 3) for x in wc['Gamma_m']]}, Gamma_eff = {wc['Gamma_eff']:.4f}, "
              f"kappa = {wc['kappa']:.4f} (nu+4 = {nu + 4}), kappa2 = {wc['kappa2']:.3f} (kappa^2 = {wc['kappa']**2:.2f})", flush=True)
        print(f"   Ghat_jj no L4 {['%.4g' % x for x in raw]}, with general-weight L4 {['%.4g' % x for x in l4]} (budget ok {ok}); L4 factor (22) x{l4[2]/raw[2]:.2f} "
              f"[power weight, same box: x{l4p[2]/rawp[2]:.2f}]; L4 extras general/power x{extra_ratio:.2f}; rate margin 28^4/Ghat_22 = {28.**4/l4[2]:.3f}, Gw = {1/(1 - l4[2]/28.**4):.3f}", flush=True)
# ---- inner-ramp L^4 foot terms (est.), crude family, weight (4,4,0.3): data F_j from E^A (certified rerun), footw from weights_exact (certified, power weight)
EA = json.load(open(os.path.join(HERE, 'lemmaS_A_rerun/rigor/assembly_arb.json')))['E']
c = out['cases'].get('crude|mu=4,nu=4,a=0.3')
if c and c['closed']:
    _, Wl = wfuncs(4, 4, 0.3); xx = np.linspace(-15, 15, 300001); IW = 4*float(np.sum(Wl(xx))*(xx[1] - xx[0])); Gw = c['Gw']; eps28 = 28.**-4
    Fj = [eps28*Gw*IW*2*e for e in EA]                                   # space-time cutoff-free ramp data int int eps Theta z_j^2 on [28,56]
    footw = {1: 1.5707250797145307, 2: 2862.4867820255}; r = c['L4_extras_ratio_general_over_power']
    foot = {j: footw[j]*max(r, 1.)*Fj[2] for j in (1, 2)}               # (A): foot_j = footw_j * Jfree2 (level-2 data in both levels)
    # S33 sup-free (lemmaS.tex l.331): K_f sqrt(F_max) (J_F/sqrt(l) + sqrt2 J_F^{3/4} J_G^{1/4}); (A): K_f = 1.17e4 -> scaled by r; F_max = eps28 Gw 2 E2^A; J_G from (A)'s level-1 dissipation
    amin = (16/17)**2; gp1 = 0.24617403; JG = 2*(gp1*EA[1]/(2*amin))*IW*Gw
    s33 = 1.17e4*max(r, 1.)*math.sqrt(eps28*Gw*2*EA[2])*(Fj[2]/math.sqrt(math.log(2)) + math.sqrt(2)*Fj[2]**0.75*JG**0.25)
    out['ramp_foot_est'] = dict(IW=IW, Gw=Gw, F=Fj, foot=foot, s33=s33, total_level2=foot[2] + s33)
    print(f"inner-ramp L4 foot terms (est., crude, mu=nu=4, a=0.3): F_j = {['%.3g' % x for x in Fj]}; foot_1 = {foot[1]:.3g}, foot_2 = {foot[2]:.3g}, S33 = {s33:.3g}; level-2 total {foot[2] + s33:.3g}")
json.dump(out, open(os.path.join(HERE, 't8_L4_cost.json'), 'w'), indent=1)
