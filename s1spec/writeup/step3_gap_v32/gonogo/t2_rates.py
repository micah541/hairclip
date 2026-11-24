"""GO/NO-GO test T2 (+T4, rough): rate constants Ghat_jj for a general weight with |W^(m)| <= Gamma^m W.
The tables of lemmaS_prop1_bounds.table already store ABSOLUTE values of the coefficients of om^(m) A_i^(r), and lemmaS_L4terms uses |coef| gam^m
as well, so replacing gam by Gamma is the general-weight bound (no further loss of cancellations).  L^4 'extras' as in lemmaS_Bbox.py (s33 = 'ibp'):
the general-weight L^4 identities are those of L4_general_weight.md: the power-weight ones with |Theta'| <= kappa Theta, kappa = sup|Theta'/Theta|
(v1 of this test left kappa at its default 5; corrected here).
Criterion: damping from the factor (1 + g y^-4) needs g/(1 + g y_a^-4) >= Ghat, i.e. Ghat < y_a^4 at the inner end y_a of the support (28^4 = 6.15e5; 40^4 = 2.56e6).
Also prints Gamma_m for the two-sided weight W = 1/(e^{-mu s} + 1 + e^{nu (s - ln 1.5)}).   FLOAT.  run from writeup/lemmaS_v32."""
import sys, json, math
sys.path.insert(0, '.')
import numpy as np, sympy as sp
import lemmaS_F1 as F, lemmaS_L4terms as L4, run_F1 as R
out = {}
s = sp.Symbol('s')
for mu, nu in ((4, 2), (6, 4), (4, 4)):
    W = 1/(sp.exp(-mu*s) + 1 + sp.exp(nu*(s - sp.log(sp.Rational(3, 2)))))
    xs = np.linspace(-12, 12, 48001); G = []
    for m in range(1, 5):
        f = sp.lambdify(s, sp.diff(W, s, m)/W, 'numpy'); G.append(float(np.max(np.abs(f(xs)))))
    geff = max(g**(1/(m + 1)) for m, g in enumerate(G)); out[f'weight_mu{mu}_nu{nu}'] = dict(Gamma_m=G, Gamma_eff=geff)
    print(f'weight mu={mu} nu={nu}: Gamma_m = {[round(g, 2) for g in G]}, Gamma_eff = max Gamma_m^(1/m) = {geff:.3f}; with the damping factor add <= 4', flush=True)
T = F.build_tables()
for name, M0, M1 in (('crude box (0.11, 0.52)', 0.11, 0.52), ('crude box, M1 = 0.8', 0.11, 0.8), ('fine box [150,250] (0.21, 0.21)', 0.21, 0.21)):
    for gam, kap in ((1.0, 5.0), (4.0, 6.0), (4.0, 8.0), (6.0, 10.0), (8.0, 12.0)):   # kap = sup|Theta'/Theta| for Theta = eps W: nu + 4 on the outer tail (|mu - 4| on the inner)
        bc = R.box(0, 0, 0, M0=M0, M1=M1); amin = (1 + M0**2)**-2
        cc, G = F.constants(gam, bc, 'quasi', tables=T); row = []; ok = True; raw = []
        for j in range(3):
            YC0 = {int(k): float(v) for k, v in G[j]['YC'].items()}; raw.append(L4.ghat_from_YC(YC0, j, amin, kap=kap)[j])
            try:
                add, log, dis, spare = L4.l4_extra_YC(T, j, M0, M1, gam, kap=kap, s33='ibp'); ok = ok and dis <= spare*(1 + 1e-9)
                YC = {k: YC0[k] + add.get(k, 0) for k in YC0}; row.append(L4.ghat_from_YC(YC, j, amin, kap=kap)[j])
            except Exception as e:
                row.append(float('nan')); ok = False
        out[f'{name}|Gamma={gam}|kap={kap}'] = dict(Ghat_no_extras=raw, Ghat_with_L4=row, dissipation_budget_ok=bool(ok), ratio_to_28p4=[g/28.**4 for g in row], ratio_to_40p4=[g/40.**4 for g in row])
        print(f'{name}, Gamma = {gam}, kappa = {kap}: Ghat_jj (no extras) = {["%.3g" % g for g in raw]}; with L4 extras = {["%.3g" % g for g in row]} (budget ok: {ok}); Ghat_22/28^4 = {row[2]/28.**4:.3g}, /40^4 = {row[2]/40.**4:.3g}', flush=True)
json.dump(out, open('../step3_gap_v32/gonogo/t2_rates.json', 'w'), indent=1, default=float)
