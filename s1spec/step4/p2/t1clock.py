"""Theorem 4(a) (review B7): the hand-off time t_1 (since Gamma0) and k_max of the base at t_1, for tau* anywhere in slab SLAB.
Physical time of the TRUE solution: dt = lam^4 (1 - theta) dtau (Theorem P2; lam the true scale, theta the clock modulation), and
lam^4 = lam~^4 (lam/lam~)^4 with lam~ = e^{ln lam~} the base scale.  Per slab, tau runs over 0.15 (t_c in [-1, 1], dtau/dt_c = 0.075):
  int_slab lam~^4 dtau in Dtau_s [e^{4 lo_s}, e^{4 hi_s}]   (T6 range_lnlam, outward);  slab 2 from the switch t_c,a on, slab SLAB a fraction in [0, 1].
Clock factor: |theta| = |x_3|/n <= theta_s := xbar_s / nV_s, with (M) |x| <= xbar_s := C_x Ybar + x_r,s + x_NL (T5 (M), T7 4c(M)):
  C_x from the t7 block of CONSTS_JSON, x_r,s = mbar_s |c_phi|_2 R2_s per slab (as n3far.xr_per_slab: T7 C3 mbar, the largest 2-norm of
  the T7 certificates' c_phi, T6 R2), x_NL = T7 meta.modulation_nl.xNL_box_max; nV_s from CONSTS_JSON.
Scale factor (review O6a): the base identity gives d ln lam~/dtau = a~ - 1/4 = -c and the true one d ln lam/dtau = a - 1/4, and
  x_1 = (a - a~) - theta c, so d ln(lam/lam~)/dtau = x_1 + theta c and, per slab s (cumulative to the end of the slab),
  |ln(lam/lam~)| <= rho_s := rho_a + sum_{s' <= s} Dtau_s' (xbar_s' + theta_s' cbar_s'),   cbar_s = max |T6 range_c| on slab s;
  rho_a = mu_1 + s_clk cbar_2 bounds ln(lam/lam~) at tau_a from the switch's IFT (switch_s2.json 'ift': scale shift mu1, clock shift s),
  whether the switch compares with the base at tau_a or at tau_a + s_clk.  The latter also moves the start of the integral: the term
  [-1, 1] s_clk sup_slab2 lam^4 (1 - theta) is added (it covers both conventions).
t_1 = t_96 + sum_s int_slab lam~^4 dtau * e^{[-4 rho_s, 4 rho_s]} * (1 + [-theta_s, theta_s]) (+ the start term); every factor is positive.
k_max: the base tip curvature hull over slab SLAB divided by lam~ (the base curve, in its own scale).
usage: python3 -I t1clock.py CONSTS_JSON [SLAB] [YBAR]   (env T6DIR: residuals and jets; env SWITCHJSON: default p2/out/switch_s2.json)"""
import sys, os, json, math
from flint import arb, ctx
ctx.prec = 128
S4 = '/home/claude/s1spec/step4'
T6D = os.environ.get('T6DIR', f'{S4}/p2/t6')                                # env T6DIR: residuals and jets of one T6 data set
T6 = json.load(open(f'{T6D}/residuals.json')); cons = json.load(open(sys.argv[1]))
SLAB = int(sys.argv[2]) if len(sys.argv) > 2 else 180; YB = arb(sys.argv[3] if len(sys.argv) > 3 else '1e-12')
T7P = f'{S4}/t7/constants.json'; t7c = json.load(open(T7P))
SWP = os.environ.get('SWITCHJSON', f'{S4}/p2/out/switch_s2.json'); sw = json.load(open(SWP))
U = lambda v: math.nextafter(float(v.upper()), math.inf); Lw = lambda v: math.nextafter(float(v.lower()), -math.inf)
t96 = arb('4.597850539485830831', '7.8e-19'); tca = arb(982986695892350) / arb(2) ** 50
assert sw['raw']['meta']['slab'] == 2 and sw['raw']['meta']['tc'] == float(tca.mid()) and sw['ift']['inside_apriori']   # this switch run
rec = {r['slab']: r for r in T6['per_slab']}
pc = {p['slab']: p for p in cons['per_slab']}
def lam4(s):
    lo, hi = rec[s]['scalars']['range_lnlam']; return (4 * arb(lo)).exp().union((4 * arb(hi)).exp())
CX = arb(repr(cons['t7']['Cx'])); XNL = arb(repr(t7c['meta']['modulation_nl']['xNL_box_max']))
CPHI = max((arb(repr(sum(v * v for v in c['c_phi']))).sqrt() for c in t7c['certificates']), key=lambda v: float(v.upper()))
def xr(s): return arb(repr(t7c['slabs'][str(s)]['C3']['mbar'])) * CPHI * arb(repr(rec[s]['norms']['R2']))
def cbar(s): return arb(max(abs(v) for v in rec[s]['scalars']['range_c']))
MU1, SCLK = arb(repr(sw['ift']['mu1'])), arb(repr(sw['ift']['s']))
rho = MU1 + SCLK * cbar(2); rho_a = rho
dt_base = arb(0); dt = arb(0); th_max = 0.0; xb_max = 0.0; xr_max = 0.0; per = []
for s in range(2, SLAB + 1):
    dtau = arb('0.075') * (1 - tca) if s == 2 else arb('0.15')                 # tau-length of the part of slab s on [tau_a, tau*] (at most)
    xb = CX * YB + xr(s) + XNL; th = xb / arb(repr(pc[s]['norms']['nV']))
    rho = rho + dtau * (xb + th * cbar(s))                                      # |ln(lam/lam~)| on slab s (to its end)
    frac = arb('0.5', '0.5') if s == SLAB else arb(1)                           # tau* anywhere in slab SLAB
    ib = dtau * lam4(s) * frac
    dt_base += ib; dt += ib * (4 * arb(0, rho.upper())).exp() * (1 + arb(0, th.upper()))
    th_max = max(th_max, U(th)); xb_max = max(xb_max, U(xb)); xr_max = max(xr_max, U(xr(s)))
    if s == 2: start = SCLK * (4 * arb(rec[2]['scalars']['range_lnlam'][1])).exp() * (4 * arb(0, rho.upper())).exp() * (1 + arb(0, th.upper()))
    if s in (2, 3, SLAB) or s % 40 == 0: per.append(dict(slab=s, theta=U(th), rho=U(rho), xbar=U(xb)))
t1 = t96 + dt + start * arb(0, 1)
t1_old = t96 + dt_base * (1 + arb(0, th_max))                                   # the former enclosure (no lam/lam~, no x_NL), for comparison
lo, hi = rec[SLAB]['scalars']['range_lnlam']
sys.path.insert(0, f'{S4}/t6'); import reader as RD
J = RD.jets(SLAB, d=f'{T6D}/jets'); k0 = J['k0']
ktip = arb(float(k0[0, 0, 0])).union(arb(float(k0[0, 0, 1]))); kmax = arb(float(abs(k0).max()))
out = dict(args=dict(argv=sys.argv[1:], consts=sys.argv[1], slab=SLAB, ybar=float(YB.mid()), t6=T6D, t7=T7P, switch=SWP),
           slab=SLAB, t1=[Lw(t1), U(t1)], dt=[Lw(dt_base), U(dt_base)], dt_true=[Lw(dt), U(dt)], start_term=U(start),
           theta_max=th_max, xbar_max=xb_max, x_r_max=xr_max, x_NL=U(XNL), C_x=U(CX), c_phi=U(CPHI),
           lnlam_ratio=dict(at_tau_a=U(rho_a), on_slab=U(rho), mu1=U(MU1), s_clk=U(SCLK)), per_slab=per,
           t1_without_scale_factor=[Lw(t1_old), U(t1_old)],
           lnlam=[lo, hi], base_kmax=[Lw(ktip / arb(hi).exp()), U(kmax / arb(lo).exp())])
json.dump(out, open(f'{S4}/p2/out/t1clock.json', 'w'), indent=1); print(json.dumps(out))
