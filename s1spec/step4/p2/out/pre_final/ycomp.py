"""P2 (Y) comparison majorant, certified in arb (outward), per slab of the v2 base.
(Y): 1/2 Lambda' <= -rho Lambda - kD*Dres + a_c Y + Xi'^2  (kD = 1/2: half of T7's reserve D is kept for Lemma J),
residual terms of T5 (E0),(E2),(I0),(I2) [see STEP4_RESULT.md, P2 section] with interpolations W1 <= (1/2+sqrt c1)A + B/2 etc.,
rho_m = rho_c(interval) - dT6 - dN.  Constant majorant: Ybar := Y_c >= max(Y(tau0), max_m(a_c/rho_m + Xi'_m/sqrt(rho_m)), 1e-12).
usage: python3 ycomp.py [Y0]   -> out/ycomp.json"""
import math, sys, json, os
from flint import arb, ctx
ctx.prec = 128
S4 = '/home/claude/s1spec/step4'
T6 = json.load(open(os.environ.get('T6RES', f'{S4}/t6/residuals.json'))); T7 = json.load(open(f'{S4}/t7/constants.json'))
SLAB0 = int(os.environ.get('SLAB0', '0'))
def up(x): return math.nextafter(float(arb(x).upper()), math.inf)   # float upper bound, rounded outward
def U(v): return arb(v) if isinstance(v, arb) else arb(repr(float(v)))  # stored floats are certified sups: take them exactly
EPS, EPSI, C1 = arb('1e-5'), arb('3e-5'), arb('0.5804')
PHI1, PHI2 = arb('1.083'), arb('1.3')
CHI1, CHI2 = arb('0.73'), arb('0.84')             # sup|chi_32'|, sup|chi_32''| (S'max 2.1875/3, S''max 7.53/9)
cl = T7['meta']['closure']; MUE = U(cl['mu_e']); TH0 = arb('0.025')
c1W = U(T7['meta']['weight_consts']['c1W']); PSI1 = U(T7['meta']['weight_consts']['Psi1']); PSI2 = 2 * c1W
cphi = None
for c in T7['certificates']:
    v = c['c_phi']; cphi = v if cphi is None else [max(a, b) for a, b in zip(cphi, v)]
CPHI = sum((U(v) ** 2 for v in cphi), arb(0)).sqrt()
DN = arb('1e-6')                                  # nonlinear allowance (STEP4_RESULT, P2 section): actual box |e|,|e'| <= 1e-8
VARTH = arb('1.001')                              # |1 - vartheta| <= 1 + 1e-3
sc1, sc1W = C1.sqrt(), c1W.sqrt()
kW1A, kW1B = arb(1) / 2 + sc1, arb(1) / 2         # W1 <= kW1A A + kW1B B
kV1 = arb(1) / 2 + sc1W                           # V1 <= kV1 V0 + V2/2 ; V3 <= kV1 V2 + V4/2 ; W3 <= kW1A B + D4/2
# (rho chi u)''/rho = chi u'' + 2(chi' + chi rho'/rho) u' + (chi'' + 2 chi' rho'/rho + chi rho''/rho) u
cW1 = 2 * (CHI1 + PHI1); cA = CHI2 + 2 * CHI1 * PHI1 + PHI2
D4CUB = arb('3.0e-9')                                        # gap 9 (T7 Lemma 4a): eps D4^2 slice reserved for m23 |e''|^3
shares = dict(B=TH0 / 8, D4=EPS * (arb(1) / 16 - D4CUB), V2=MUE / 16, V4=MUE * EPSI / 16)
farg = {int(m): v['farg'] for m, v in T7['slabs'].items()}          # gap 10: raw farg values republished in t7/constants.json (SE-corrected base)
Y0 = arb(sys.argv[1]) if len(sys.argv) > 1 else arb(0)
rows = []; worst = arb(0)
for rec in T6['per_slab']:
    m = rec['slab']; n = rec['norms']
    if m < SLAB0: continue
    s7 = T7['slabs'][str(m)]; cid = s7['cert']['id']
    R2, Rs, R1, RIs = (VARTH * U(n[k]) for k in ('R2', 'Rs', 'R1', 'RIs'))
    mbar = U(s7['C3']['mbar']); g = farg[m]
    g0, g2, g2E = U(g['sqrt_lmax_G0']), U(g['sqrt_lmax_G2']), U(g['sqrt_lmax_G2E'])
    xr = mbar * CPHI * R2                          # residual part of the modulation, |x_r| <= mbar |c_phi| R2
    rB = (1 + cW1 * kW1B) * R2 + EPS * Rs * (2 * PHI1 * kW1A + PHI2) + EPS * g2E * xr
    rD4 = EPS * Rs * (1 + PHI1)
    rI2 = MUE * EPSI * (RIs + g2 * xr)
    rV2 = MUE * R1 * (1 + PSI1) + rI2 * (2 * PSI1 * kV1 + 2 * c1W)
    rV4 = rI2 * (1 + PSI1)
    rA = (cA + cW1 * kW1A) * R2 + arb(3).sqrt() * xr
    rV0 = MUE * (R1 * (2 * PSI1 * kV1 + PSI2) + g0 * xr)
    Xi2 = rB ** 2 / (4 * shares['B']) + rD4 ** 2 / (4 * shares['D4']) + rV2 ** 2 / (4 * shares['V2']) + rV4 ** 2 / (4 * shares['V4'])
    ac = rA + rV0 / MUE.sqrt()
    rho = U(cl['rho_c'][str(cid)]) - U(T7['meta']['T6_sensitivity']['per_interval'][str(cid)]) - DN
    lev = ac / rho + Xi2.sqrt() / rho.sqrt()
    worst = max(worst, lev, key=lambda z: up(z))
    rows.append(dict(slab=m, k=rec['k'], cert=cid, rho=float(rho.lower()), a_c=up(ac), Xi=up(Xi2.sqrt()), level=up(lev), xr=up(xr)))
Yc = max(worst, Y0, arb('1e-12'), key=lambda z: up(z))
out = dict(Y0=up(Y0), Yc=up(Yc), level_max=up(worst), argmax=max(rows, key=lambda r: r['level'])['slab'], rows=rows,
           consts=dict(mu_e=float(MUE.mid()), theta0=0.025, cphi=up(CPHI), dN=1e-6, Psi1=up(PSI1), Psi2=up(PSI2), c1W=up(c1W)))
out['T6RES'] = os.environ.get('T6RES', 't6/residuals.json'); out['SLAB0'] = SLAB0
json.dump(out, open(os.environ.get('YOUT', f'{S4}/p2/out/ycomp.json'), 'w'), indent=1)
print(f"Y_c = {out['Yc']:.4e}  (residual level max {out['level_max']:.4e} at slab {out['argmax']}; Y0 {out['Y0']:.3e}); "
      f"level at first/10/100/last: " + ', '.join(f"{rows[i]['level']:.3e}" for i in (0, 10, 100, -1)))
