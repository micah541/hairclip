"""G0_INDEP: arb (outward) evaluation of the sizes quoted in G0_INDEP.md sec. 2.  The COEFFICIENTS fed in are hand-derived from certified
sups of the proof tree (tags in the report); the arithmetic here is certified (python-flint arb, 200 bit), sup over t by ball subdivision."""
import json
from flint import arb, ctx
ctx.prec = 200
A = lambda s: arb(str(s))
up = lambda x: float(x.upper())
MU, EPSI, MU1, MU2, C1W, PSI1 = A('0.08'), A('3e-5'), A('0.1552'), A('0.6555'), A('0.8477'), A('1.305')
out = {}
# (1) beta coupling at level 0:  mu*cb*(b0+b1 al) A V0 ; A = s a_, V0 = s v_/sqrt(mu)  ->  cross coefficient kappa a_ v_, max over the circle kappa/2
def kap(al):
    al = A(al); return (A('0.2596') + 2 * al * A('1.796')) * (A('0.00244') + A('66.53') * al) * MU.sqrt()
rho = A('0.0026930')
for al in ('1.5e-4', '3.6849e-5', '0'):
    k = kap(al); out[f'beta_cross_kappa_al={al}'] = up(k); out[f'beta_cross_half_over_rho_al={al}'] = up(k / 2 / rho)
    out[f'beta_second_order_shift_al={al}'] = up(k * k / (4 * (MU1 - A('0.0152'))))      # eigenvalue shift of diag(-r_A,-r_V)+cross, r_V-r_A >= mu1-0.0152
# degree-0 part of beta: 0.657 Fn(al) + STRIP_BETA (1.06e-9) -> forcing mu*cb*(...)*V0 ; Young against mu*MU1*V0^2 -> absolute term per unit tau
def Fn(al): al = A(al); return A('7.9326e-7') * al + A('2.7694e-6') * al**2 + A('12.1044') * al**3
def FnW(al): al = A(al); return A('1.093e-7') * al + A('6.5e-9') * al**2 + A('0.1823') * al**3
for al in ('1.5e-4', '3.6849e-5'):
    cb = A('0.2596') + 2 * A(al) * A('1.796')
    out[f'I0_forcing_coded_FnW*1.67_al={al}'] = up(FnW(al) * A('1.67')); out[f'I0_forcing_missing_beta_Fn_al={al}'] = up(cb * A('0.657') * Fn(al))
strip = A('0.2596') * A('1.06e-9'); out['I0_forcing_missing_beta_strip_const'] = up(strip)
F0add = MU * strip**2 / (4 * MU1); out['floor_add_from_beta_strip (compare F0=1.312e-21)'] = up(F0add); out['ratio_to_F0'] = up(F0add / A('1.312192e-21'))
out['sqrtLF_with_beta_strip'] = up((2 * (A('1.312192e-21') + F0add) / rho).sqrt())
# (2) level-0 top-order cross terms: [(xa1+xb1 s) V1 V2 + (xa0+xb0 s) V0 V2] / (mu1 V0^2 + V2^2/2), V2 = t V0, V1 <= V0 sqrt(t + c1W)
def sup_ratio(c1, c0, den0, T=64, n=6400):
    best = arb(0)
    for i in range(n):
        t = arb((i + 0.5) * T / n, 0.5 * T / n)
        r = (c1 * t * (t + C1W).sqrt() + c0 * t) / (den0 + t * t / 2)
        if r.upper() > best.upper(): best = r.upper()
    tail = 2 * (c1 * (1 + C1W / T).sqrt() / A(T).sqrt() + c0 / T)       # t >= T
    return max(up(best), up(tail))
for tag, d001, s in (('typed box (alpha<=1e-3,|w_x|<=0.02): d001=0.0227, sup|w_xx|=0.01875', '0.0227', '0.01875'),
                     ('est. actual box (alpha<=1.5e-4, |w_x|<=0.0153): d001~2.2e-3', '0.0022', '0.01875')):
    d = A(d001); s = A(s); cu1 = 4 * (A('0.085') * (A('1') if d001 == '0.0227' else A('0.15')) + (A('0.02') if d001 == '0.0227' else A('0.0153')))   # |d c/d u_x| <= 4|u_x|
    c1 = 2 * d * PSI1 + cu1 * (A('1.54e-3') * (A('1') if d001 == '0.0227' else A('0.15')) + s); c0 = d * A('1.6954') + PSI1 * cu1 * (A('1.54e-3') + s)
    out[f'lvl0_cross_coeffs [{tag}]'] = [up(c1), up(c0)]; out[f'lvl0_cross_over_dissipation [{tag}]'] = sup_ratio(c1, c0, MU1)
# (3) level-2 cross: d001 V4 (2.61 V3 + 1.7 V2) / (V4^2/2 + mu2 V2^2), V4 = t V2, V3 <= V2 sqrt(t + c1W)
for d001 in ('0.0227', '0.0022'):
    d = A(d001); best = 0.0
    for i in range(6400):
        t = arb((i + 0.5) * 64 / 6400, 0.5 * 64 / 6400)
        r = d * t * (A('2.61') * (t + C1W).sqrt() + A('1.7')) / (t * t / 2 + MU2); best = max(best, up(r))
    out[f'lvl2_cross_over_dissipation d001={d001}'] = best
# (4) level-2 forcing: coded |al|^3 vs needed 2.701e-8|al| + 1.156e-7 al^2 + 2.603e-3|al|^3
for al in ('1.5e-4', '3.6849e-5', '6.9e-7'):
    a = A(al); out[f'lvl2_forcing_needed_over_coded_al={al}'] = up((A('2.701e-8') * a + A('1.156e-7') * a * a + A('2.603e-3') * a**3) / a**3)
out['lvl2_forcing_abs_Young (mu epsI l^2/(2*0.95)) al=1.5e-4'] = up(MU * EPSI * (A('2.701e-8') * A('1.5e-4'))**2 / A('1.9'))
# (5) relative down-rounding of typed constants (typed/certified - 1)
pairs = dict(AF=('6.9555e-2', '0.06955503763448806'), CF0=('0.4421', '0.44214268661599576'), CF2=('0.4235', '0.4235276571468286'), CF3=('0.4250', '0.42501723188074436'),
             DS1=('1.5418e-2', '0.015418374795665034'), DS2=('2.3648e-2', '0.02364817888631596'), SB=('3.5955e-5', '3.5955166253940984e-05'), SA=('1.3372e-4', '0.0001337202808774133'),
             IA=('9.5997e-12', '9.599705423877498e-12'), d010=('1.67e-3', '0.0016711395692254884'), d110=('1.086', '1.0864220764487984'), d300=('3.314', '3.314008813351393'),
             FN0=('7.93e-7', '7.9325012815966e-07'), FN1=('2.769e-6', '2.7693212788311072e-06'), C3=('12.104', '12.104385241866112'), FNW0=('1.09e-7', '1.092e-7'), d100_ramp=('1.1e-4', '1.129e-4'))
out['typed_over_certified_minus_1'] = {k: float((A(t) / A(c) - 1).mid()) for k, (t, c) in pairs.items()}
json.dump(out, open('g0_numbers.json', 'w'), indent=1)
for k, v in out.items(): print(k, '=', v)
