"""Theorem 4(a) (review B7): the hand-off time t_1 (since Gamma0) and k_max of the base at t_1, for tau* anywhere in slab SLAB.
Physical time: dt = lam^4 dtau (rescaled clock tau, lam = e^{ln lam}); per slab tau runs over 0.15 (t_c in [-1, 1], dtau/dt_c = 0.075) and the
certified range of ln lam (T6 range_lnlam, outward) gives int lam^4 dtau in 0.15 [e^{4 lo}, e^{4 hi}].  The switch is at t_c,a in slab 2,
where the physical time is t_96 (P1).  The clock modulation of the true solution, |theta| <= (C_x Ybar + x_r)/nV per slab (T5 (M), as
n3far.py), moves t_1 by at most theta_max (t_1 - t_96).  k_max: the base tip curvature hull over slab SLAB divided by lam.
usage: python3 -I t1clock.py CONSTS_JSON [SLAB] [YBAR]"""
import sys, json, math
from flint import arb
S4 = '/home/claude/s1spec/step4'
T6 = json.load(open(f'{S4}/p2/t6/residuals.json')); cons = json.load(open(sys.argv[1]))
SLAB = int(sys.argv[2]) if len(sys.argv) > 2 else 180; YB = arb(sys.argv[3] if len(sys.argv) > 3 else '1e-12')
U = lambda v: math.nextafter(float(v.upper()), math.inf); Lw = lambda v: math.nextafter(float(v.lower()), -math.inf)
t96 = arb('4.597850539485830831', '7.8e-19'); tca = arb(982986695892350) / arb(2) ** 50
rec = {r['slab']: r for r in T6['per_slab']}
def lam4(s):
    lo, hi = rec[s]['scalars']['range_lnlam']; return (4 * arb(lo)).exp().union((4 * arb(hi)).exp())
part2 = arb('0.075') * (1 - tca) * lam4(2)                                  # slab 2 after the switch
full = sum((arb('0.15') * lam4(s) for s in range(3, SLAB)), arb(0))
last = arb('0.15') * lam4(SLAB) * arb('0.5', '0.5')                        # tau* anywhere in slab SLAB: a fraction in [0, 1] of it
dt = part2 + full + last
XR = arb('2.0e-16')                                                         # x_r = mbar c_phi R2 <= 2.0e-16 on slabs 2-180 (T8, review B6)
th = max(U((arb(repr(cons['t7']['Cx'])) * YB + XR) / arb(repr(p['norms']['nV']))) for p in cons['per_slab'] if 2 <= p['slab'] <= SLAB)
t1 = t96 + dt * (1 + arb(0, th))
lo, hi = rec[SLAB]['scalars']['range_lnlam']
sys.path.insert(0, f'{S4}/t6'); import reader as RD
J = RD.jets(SLAB, d=f'{S4}/p2/t6/jets'); k0 = J['k0']
ktip = arb(float(k0[0, 0, 0])).union(arb(float(k0[0, 0, 1]))); kmax = arb(float(abs(k0).max()))
kb = (ktip / arb(hi).exp()).union(kmax / arb(lo).exp())
out = dict(args=dict(argv=sys.argv[1:], consts=sys.argv[1], slab=SLAB, ybar=float(YB.mid()), t6='p2/t6/residuals.json, p2/t6/jets'), slab=SLAB, t1=[Lw(t1), U(t1)], dt=[Lw(dt), U(dt)], theta_max=th, lnlam=[lo, hi], base_kmax=[Lw(ktip / arb(hi).exp()), U(kmax / arb(lo).exp())])
json.dump(out, open(f'{S4}/p2/out/t1clock.json', 'w'), indent=1); print(json.dumps(out))
