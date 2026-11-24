"""Fully rigorous (but crude) bounds for int V_i^2 that use ONLY the upper barrier Ly<=L(tau): mu V0^2 <= Ly, mu epsI V2^2 <= Ly.
Lambda_int = sum_k Lup_k D_k + L*/rho (tail: Ly <= L* exp(-rho tau)).  Usage: python3 crude_int.py rate_table_X.json"""
import sys, json
from flint import arb, ctx
ctx.prec = 120
import constants as K
A = lambda v: arb(str(v))
r = json.load(open(sys.argv[1]))
lam = A(0)
for st in r['steps']: lam = lam + A(st['Ly_hi']) * A(st['D'])
t = r['tail']; lam = lam + A(t['L_star']) / A(t['rho']) * (1 + A('1e-9'))
mu = A(K.MU); eI = A(K.EPSI); C1W = A(K.C1W)
I0 = lam / mu; I2 = lam / (mu * eI); l = (I2 / I0).sqrt()
I1 = C1W * I0 + (l * I0 + I2 / l) / 2
f = lambda x: float(x.upper())
print(f"{sys.argv[1]}: Lambda_int<={f(lam):.4e}  int V0^2<={f(I0):.4e}  int V1^2<={f(I1):.4e}  int V2^2<={f(I2):.4e}")
