"""Sympy check of every pointwise identity behind the GENERAL-WEIGHT L^4 treatment (derivation: L4_general_weight.md).
z_{k+1} = d z_k / ds, Phi(s) an ARBITRARY C^2 weight, eta an arbitrary cutoff, Xi an arbitrary weight (Phi = Xi eta^2 on cutoff ramps).
 (I1)  Phi z2^4 = d/ds(Phi z1 z2^3) - 3 Phi z1 z2^2 z3 - Phi' z1 z2^3                                   [L^4 identity]
 (I2)  Phi z3 phi' = d/ds(Phi z3 phi) - Phi' z3 phi - Phi z4 phi,  phi = z2|z2|/2, phi' = |z2| z3            [{2,3,3} by one IBP]
 (I3)  (Xi eta^2)' = Xi' eta^2 + 2 Xi eta eta'                                                            [foot split of Phi']
 (I4)  -2 Xi eta eta' z3 phi = (1/3)(Xi eta eta')' |z2|^3 - (1/3) d/ds(Xi eta eta' |z2|^3)                 [exact {2,3,3} foot term; both signs of z2]
 (I5)  Phi z_{k+1}^2 = d/ds(Phi z_k z_{k+1} - Phi' z_k^2/2) - Phi z_k z_{k+2} + (1/2) Phi'' z_k^2              [interpolation step]
 (I6)  Xi = e^{-4s} W(s - c) (1 + g e^{-4s}):  Xi'/Xi = -4 + W'/W - 4 theta,  theta = g e^{-4s}/(1 + g e^{-4s}),  theta' = -4 theta (1 - theta),
       Xi''/Xi = L^2 + L'  (L = Xi'/Xi)                                                                    [log-derivatives of the Lemma-U weight]
 (I7)  W = 1/(e^{-mu x} + 1 + e^{nu x}):  W'/W = (mu p - nu q),  (ln W)'' = -[mu^2 p + nu^2 q - (mu p - nu q)^2] <= 0,  p = e^{-mu x}/Den, q = e^{nu x}/Den
       so -nu < W'/W < mu and W is log-concave                                                             [two-sided weight]
 (I8)  Theta = W D, D = 1 + g e^{-4s}:  Theta^(m)/Theta = sum_k C(m,k) (W^(k)/W) (D^(m-k)/D),  D^(k)/D = (-4)^k theta (k >= 1), m = 1..3
Inequalities used afterwards (not identities): Cauchy-Schwarz and |z1| <= M1; see L4_general_weight.md."""
import sympy as sp
s, c, g, mu, nu = sp.symbols('s c g mu nu', positive=True)
Ph = sp.Function('Phi')(s); Xi = sp.Function('Xi')(s); eta = sp.Function('eta')(s)
z1 = sp.Function('z1')(s); ph = sp.Function('phi')(s)
z = {1: z1}
for k in range(2, 6): z[k] = sp.diff(z1, s, k - 1)
D = lambda e, k=1: sp.diff(e, s, k)
ok = {}
ok['I1'] = sp.expand(Ph*z[2]**4 - (D(Ph*z[1]*z[2]**3) - 3*Ph*z[1]*z[2]**2*z[3] - D(Ph)*z[1]*z[2]**3)) == 0
ok['I2'] = sp.expand(Ph*z[3]*D(ph) - (D(Ph*z[3]*ph) - D(Ph)*z[3]*ph - Ph*z[4]*ph)) == 0
x = sp.Symbol('x', real=True); ok['phi_prime'] = sp.simplify(sp.diff(x*sp.Abs(x)/2, x) - sp.Abs(x)) == 0
ok['I3'] = sp.expand(D(Xi*eta**2) - (D(Xi)*eta**2 + 2*Xi*eta*D(eta))) == 0
for sign, name in ((1, 'I4(z2>0)'), (-1, 'I4(z2<0)')):          # |z2| = sign*z2 on the branch, phi = sign*z2^2/2
    a3 = sign*z[2]**3; phis = sign*z[2]**2/2; F = Xi*eta*D(eta)
    ok[name] = sp.expand(-2*F*z[3]*phis - (sp.Rational(1, 3)*D(F)*a3 - sp.Rational(1, 3)*D(F*a3))) == 0
for k in (1, 2, 3):
    zk, zk1, zk2 = z[k], z[k + 1], z[k + 2]
    ok[f'I5(k={k})'] = sp.expand(Ph*zk1**2 - (D(Ph*zk*zk1 - D(Ph)*zk**2/2) - Ph*zk*zk2 + D(Ph, 2)*zk**2/2)) == 0
Wf = sp.Function('W')
XiW = sp.exp(-4*s)*Wf(s - c)*(1 + g*sp.exp(-4*s)); th = g*sp.exp(-4*s)/(1 + g*sp.exp(-4*s))
L = sp.diff(XiW, s)/XiW
ok['I6a'] = sp.simplify(L - (-4 + sp.diff(Wf(s - c), s)/Wf(s - c) - 4*th)) == 0
ok['I6b'] = sp.simplify(sp.diff(th, s) + 4*th*(1 - th)) == 0
ok['I6c'] = sp.simplify(sp.diff(XiW, s, 2)/XiW - (L**2 + sp.diff(L, s))) == 0
Den = sp.exp(-mu*x) + 1 + sp.exp(nu*x); W = 1/Den; p, q = sp.exp(-mu*x)/Den, sp.exp(nu*x)/Den
ok['I7a'] = sp.simplify(sp.diff(W, x)/W - (mu*p - nu*q)) == 0
ok['I7b'] = sp.simplify(sp.diff(sp.log(W), x, 2) + (mu**2*p + nu**2*q - (mu*p - nu*q)**2)) == 0
# -(ln W)'' = Var of the exponent X in {-mu, 0, nu} with probabilities (p, 1-p-q, q)  ->  >= 0 (log-concave); explicit: = p(1-p)mu^2 + q(1-q)nu^2 + 2pq mu nu >= 0
ok['I7c'] = sp.simplify((mu**2*p + nu**2*q - (mu*p - nu*q)**2) - (p*(1 - p)*mu**2 + q*(1 - q)*nu**2 + 2*p*q*mu*nu)) == 0
Dd = 1 + g*sp.exp(-4*s); Wg = Wf(s - c); Th = Wg*Dd
for m in (1, 2, 3):
    rhs = sum(sp.binomial(m, k)*(sp.diff(Wg, s, k)/Wg)*(sp.diff(Dd, s, m - k)/Dd) for k in range(m + 1))
    ok[f'I8(m={m})'] = sp.simplify(sp.diff(Th, s, m)/Th - rhs) == 0
    if m >= 1: ok[f'I8D(k={m})'] = sp.simplify(sp.diff(Dd, s, m)/Dd - (-4)**m*th) == 0
for k, v in ok.items(): print(f'{k:10s} {"verified" if v else "FAILED"}')
assert all(ok.values())
print('all general-weight identities verified (I1-I8)')
