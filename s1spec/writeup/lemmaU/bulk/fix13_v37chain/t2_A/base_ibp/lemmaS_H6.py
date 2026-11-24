#!/usr/bin/env python3
"""H6 on the physically fixed strip [50,150]: K2 = L^2-normalised curvature, K3 = third derivative, from hand-off shell energies (levels 1..3) propagated by (B) (no inflow).
Local windows [y,1.5y] (conv=(1.5y/28)^gam, ell=ln1.5).  K2(y)=(L/y)(1.05U+S1), K3(y)=(L/y)^2(S2+S1+|U_ss-U_s|), L=50.   limit 1.5M=0.18.
usage: lemmaS_H6.py E1 E2 E3"""
import sys, math
E1, E2, E3 = [float(x) for x in sys.argv[1:4]]
al = 5e-4; p = 1.0246; L = 50.; eC = math.exp(0.16)   # e^{C yE^-4}, C=Ghat_33 (1e6) at yE=50 -> 0.16 (levels 1,2 smaller)
U = lambda y: 1.3056*al*y**p
def sup(Eu_a, Eu_b, ell): return math.sqrt(Eu_a/ell + 2*math.sqrt(Eu_a*(Eu_a + Eu_b)))
worst2 = worst3 = 0
for y in range(50, 151, 5):
    yb = 1.5*y; cv = yb/28.; ell = math.log(1.5)
    Eu1, Eu2, Eu3 = cv*E1*eC, cv*E2*eC, cv*E3*eC
    S1 = sup(Eu1, Eu2, ell); S2 = sup(Eu2, Eu3, ell)
    K2 = (L/y)*(1.05*U(yb) + S1); K3 = (L/y)**2*(S2 + S1 + U(yb)*abs(p*p - p))
    worst2 = max(worst2, K2); worst3 = max(worst3, K3)
print(f'hand-off E1={E1:g} E2={E2:g} E3={E3:g}: max K2 = {worst2:.3f}, max K3 = {worst3:.3f} (limit 0.18)')
