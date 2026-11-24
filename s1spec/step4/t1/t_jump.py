import numpy as np
from flint import arb
from base2 import Base; from arbops import SlabPolish
B = Base('P11'); i = 12; R, D = B.rep(i); q = B.slabs[i]
SP = SlabPolish(R, D, B.RD, q['nodes'], q['dt'])
Y, Yd = SP.Y_at(0.0); els, gl = SP.op.elements(Y, Yd)
jP, jk = [], [[] for _ in range(5)]
for e in range(len(els) - 1):
    a, b = els[e], els[e + 1]
    jP.append(float(abs(a['Phi'](arb(1)) - b['Phi'](arb(0))).mid()))
    pa, pb = a['kp'], b['kp']
    for m in range(5):
        va = pa(arb(1)) / a['h']**m; vb = pb(arb(0)) / b['h']**m
        jk[m].append(float(abs(va - vb).mid())); pa = pa.derivative(); pb = pb.derivative()
print('max jump of Phi at nodes:', max(jP), ' at node', int(np.argmax(jP)) + 1)
for m in range(5): print(f'max jump of kappa^({m}):', max(jk[m]), ' at node', int(np.argmax(jk[m])) + 1)
print('kappa\'(0), kappa\'\'\'(0):', float((els[0]['ku'](arb(0)) / els[0]['h']).mid()), float((els[0]['ku'].derivative().derivative()(arb(0)) / els[0]['h']**3).mid()))
