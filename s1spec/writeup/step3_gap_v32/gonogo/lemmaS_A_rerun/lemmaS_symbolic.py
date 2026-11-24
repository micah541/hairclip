#!/usr/bin/env python3
"""Symbolic check: S1-ray reduction of step3nl/graph_eq.pkl gives V_tau=-(s/4)V_s-(a V_ss+b V_s^2)_ss."""
import pickle, sympy as sp
G = pickle.load(open('/home/claude/s1spec/step3nl/graph_eq.pkl', 'rb')); rhs = G['rhs']
s = sp.Symbol('sigma')
sub = {f: {'p': 0, 'q': s, 'k': 0, 'm': 0}[f.func.__name__] for f in rhs.atoms(sp.Function) if f.func.__name__ in ('p','q','k','m')}
far = sp.simplify(rhs.subs(sub))
for e in far.free_symbols:
    if e.name == 'epsilon': far = far.subs(e, 1)
U = [f for f in far.atoms(sp.Function) if f.func.__name__ == 'U'][0]
W = sp.Function('Wf')(s); V = sp.Function('V')(s)
d1 = sp.diff(sp.simplify(far.subs(U, W).doit()), s)
for n in range(5, 1, -1): d1 = d1.subs(sp.Derivative(W, (s, n)), sp.diff(V, s, n-1))
d1 = d1.subs(sp.Derivative(W, s), V)
a = lambda x: (1+x**2)**-2; b = lambda x: -3*x*(1+x**2)**-3
print('residual:', sp.simplify(d1 - (-(s/4)*sp.diff(V, s) - sp.diff(a(V)*sp.diff(V, s, 2) + b(V)*sp.diff(V, s)**2, s, 2))))
