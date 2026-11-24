import numpy as np
from scipy.optimize import brentq
from evans import evans
for which, idx in (('even',0),('odd',1)):
    xs = np.linspace(-2.0, 3.0, 201)
    vals = [evans(x)[idx].real for x in xs]
    roots = []
    for a,b,fa,fb in zip(xs[:-1],xs[1:],vals[:-1],vals[1:]):
        if fa == 0: roots.append(a)
        elif fa*fb < 0:
            roots.append(brentq(lambda x: evans(x)[idx].real, a, b, xtol=1e-14))
    print(which, "real zeros in [-2,3]:", [f"{r:+.10f}" for r in roots], flush=True)
    # imaginary part check on the real axis
    print(which, "max |Im D|/|D| on real axis:", max(abs(evans(x)[idx].imag)/abs(evans(x)[idx]) for x in xs[::20]))
