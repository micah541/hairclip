import sys, numpy as np
from handoff2 import *
for f in sys.argv[1:]:
    d = analyse(f); u = d['u']; n = len(u); a = d['alpha']; b = d['beta2']
    w = u - a*PSI3[:n] - b*H2[:n]
    print(f, 'kmax', d['kmax'], 'alpha', a, 'r2', d['r2'], 'xe', d['xe'])
    cum = np.cumsum(wh[:n]*w**2)
    for x in [1,2,3,4,6,8,10,12,15,20,25]:
        i = int(x/HX)
        if i < n: print(f"  x={x:5.1f} u={u[i]:+.3e} w={w[i]:+.3e}  w*(1/a^2)={w[i]/a**2:+8.2f}  cum r2 frac={cum[i]/cum[-1]:.3f}")
