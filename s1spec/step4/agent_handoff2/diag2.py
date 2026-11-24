import numpy as np
from handoff2 import *
res = {}
for k in [1120, 2240, 4480]:
    d = analyse(f'snap_extA_{k}.npz'); u = d['u']; n = len(u); w = u - d['alpha']*PSI3[:n] - d['beta2']*H2[:n]
    res[k] = (d, w)
    print(k, 'alpha', d['alpha'], 'beta2', d['beta2'], 'r2', d['r2'])
xsel = [0, 0.5, 1, 2, 3, 4, 6, 8, 10, 12, 15]
print('x   ' + '  '.join(f'w({k})/kmax' for k in res))
for x in xsel:
    i = int(round(x/HX)); print(f'{x:5.1f}', '  '.join(f"{res[k][1][i]*res[k][0]['kmax']:+.4e}" for k in res))
