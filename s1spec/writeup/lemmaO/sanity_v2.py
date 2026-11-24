#!/usr/bin/env python3
"""PDE sanity check for lemmaO_consts_v2 (cylinder [-X,X], periodic box [-6,6), Fourier+BDF)."""
import math, numpy as np
from scipy.integrate import solve_ivp
from lemmaO_consts_v2 import constants
L = 12.0; N = 192
xg = -L/2 + L*np.arange(N)/N; dx = L/N
kk = 2*np.pi*np.fft.fftfreq(N, d=dx)
D = lambda u, m=1: np.real(np.fft.ifft((1j*kk)**m*np.fft.fft(u)))
a = lambda p: (1+p**2)**-2; b = lambda p: -3*p*(1+p**2)**-3
rhs = lambda t, h: -D(a(D(h))*D(h,3) + b(D(h))*D(h,2)**2)
def run(X, modes, Tend, slope):
    h0 = sum(c*np.cos(2*np.pi*k*xg/L + ph) for k, c, ph in modes); h0 *= slope/np.max(np.abs(D(h0)))
    ins = np.abs(xg) <= X
    M = max([np.max(np.abs(D(h0, j)[ins])) for j in range(1, 5)] + [math.sqrt(np.sum(D(h0, 5)[ins]**2)*dx)])
    ts = np.concatenate([[0], np.geomspace(Tend*1e-4, Tend, 40)])
    H = solve_ivp(rhs, (0, Tend), h0, method='BDF', t_eval=ts, rtol=1e-11, atol=1e-13).y.T
    K2 = max(np.max(np.abs(D(h, 2)[ins])) for h in H); K3 = max(math.sqrt(np.sum(D(h, 3)[ins]**2)*dx) for h in H)
    assert max(np.max(np.abs(D(h)[ins])) for h in H) <= 0.25
    c = constants(K2, K3, M, X, t1=Tend, w=0.2 if X < 3 else 0.5); inn = np.abs(xg) <= 1
    worst = 0
    for t, h in zip(ts[1:], H[1:]):
        g = h - h0
        for obs, bnd in [(np.max(np.abs(g[inn])), c['C0']*math.sqrt(t)), (np.max(np.abs(D(g)[inn])), c['C1']*math.sqrt(t)),
                         (np.max(np.abs(D(g, 2)[inn])), c['C2']*math.sqrt(t)), (np.sum(D(g, 3)[inn]**2)*dx, c['C3']*t)]:
            worst = max(worst, obs/bnd)
    print(f'X={X} M={M:.3f} K2={K2:.3f} K3={K3:.3f}: worst observed/bound = {worst:.2e}'); assert worst < 1
run(2.4, [(1, 1.0, .3), (2, .5, 1.), (3, .3, 2.)], 1e-3, 0.1)
run(5.0, [(2, 1.0, .3), (3, .5, 1.), (4, .3, 2.)], 1e-3, 0.1)
run(2.4, [(k, k**-5.0, 1.7*k) for k in range(1, 30)], 1e-4, 0.1)
print('OK')
