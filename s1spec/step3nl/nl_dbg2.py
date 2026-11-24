import numpy as np, sys
import nlflow
from nlflow import Flow
fl = Flow(X=34.0, h=0.05); x = fl.x
cut = lambda x: np.where(x < 8, 1.0, np.where(x < 10, np.cos(np.pi / 4 * (x - 8))**2, 0.0))
u0 = 0.2 * np.exp(-x**2 / 4) * cut(x)
# monkeypatch: record diagnostics
orig = fl.rhs_parts
def rp(v):
    out = orig(v)
    return out
H, v = nlflow.run(fl, u0, dt=2e-3, T=3.5, every=0.5, verbose=False)
J0, J1, J2, J3, J4 = fl.jet(v)
i = np.argmax(np.abs(J0)); print("sup|u| =", np.abs(J0).max(), "at x =", x[i])
print("max|u'| =", np.abs(J1).max(), "at x =", x[np.argmax(np.abs(J1))])
print("max|u''| =", np.abs(J2).max(), "at x =", x[np.argmax(np.abs(J2))])
print("min alpha =", (1 - fl.P[2] * J0).min())
for xx in [0, 2, 5, 10, 15, 20, 25, 30, 33]:
    j = int(round(xx / fl.h)); print(f"x={xx:5.1f} u={J0[j]:+.4e} u'={J1[j]:+.4e} u''={J2[j]:+.4e} u''''={J4[j]:+.3e}")
