"""[est, float] size of the terms of the E0 identity that come from the NONLINEAR DRIFT  -(q/4) u' (1/alpha_g - 1),  alpha_g = 1 - k u,
which is not part of T_c (nlin_E0.py: curvature part only) and hence not in alin_E0 / cubic_E0fine.
alpha-linear part:  -(qk/4)(u0 w)'  paired with rho w  =  (alpha/8) int rho w^2 d(x),   d = ((qk)' - phi' qk) psi - qk psi',   (qk)' = (1+kp)k + q m.
Also the w-dependent parts of a A1(u) + b B1(u):  |a| (1 + sup|q' - phi' q|/2) A^2,  |b| sup|phi' cos + k sin|/2 A^2  (coefficients only).
Float profile: step3nl_v4/agent_mod/prof_float.npz, psi3^: step3nl_v4/psi3md2_14.json.  Not certified."""
import json, os, numpy as np
from numpy.polynomial import chebyshev as Ch
V4 = '/home/claude/s1spec/step3nl_v4'
d = np.load(V4 + '/agent_mod/prof_float.npz'); x, p, q, k, m, th = [d[n] for n in ('x', 'p', 'q', 'k', 'm', 'th')]
Jf = json.load(open(V4 + '/psi3md2_14.json'))
def ev(xv, dn):
    o = np.zeros_like(xv)
    for pc in Jf['pieces']:
        a_, b_ = pc['a'], pc['b']; msk = (xv >= a_) & (xv <= b_); cc = np.array(pc['coef'])
        for _ in range(dn): cc = Ch.chebder(cc)
        o[msk] = Ch.chebval(1 - 2 * (xv[msk] - a_) / (b_ - a_), cc) * (-2 / (b_ - a_)) ** dn
    return o
msk = x <= 31.25
psi, psi1 = ev(x, 0), ev(x, 1)
E_ = (1 + x * x) ** (-1 / 3); ph1 = x * E_ / 3
qk = q * k; qk1 = (1 + k * p) * k + q * m
dd = ((qk1 - ph1 * qk) * psi - qk * psi1) / 8
rho = np.exp(-(1 + x * x) ** (2 / 3) / 4)
i = np.argmax(np.abs(dd[msk]))
res = dict(sup_abs_d=float(np.abs(dd[msk]).max()), argmax_x=float(x[msk][i]), sup_qk=float(np.abs(qk[msk]).max()), psi_tip=float(psi[0]),
           a_w_coef=float(1 + 0.5 * np.abs((1 + k * p) - ph1 * q)[msk].max()), b_w_coef=float(0.5 * (np.abs(ph1) + np.abs(k))[msk].max()),
           cubic_drift_coef=float(np.abs((qk1 - ph1 * qk) / 4)[msk].max() / 3))
res['alpha_form_allowance'] = float(np.ceil(2 * res['sup_abs_d'] * 10) / 10)      # x2 safety on a float sup
res['note'] = '[est] float; alpha_form_allowance = 2 x sup|d| rounded up to 0.1, used as the drift allowance in e0_cert.py (7)'
json.dump(res, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'e0_drift_est.json'), 'w'), indent=1)
print(json.dumps(res, indent=1))
