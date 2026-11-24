"""high-precision direct residual r = V~ + k_xx - (1/4-a) p - b n at a few x (acb), vs the kernel value"""
import numpy as np
from flint import arb, acb, arb_poly, ctx
from base2 import Base; from arbops import SlabPolish, A; import resid
from rep import DS1
ctx.prec = 128
B = Base('P11'); i = 12; R, D = B.rep(i); q = B.slabs[i]
SP = SlabPolish(R, D, B.RD, q['nodes'], q['dt'])
Y, Yd = SP.Y_at(0.0); els, gl = SP.op.elements(Y, Yd)
a = gl['a']; b = gl['k2t'] + (arb(1) / 4 - a) * A(DS1)
NG = 60; gr = [arb.legendre_p_root(NG, k, weight=True) for k in range(NG)]
gu = [(g + 1) / 2 for g, w in gr]; gw = [w / 2 for g, w in gr]
th = arb(-1) * arb.pi() / 2; tht = arb(0); Gm = acb(-A(DS1)); Gt = acb(0); Kz = acb(0); rho_acc = arb(0)
sm = resid.samples(SP, 0.0); o, f = resid.norms(sm)
targets = [2.0, 5.0, 10.0, 20.0]; ti = 0
for el in els:
    he = el['h']; xe = el['x']
    thp = arb_poly([th]) + el['kp'].integral() * he            # theta on the element
    ttp = arb_poly([tht]) + el['kt'].integral() * he           # theta_tau on the element
    rkp = el['kt'] - el['Phi'].derivative() * (1 / he); rtp = arb_poly([rho_acc]) + rkp.integral() * he   # rho_theta
    while ti < len(targets) and float((xe + he).mid()) >= targets[ti]:
        us = (A(targets[ti]) - xe) / he
        # integrate from u=0 to us with a scaled Gauss rule
        Gx = Gm; Gtx = Gt; Kx = Kz
        for uq, wq in zip(gu, gw):
            uu = us * uq; T = acb.exp_pi_i(acb(thp(uu)) / arb.pi()) if False else acb(thp(uu).cos(), thp(uu).sin())
            Gx += he * us * wq * T; Gtx += he * us * wq * acb(0, 1) * ttp(uu) * T; Kx += he * us * wq * rtp(uu) * acb(thp(uu).cos(), -thp(uu).sin())
        tx = thp(us); nu = acb(0, 1) * acb(tx.cos(), tx.sin())
        dot = lambda z, w: z.real * w.real + z.imag * w.imag
        V = dot(Gtx, nu); p = dot(Gx, nu); n = nu.real
        kxx = el['ku'].derivative()(us) / he**2
        r = V + kxx - (arb(1) / 4 - a) * p - b * n
        rK = (acb(tx.cos(), tx.sin()) * Kx).real
        X = sm['X']; jj = np.unravel_index(np.argmin(np.abs(X - targets[ti])), X.shape)
        print(f"x={targets[ti]}: r_direct(arb) = {float(r.mid()):.4e}  r_kernel(arb) = {float(rK.mid()):.4e}  r_kernel(nearest Gauss pt x={X[jj]:.3f}) = {f['r'][jj]:.4e}")
        ti += 1
    # advance to the element end
    for uq, wq in zip(gu, gw):
        T = acb(thp(uq).cos(), thp(uq).sin()); Gm += he * wq * T; Gt += he * wq * acb(0, 1) * ttp(uq) * T; Kz += he * wq * rtp(uq) * acb(thp(uq).cos(), -thp(uq).sin())
    th = thp(arb(1)); tht = ttp(arb(1)); rho_acc = rtp(arb(1))
    if ti >= len(targets): break
