"""Task 4: ramp geometry of Lemma U on [28, 56] for the weight Theta*phi (same formulas as u2_K.py; import-safe).
Set module globals G_LO, G_HI, KB, A_W before use (defaults: primary g-range, K_b = 1e-3, a_W = 0.3)."""
import sys, os, math
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import t4lib as L
sys.path.insert(0, L.P2)
from flint import arb
import ulib as U
from ulib import hull, LN2, LN28, LN56
from t4lib import A, Q_
G_LO, G_HI = Q_(300000), Q_(640000); A_W = 0.3; KB = A(Q_(1, 1000))
WT = 16.0
def ramp_q(box, point=False):
    tb = box[0]
    t = arb((tb[0] + tb[1])/2) if point else hull(tb[0], tb[1])
    E2 = [U.e2(t, l) for l in range(5)]
    (wlo, whi), (alo, ahi), (gl, gh) = box[1], box[2], box[3]
    a = arb((alo + ahi)/2) if point else hull(alo, ahi)
    if whi == math.inf:
        p, q = U.pq_lower_tail((arb(tb[1]) - 1)*LN2 - wlo if not point else (t - 1)*LN2 - wlo, a)
    else:
        w = arb((wlo + whi)/2) if point else hull(wlo, whi)
        p, q = U.pq_from_sigma((t - 1)*LN2 - w, a)
    Rk = U.Rk(p, q)
    gb = arb((gl + gh)/2) if point else hull(gl, gh)
    s = LN28 + t*LN2; th = U.theta_s(s, gb)
    RT = {0: arb(1), 1: Rk[1] - 4*th, 2: Rk[2] + th*(16 - 8*Rk[1]), 3: Rk[3] + th*(-12*Rk[2] + 48*Rk[1] - 64)}
    be = (KB/28)*(-t*LN2).exp()
    f1, f2, f3 = -be, be + be*be, -(be + 3*be*be + be*be*be)
    RP = {0: arb(1), 1: RT[1] + f1, 2: RT[2] + 2*RT[1]*f1 + f2, 3: RT[3] + 3*RT[2]*f1 + 3*RT[1]*f2 + f3}
    rem = {m: sum((math.comb(m, l)*E2[l]*RP[m - l] for l in range(1, m + 1)), arb(0)) for m in (1, 2, 3)}
    Lp = -4 + RP[1]
    rho = (E2[2] + 2*Lp*E2[1])/2
    et = [U.eta(t, l) for l in range(3)]
    ceta = Lp*et[0]*et[1] + et[1]*et[1] + et[0]*et[2]
    return dict(rem=rem, rho=rho, ceta=ceta, L=Lp, RP=RP, W1=Rk[1], th=th)
def ramp_boxes():
    ts = [(i/16, (i + 1)/16) for i in range(16)]
    ws = [(0.0, 0.25), (0.25, 0.5), (0.5, 1.0), (1.0, 2.0), (2.0, 4.0), (4.0, 8.0), (8.0, WT), (WT, math.inf)]
    return [(tt, ww, (A_W, A_W), (float(G_LO), float(G_HI))) for tt in ts for ww in ws]
def mkramp(coef=None, Y=None, which='K'):
    def f(box, point=False):
        q = ramp_q(box, point)
        if which == 'K':
            tot = arb(0)
            for m, c in coef.items(): tot += c*abs(q['rem'][m])
            if Y is not None: tot += Y*abs(q['rho'])
            return tot
        if which in ('rho', 'ceta', 'L', 'W1'): return q[which]
        if which.startswith('RP'): return q['RP'][int(which[2])]
        raise ValueError(which)
    f.scales = [0.6931, 1.0, 1.0, 1e-5]
    return f
