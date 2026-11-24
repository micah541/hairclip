"""Certified S1 chart v2 (step-1 objects only): prof_tight (sigma in [tip, 8]) + ext813 (sigma in [8, 13], /home/claude/s1_tight/ext813)
+ tail envelopes beyond 13.  tau = 8 - sigma; x = sigma - sigma_tip = TT - tau, TT = T0 + dT (dT = the tip-position ball).
Correlation in TT is kept: every quantity is evaluated at tau = T0 - x (x exact) and moved by dT * (its tau-derivative); the
angle normalisation theta_abs(x) = theta(T0 - x) - theta(T0) - pi/2 + dT (k_tip - k(x)).
Theta (asymptotic angle) = theta_abs(13) + int_13^oo k, by the identity (remainder 2880 kb(13)/13^7 + N-terms)."""
import sys, math, os
import numpy as np
from flint import arb, acb, ctx
ctx.prec = 128
sys.path.insert(0, '/home/claude/s1spec/step3'); sys.path.insert(0, '/home/claude/s1spec/step4/p2/t9h')
import hstail as HT
from s1chart import TT
def iv(a): return arb(a[0]).union(arb(a[1]))
P1 = os.environ.get('S1PROF', '/home/claude/s1_tight/prof/prof_tight.txt'); P2 = '/home/claude/s1_tight/ext813/ext813_prof.txt'
I = acb(0, 1)
T0 = arb(TT.mid()); dT = arb(0, TT.rad())

def load(path):
    rows = []
    for line in open(path):
        if line.startswith('END'): break
        rows.append(list(map(float, line.split())))
    return np.array(rows)

def d1(p, q, k, m): return [k * q, -1 - k * p, -m, -p / 4, -k]                         # d/dtau
def d2(p, q, k, m): return [-m * q - k - k * k * p, m * p - k * k * q, p / 4, -k * q / 4, m]
def d3(p, q, k, m): return [p * q / 4 + 2 * m + 3 * k * m * p - k * k * k * q, -p * p / 4 + 3 * k * m * q + k * k + k * k * k * p, k * q / 4,
                            (m * q + k + k * k * p) / 4, -p / 4]

class S1:
    def __init__(self):
        self.A1, self.A2 = load(P1), load(P2)
        self.hb = 2e-4
        self.hstep = 0.99 * min(self.hb, (self.A1[-1, 1] - self.A1[0, 0]) / len(self.A1))   # never skips a row
        st = self.state(T0); self.thP0, self.ktip = st[4], st[2]
        Cfull = -arb.pi() / 2 - self.thP0 + dT * self.ktip          # theta_abs = theta_prof(T0 - x) - dT k(x) + C
        import os
        if os.environ.get('COHERENT') == '1':                       # G1 fix (a): the tip offset C as one scalar, evaluated at its midpoint
            self.C, self.C_rad = arb(Cfull.mid()), float(Cfull.rad())
        else:
            self.C, self.C_rad = Cfull, 0.0
        s13 = self.state(arb(-5)); p, q, k, m, th = s13
        self.th13 = th + self.C
        self.G13 = acb(q, p) * (I * self.th13).exp()
        kb, mb, pb, qs = HT.profile_bounds(arb(13)); s = arb(13)
        lead = p / s + 4 * m / s ** 2 + 8 * k / s ** 3 - 24 * p / s ** 5 - 480 * m / s ** 6 - 2880 * k / s ** 7
        Nb = kb * (kb * mb + qs / 4); Aint = sum((HT.profile_bounds(s + n)[0] for n in range(60)), arb(0)) * 2
        rem = 2880 * kb / s ** 7 + (4 / s + 96 / s ** 5) * Nb * Aint / kb
        self.I13 = lead + arb(0, rem.upper())
        self.Theta = self.th13 + self.I13
        self.kb13 = kb

    def row(self, tau):
        t = float(tau.mid())
        if t >= 0:
            A = self.A1; h = (A[-1, 1] - A[0, 0]) / len(A); i = int(min(max(t / h, 0), len(A) - 1))
            while i > 0 and t < A[i, 0]: i -= 1
            while i < len(A) - 1 and t > A[i, 1]: i += 1
        else:
            A = self.A2; i = int(min(max(-t / self.hb, 0), len(A) - 1))
            while i > 0 and t > A[i, 1]: i -= 1
            while i < len(A) - 1 and t < A[i, 0]: i += 1
        return A[i]

    def state(self, tau):
        r = self.row(tau)
        tc = arb(r[2]); c = [iv([r[3 + 2 * j], r[4 + 2 * j]]) for j in range(5)]; h = [iv([r[13 + 2 * j], r[14 + 2 * j]]) for j in range(5)]
        D = tau - tc
        assert float(abs(D).upper()) <= self.hb * 0.51 + 1e-12, (float(tau.mid()), r[:3])
        a1, a2, a3 = d1(*c[:4]), d2(*c[:4]), d3(*h[:4])
        return [c[j] + a1[j] * D + a2[j] * D * D / 2 + a3[j] * D * D * D / 6 for j in range(5)]

    def hull(self, tau_lo, tau_hi):
        """state hulls over tau in [tau_lo, tau_hi] (union of row hulls)"""
        out = None; t = tau_lo
        while t <= tau_hi + 1e-15:
            r = self.row(arb(t)); h = [iv([r[13 + 2 * j], r[14 + 2 * j]]) for j in range(5)]
            out = h if out is None else [a.union(b) for a, b in zip(out, h)]
            t += self.hstep
        r = self.row(arb(tau_hi)); h = [iv([r[13 + 2 * j], r[14 + 2 * j]]) for j in range(5)]
        return [a.union(b) for a, b in zip(out, h)]

    def chart(self, x):
        """sigma <= 13 (x <= T0 + 5): Gamma1, theta1, k, k', k'', k''' (sigma-derivatives; x exact arb point)"""
        st = self.state(T0 - x); p, q, k, m, th = st
        g = d1(p, q, k, m)
        p, q = p + dT * g[0], q + dT * g[1]; k2, m2 = k + dT * g[2], m + dT * g[3]
        th_abs = th - dT * k + self.C
        return dict(G=acb(q, p) * (I * th_abs).exp(), qp=acb(q, p), th=th_abs, k=k2, k1=m2, k2=p / 4, k3=-k2 * q / 4)

    def chart_ball(self, x):
        """as chart(x) for an x-ball (sigma <= 13): the state hulls over tau in T0 - x (union of row hulls), then the same corrections"""
        tb = T0 - x
        p, q, k, m, th = self.hull(float(tb.lower()) - 1e-12, float(tb.upper()) + 1e-12)[:5]
        g = d1(p, q, k, m)
        p, q = p + dT * g[0], q + dT * g[1]; k2, m2 = k + dT * g[2], m + dT * g[3]
        th_abs = th - dT * k + self.C
        return dict(G=acb(q, p) * (I * th_abs).exp(), qp=acb(q, p), th=th_abs, k=k2, k1=m2, k2=p / 4, k3=-k2 * q / 4)

    def Dset(self, a, b):
        """theta1(t) - Theta for sigma = t in [a, b] (13 <= a <= b)"""
        A = self.kb13 * (arb(b) - 13)
        kb, mb, pb, qs = HT.profile_bounds(arb(a)); a_ = arb(a)
        B = pb / a_ + 4 * mb / a_ ** 2 + 16 * kb / a_ ** 3
        d1_ = (self.th13 - self.Theta) + arb(0, A.upper()); d2_ = arb(0, B.upper())
        lo = max(float(d1_.lower()), float(d2_.lower())); hi = min(float(d1_.upper()), float(d2_.upper()))
        assert lo <= hi
        return iv([lo, hi]), kb, mb, pb

    def far_chart(self, x, step=0.05):
        """sigma >= 13: Gamma1 = G13 + e^{i Theta} int_13^sigma e^{i D}, theta1 = Theta + D"""
        s = x - T0 + 8                                            # sigma (dT moves the point along the ray only)
        su = float(s.upper()); n = max(2, int((su - 13.0) / step) + 1); stp = (su - 13.0) / n
        tot = acb(0)
        for j in range(n):
            D, _, _, _ = self.Dset(13.0 + j * stp, 13.0 + (j + 1) * stp)
            tot += arb(stp) * acb(D.cos(), D.sin())
        tot += acb(arb(0, (dT.rad() * 2 + s.rad() * 2)), 0)
        D, kb, mb, pb = self.Dset(max(float(s.lower()), 13.0), su)
        return dict(G=self.G13 + (I * self.Theta).exp() * tot, th=self.Theta + D, D=D, tot=tot, k=arb(0, kb.upper()), k1=arb(0, mb.upper()),
                    k2=arb(0, (pb / 4).upper()), k3=arb(0, (kb * (s + 1) / 4).upper()))

if __name__ == '__main__':
    S = S1()
    print('Theta', S.Theta.mid(), 'rad', S.Theta.rad(), ' th13', S.th13.rad(), ' I13', S.I13.mid(), S.I13.rad())
    for x in (0.0, 1.0, 5.0, 10.0, 12.0, 15.0):
        c = S.chart(arb(x)); print(x, 'th rad', float(c['th'].rad()), 'G rad', float(c['G'].real.rad()), float(c['G'].imag.rad()), 'k rad', float(c['k'].rad()))
    for x in (16.0, 30.0, 60.0):
        c = S.far_chart(arb(x)); print(x, 'far th rad', float(c['th'].rad()), 'G rad', float(c['G'].real.rad()), float(c['G'].imag.rad()))
