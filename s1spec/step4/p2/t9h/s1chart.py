"""Certified S1 chart for T9 (no new S1 certification: step-1/3 objects only).
x = sigma - sigma_tip (S1 arclength from the tip), tau = 8 - sigma = TT - x.
(i)  sigma <= 8: tight step-1 CAPD profile /home/claude/s1_tight/prof/prof_tight.txt (centre enclosures at t_c, hulls per box);
     third-order Taylor in tau from the box centre (tau-ODE p' = kq, q' = -1-kp, k' = -m, m' = -p/4, th' = -k), hull remainder.
(ii) sigma >= 8: envelope chart from the step-1 tail theorem bounds (step3 hstail.profile_bounds: |k|,|m|,|p| <= kb,mb,pb for
     sigma >= s1) and the asymptotic angle Theta, certified from the CAPD state at sigma = 8 by the exact identity
       int_s^oo k = p/s + 4m/s^2 + 8k/s^3 - 24p/s^5 - 480m/s^6 - 2880k/s^7 + 20160 int k/t^8 + 4 int N/t - 96 int N/t^5
     (k''' + (s/4)k = N, N = -k(k k' + J)), remainder <= 2880 kb/s^7 + N-terms.
Frame: Gamma1 = (q + i p) e^{i theta}, theta(tip) = -pi/2, Gamma1(tip) = (-2.5008.., 0) (the T1/T6 base convention)."""
import sys, math
import numpy as np
from flint import arb, acb, ctx
ctx.prec = 128
sys.path.insert(0, '/home/claude/s1spec/step3')
import hstail as HT
PROF = '/home/claude/s1_tight/prof/prof_tight.txt'
TT = arb((10.249077650575472 + 10.249077651519443) / 2, (10.249077651519443 - 10.249077650575472) / 2 * 1.000001)
TT = TT.union(arb(10.249077650575472)).union(arb(10.249077651519443))
import os as _os
if _os.environ.get('S1TT'):                     # conditioned tau_tip (s1_tight/ext_tip)
    _lo, _hi = map(float, _os.environ['S1TT'].split(',')); TT = arb(_lo).union(arb(_hi))
SIG_TIP = 8 - TT
I = acb(0, 1)

def iv(lo, hi): return arb(lo).union(arb(hi))

class S1:
    def __init__(self, path=PROF):
        rows = []
        for line in open(path):
            if line.startswith('END'): break
            rows.append(list(map(float, line.split())))
        self.A = np.array(rows); self.N = len(rows)
        self.Tend = self.A[-1, 1]; self.hb = self.A[0, 1] - self.A[0, 0]
        st = self.state(TT)                       # theta_prof at the tip
        self.th_off = -arb.pi() / 2 - st[4]
        st8 = self.state(arb(0))                  # sigma = 8
        p, q, k, m, th = st8
        self.st8 = st8; self.th8 = th + self.th_off
        self.G8 = acb(q, p) * (I * self.th8).exp()
        kb, mb, pb, qs = HT.profile_bounds(arb(8))
        Jb = 4 * kb * mb                          # crude: qs = 4 kb mb + 4 J (J <= qs/4)
        s = arb(8)
        lead = p / s + 4 * m / s ** 2 + 8 * k / s ** 3 - 24 * p / s ** 5 - 480 * m / s ** 6 - 2880 * k / s ** 7
        Nb = kb * (kb * mb + qs / 4)               # |N| <= kb (kb mb + J), J <= qs/4
        Aint = self.kint(arb(8))                  # int_8^oo kb
        rem = 2880 * kb / s ** 7 + (4 / s + 96 / s ** 5) * Nb * Aint / kb
        self.I8 = lead + arb(0, rem.upper())      # int_8^oo k
        self.Theta = self.th8 + self.I8
        self.rem8 = rem

    def kint(self, s0, L=60):
        """int_{s0}^oo kb(s) ds (kb decreasing, unit steps) + tail beyond s0 + L by the last value times a generous decay"""
        tot = arb(0); s = arb(s0)
        for n in range(L):
            tot += HT.profile_bounds(s + n)[0]
        return tot + HT.profile_bounds(s + L)[0] * 10

    def box(self, tau):
        i = int(min(max(float(tau.mid()) / self.hb, 0), self.N - 1))
        return i

    def state(self, tau):
        """(p, q, k, m, theta_prof) at tau (arb ball inside one box or at its edge)"""
        i = self.box(tau); r = self.A[i]
        tc = arb(r[2]); c = [iv(r[3 + 2 * j], r[4 + 2 * j]) for j in range(5)]; h = [iv(r[13 + 2 * j], r[14 + 2 * j]) for j in range(5)]
        D = tau - tc
        assert float(abs(D).upper()) <= self.hb * 1.01, (float(tau.mid()), i)
        def d1(p, q, k, m): return [k * q, -1 - k * p, -m, -p / 4, -k]
        def d2(p, q, k, m): return [-m * q - k - k * k * p, m * p - k * k * q, p / 4, -k * q / 4, m]
        def d3(p, q, k, m): return [p * q / 4 + 2 * m + 3 * k * m * p - k ** 3 * q, -p * p / 4 + 3 * k * m * q + k * k + k ** 3 * p, k * q / 4,
                                    (m * q + k + k * k * p) / 4, -p / 4]
        a1, a2 = d1(*c[:4]), d2(*c[:4]); a3 = d3(*h[:4])
        return [c[j] + a1[j] * D + a2[j] * D * D / 2 + a3[j] * D * D * D / 6 for j in range(5)]

    def tip_chart(self, x):
        """sigma <= 8 (x <= TT): Gamma1, theta1, k, k' = dk/dsigma, k'' (x arb)"""
        p, q, k, m, th = self.state(TT - x)
        th = th + self.th_off
        return dict(G=acb(q, p) * (I * th).exp(), th=th, k=k, k1=m, k2=p / 4)

    def Dset(self, a, b):
        """enclosure of theta1(t) - Theta for t in [a, b] (8 <= a <= b): intersect(theta8 - Theta +- kb(8)(b - 8), +- B(a))"""
        kb8 = HT.profile_bounds(arb(8))[0]
        A = kb8 * (arb(b) - 8)
        kb, mb, pb, qs = HT.profile_bounds(arb(a))
        a_ = arb(a)
        B = pb / a_ + 4 * mb / a_ ** 2 + 16 * kb / a_ ** 3            # |int_t^oo k| <= B(t) <= B(a) (envelopes decreasing)
        d1 = (self.th8 - self.Theta) + arb(0, A.upper()); d2 = arb(0, B.upper())
        lo = max(float(d1.lower()), float(d2.lower())); hi = min(float(d1.upper()), float(d2.upper()))
        assert lo <= hi, (a, b, lo, hi)
        return iv(lo, hi), kb, mb, pb

    def far_chart(self, x, n=None):
        """sigma >= 8: Gamma1 = G8 + e^{i Theta} int_8^sigma e^{i D(t)} dt (piecewise hulls), theta1 = Theta + D(sigma)"""
        s = x - TT + 8
        su = float(s.upper()); assert float(s.lower()) >= 8 - 1e-9
        n = n or max(4, int((su - 8) / 0.02))
        tot = acb(0); a = 8.0; step = (su - 8.0) / n
        for j in range(n):
            pa, pb_ = a + j * step, a + (j + 1) * step
            D, _, _, _ = self.Dset(pa, pb_)
            tot += arb(pb_ - pa) * acb(D.cos(), D.sin())
        tot -= acb(arb(0, (s.upper() - s.lower())), 0)          # s-ball: length uncertainty along the ray (crude)
        D, kb, mb, pb = self.Dset(float(s.lower()), su)
        G = self.G8 + (I * self.Theta).exp() * tot
        return dict(G=G, th=self.Theta + D, D=D, tot=tot, k=arb(0, kb.upper()), k1=arb(0, mb.upper()), k2=arb(0, (pb / 4).upper()))

if __name__ == '__main__':
    S = S1()
    print('theta_off', S.th_off.str(15), ' theta(8)', S.th8.str(15), ' int_8^oo k', S.I8.str(10), ' rem', S.rem8.str(3))
    print('THETA (THINF)', S.Theta.str(15))
    for x in (0.0, 1.0, 5.0, 10.0):
        c = S.tip_chart(arb(x)); print(x, c['th'].str(12), c['G'].real.str(12), c['G'].imag.str(12), c['k'].str(10))
    for x in (10.5, 12.0, 15.0, 30.0):
        c = S.far_chart(arb(x)); print(x, 'th', c['th'].str(8), 'G', c['G'].real.str(8), c['G'].imag.str(8), 'k<=', c['k'].str(3))
