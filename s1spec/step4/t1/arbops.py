"""T1 arb (python-flint) layer: exact polynomial evaluation of the slab collocation residual, polish of one slab
(float Jacobian inverse, arb residual), and the pointwise PDE residuals rho_kappa, rho_theta at any tau of the slab.
Per element everything is a polynomial in u (kappa, Q = int kappa_x^2, Phi = theta_tau flux, ALE), so all test
integrals are exact (up to arb rounding); only the knot map x_j(ell) = XT (1 + C ell^3)^(1/3) is irrational."""
import numpy as np
from fractions import Fraction
from math import comb, factorial
from flint import arb, arb_poly, arb_mat, fmpq, ctx
import scipy.sparse.linalg as spla
from rep import K0, DS1

def A(v):
    """exact arb from a float / int / Fraction"""
    if isinstance(v, Fraction): return arb(fmpq(v.numerator, v.denominator))
    return arb(float(v)) if not isinstance(v, arb) else v

def fr_poly_mul(a, b):
    out = [Fraction(0)] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b): out[i + j] += x * y
    return out

def hermite_exact(P):
    """monomial coefficients (Fractions) of H0m, H1m (m=0..4) and bubbles u^5(1-u)^5 Leg_i(2u-1)"""
    cols = []
    om5 = [Fraction(comb(5, k) * (-1) ** k) for k in range(6)]
    for m in range(5):
        s = [Fraction(comb(4 + k, k)) for k in range(5 - m)]
        c = fr_poly_mul([Fraction(0)] * m + [Fraction(1, factorial(m))], fr_poly_mul(om5, s))
        cols.append(c + [Fraction(0)] * (P + 1 - len(c)))
    for m in range(5):
        c0 = cols[m]; c = [Fraction(0)] * (P + 1)
        for k, ck in enumerate(c0):
            if ck == 0: continue
            for i in range(k + 1): c[i] += ck * comb(k, i) * (-1) ** i
        cols.append([(-1) ** m * v for v in c])
    b5 = fr_poly_mul([Fraction(0)] * 5 + [Fraction(1)], om5)
    for i in range(P - 9):
        li = legendre_shift(i)
        bc = fr_poly_mul(b5, li); cols.append(bc + [Fraction(0)] * (P + 1 - len(bc)))
    return cols

def legendre_shift(m):
    """monomial coefficients (Fractions) of Leg_m(2u - 1)"""
    # P_m(2u-1) = sum_k (-1)^(m+k) C(m,k) C(m+k,k) u^k
    return [Fraction((-1) ** (m + k) * comb(m, k) * comb(m + k, k)) for k in range(m + 1)]


class ArbBasis:
    def __init__(self, FB, prec=128):
        """FB: float Basis (for nrm and sizes)"""
        ctx.prec = prec
        self.P = FB.P; self.nb = FB.nb; self.nbub = FB.nbub; self.nm = FB.nmom - 1
        cols = hermite_exact(self.P)
        self.nrm = [Fraction(float(v)) for v in FB.nrm]
        self.Nq = [arb_poly([A(c / self.nrm[j]) for c in cols[j]]) for j in range(self.nb)]      # normalised local basis
        self.Cmat = arb_mat(self.P + 1, self.nb, [A(cols[j][i] / self.nrm[j]) for i in range(self.P + 1) for j in range(self.nb)])
        K = 3 * self.P + 4
        # integral functionals on monomials u^k, k = 0..K:  int_0^1 u^k g(u) du = sum_i g_i / (k + i + 1)
        def W(gs):
            return arb_mat(K + 1, len(gs), [A(sum((gi / (k + i + 1) for i, gi in enumerate(g)), Fraction(0))) for k in range(K + 1) for g in gs])
        Pm = [legendre_shift(m) for m in range(self.nm)]
        dPm = [[i * c for i, c in enumerate(p)][1:] or [Fraction(0)] for p in Pm]
        N9 = [c / self.nrm[9] for c in cols[9]]; N4 = [c / self.nrm[4] for c in cols[4]]
        dN9 = [i * c for i, c in enumerate(N9)][1:]; dN4 = [i * c for i, c in enumerate(N4)][1:]
        self.WPhi = W(dPm + [dN9, dN4])         # Phi against d/du of test functions
        self.WV = W(Pm + [N9, N4])               # kappa_tau|_x against test functions
        self.K = K
        self.mfac = [2 * m + 1 for m in range(self.nm)]

    def apply_W(self, poly, Wm):
        c = poly.coeffs(); c = c + [arb(0)] * (self.K + 1 - len(c))
        row = arb_mat(1, self.K + 1, c[:self.K + 1])
        return (row * Wm).entries()


class ArbOp:
    """implicit collocation residual of one mesh (float Rep R) in arb"""
    def __init__(self, R, AB):
        self.R, self.AB = R, AB; self.n = R.n
        m = R.mesh
        self.xtip = [A(v) for v in m.xtip]; self.C = [A(v) for v in m.C]; self.XT = A(m.XT); self.JT = m.JT
        self.L = R.L; self.fixval = R.fixval; self.N = R.N
        self.i02, self.i04 = R.nidx[0, 2], R.nidx[0, 4]

    def knots(self, ell):
        third = arb(1) / 3
        xo, dxo = [], []
        for Cj in self.C:
            g = 1 + Cj * ell**3; gp = g ** third
            xo.append(self.XT * gp); dxo.append(self.XT * Cj * ell**2 * gp / g)
        x = self.xtip + xo + [ell]; dx = [arb(0)] * (self.JT + 1) + dxo + [arb(1)]
        return x, dx

    def geometry(self, ell):
        x, dx = self.knots(ell); N = self.N
        h = [x[e + 1] - x[e] for e in range(N)]; dh = [dx[e + 1] - dx[e] for e in range(N)]
        H = [h[0]] + [(h[j - 1] * h[j]).sqrt() for j in range(1, N)] + [h[N - 1]]
        g = [dh[e] / h[e] for e in range(N)]
        dHH = [g[0]] + [(g[j - 1] + g[j]) / 2 for j in range(1, N)] + [g[N - 1]]      # H'/H
        return x, dx, h, dh, H, dHH

    def local(self, y, h, H, e):
        """local coefficients of element e (arb list); y: list of arb (global unknowns)"""
        AB = self.AB; d = []
        for q in range(AB.nb):
            g = self.L[e, q]
            v = y[g] if g >= 0 else A(self.fixval[e, q])
            if q < 5: v = v * (h[e] / H[e]) ** q
            elif q < 10: v = v * (h[e] / H[e + 1]) ** (q - 5)
            d.append(v)
        return d

    def poly(self, d):
        col = self.AB.Cmat * arb_mat(len(d), 1, d)
        return arb_poly(col.entries())

    def elements(self, Y, Yd):
        """per-element polynomials for state Y and time derivative Yd (lists of arb, length n+4): yields dicts"""
        AB, n, N = self.AB, self.n, self.N
        ell = Y[n]; elld = Yd[n]
        x, dx, h, dh, H, dHH = self.geometry(ell)
        h0 = h[0]
        k2t = Y[self.i02] / (h0**2 * A(AB.nrm[2])); k4t = Y[self.i04] / (h0**4 * A(AB.nrm[4]))
        a = arb(1) / 4 + (k4t + A(K0)**2 * k2t) / A(K0)
        qa = arb(1) / 4 - a
        Qe = arb(0); out = []
        u = arb_poly([0, 1])
        for e in range(N):
            d = self.local(Y, h, H, e); dd = self.local_lin(Yd, h, H, e)
            kp = self.poly(d); ku = kp.derivative(); kuu = ku.derivative(); kuuu = kuu.derivative()
            he = h[e]; ih = 1 / he
            Q = arb_poly([Qe]) + (ku * ku).integral() * ih
            Phi = -kuuu * (ih**3) - kp * kp * ku * ih + kp * Q - (arb_poly([x[e]]) + u * he) * kp * qa
            # ALE generator: nodal rescaling + mesh motion
            gr = [arb(0)] * AB.nb
            for q in range(5):
                gr[q] = d[q] * q * (dh[e] / he - dHH[e]); gr[5 + q] = d[5 + q] * q * (dh[e] / he - dHH[e + 1])
            ALE = self.poly(gr) - ku * ih * (arb_poly([dx[e]]) + u * dh[e])
            kt = self.poly(dd) + elld * ALE              # kappa_tau at fixed x
            out.append(dict(e=e, x=x[e], h=he, Phi=Phi, kt=kt, kp=kp, ku=ku, Q=Q, ALE=ALE))
            Qe = Q(arb(1))
        return out, dict(a=a, k2t=k2t, Qend=Qe, ell=ell, x=x, h=h)

    def local_lin(self, yd, h, H, e):
        AB = self.AB; d = []
        for q in range(AB.nb):
            g = self.L[e, q]
            v = yd[g] if g >= 0 else arb(0)
            if q < 5: v = v * (h[e] / H[e]) ** q
            elif q < 10: v = v * (h[e] / H[e + 1]) ** (q - 5)
            d.append(v)
        return d

    def residual(self, Y, Yd):
        """implicit residual, same row layout as the float mass()/rhs() (moment rows, Galerkin rows), plus 4 scalars"""
        AB, n, N, nm = self.AB, self.n, self.N, self.AB.nm
        els, gl = self.elements(Y, Yd)
        h = gl['h']; hn = [(h[j - 1] + h[j]) / 2 for j in range(1, N)] + [h[N - 1]]
        Hn = None
        E = [arb(0)] * n
        _, _, hh, _, H, _ = self.geometry(Y[n])
        gal = [arb(0)] * (N + 1)
        for el in els:
            e = el['e']; he = el['h']; Phi = el['Phi']
            P1, P0 = Phi(arb(1)), Phi(arb(0))
            wP = AB.apply_W(Phi, AB.WPhi); wV = AB.apply_W(el['kt'], AB.WV)
            for m in range(nm):
                GM = AB.mfac[m] * (P1 - (-1) ** m * P0 - wP[m]) / he
                E[e * (nm + 1) + m] = AB.mfac[m] * wV[m] - GM
            sc9 = (he / H[e + 1]) ** 4; sc4 = (he / H[e]) ** 4
            gal[e + 1] += sc9 * (he * wV[nm] + wP[nm])
            gal[e] += sc4 * (he * wV[nm + 1] + wP[nm + 1])
        for j in range(1, N + 1):
            E[(j - 1) * (nm + 1) + nm] = gal[j] / hn[j - 1]
        a = gl['a']; lam = Y[n + 1].exp()
        b = gl['k2t'] + (arb(1) / 4 - a) * A(DS1)
        Es = [Yd[n] - (-gl['Qend'] + (arb(1) / 4 - a) * Y[n]), Yd[n + 1] - (a - arb(1) / 4), Yd[n + 2] + b * lam, Yd[n + 3] - lam**4]
        return E, Es, gl


def mpf_to_arb(v): return arb(str(v))


class SlabPolish:
    """polish the stages of one Radau slab: Y(tau0 + theta dt) = sum_j L_j(theta) Y_j, nodes (0, c)"""
    def __init__(self, R, D, RD, nodes, dt, prec=128):
        ctx.prec = prec
        self.R, self.D, self.RD, self.dt = R, D, RD, dt
        self.AB = ArbBasis(R.B, prec); self.op = ArbOp(R, self.AB)
        s = RD.s; self.s = s
        import mpmath as mp
        mp.mp.dps = 60
        nd = [mp.mpf(0)] + list(RD.cm)
        # derivative matrix of the Lagrange basis on nd at nd (rows 1..s)
        def dL(j, t):
            tot = mp.mpf(0)
            for l in range(s + 1):
                if l == j: continue
                term = 1 / (nd[j] - nd[l])
                for k in range(s + 1):
                    if k != j and k != l: term *= (t - nd[k]) / (nd[j] - nd[k])
                tot += term
            return tot
        self.Dm = [[mpf_to_arb(dL(j, nd[i])) for j in range(s + 1)] for i in range(s + 1)]
        self.nd = nd
        self.Yn = [[A(v) for v in row] for row in nodes]           # (s+1) x (n+4) arb, row 0 fixed
        self.dtA = A(dt)

    def stage_derivs(self):
        s = self.s; ntot = len(self.Yn[0])
        return [[sum((self.Dm[i][j] * self.Yn[j][k] for j in range(s + 1)), arb(0)) / self.dtA for k in range(ntot)] for i in range(1, s + 1)]

    def residuals(self):
        Yd = self.stage_derivs(); out = []
        for i in range(1, self.s + 1):
            E, Es, gl = self.op.residual(self.Yn[i], Yd[i - 1]); out.append((E, Es))
        return out

    def step(self, J, res=None):
        """one simplified-Newton correction with the frozen float Jacobian J (factorised in self.RD)"""
        R, D, RD, s = self.R, self.D, self.RD, self.s; n = R.n
        res = self.residuals() if res is None else res
        Rz = np.zeros((s, n + 4)); emax = 0.0
        for i, (E, Es) in enumerate(res):
            Ef = np.array([float(v.mid()) for v in E]); Esf = np.array([float(v.mid()) for v in Es])
            emax = max(emax, np.abs(Ef).max(), np.abs(Esf).max())
            ell = float(self.Yn[i + 1][n].mid())
            Mm = D.mass(ell)
            Rz[i, :n] = -spla.spsolve(Mm, Ef); Rz[i, n:] = -Esf
        Rw = RD.Ti @ Rz
        dW = np.empty((s, n + 4), complex)
        import scipy.linalg as sla
        for k in RD.rep: dW[k] = sla.lu_solve(RD.lus[k], Rw[k])
        for k, j in RD.partner.items(): dW[k] = np.conj(dW[j])
        dZ = (RD.T @ dW).real
        for i in range(s):
            row = self.Yn[i + 1]
            for k in range(n + 4): row[k] = row[k] + arb(float(dZ[i, k]))
        return emax, np.abs(dZ).max()

    def Y_at(self, theta):
        """arb state and time derivative at relative time theta (float in [0,1])"""
        import mpmath as mp
        th = mp.mpf(theta); nd = self.nd; s = self.s
        Lw, dLw = [], []
        for j in range(s + 1):
            v = mp.mpf(1)
            for k in range(s + 1):
                if k != j: v *= (th - nd[k]) / (nd[j] - nd[k])
            Lw.append(mpf_to_arb(v))
            tot = mp.mpf(0)
            for l in range(s + 1):
                if l == j: continue
                term = 1 / (nd[j] - nd[l])
                for k in range(s + 1):
                    if k != j and k != l: term *= (th - nd[k]) / (nd[j] - nd[k])
                tot += term
            dLw.append(mpf_to_arb(tot))
        ntot = len(self.Yn[0])
        Y = [sum((Lw[j] * self.Yn[j][k] for j in range(s + 1)), arb(0)) for k in range(ntot)]
        Yd = [sum((dLw[j] * self.Yn[j][k] for j in range(s + 1)), arb(0)) / self.dtA for k in range(ntot)]
        return Y, Yd


def arb_refine(Y, R0, R1, AB):
    """exact (arb) transfer of a state Y (list, length n0+4) from mesh R0 to the refined mesh R1 (one knot more)"""
    op0, op1 = ArbOp(R0, AB), ArbOp(R1, AB); n0 = R0.n
    ell = Y[n0]
    x0, _, h0, _, H0, _ = op0.geometry(ell); x1, _, h1, _, H1, _ = op1.geometry(ell)
    C0, C1 = list(R0.mesh.C), list(R1.mesh.C)
    jn = R0.mesh.JT + 1 + next(i for i in range(len(C1)) if i >= len(C0) or C1[i] != C0[i])   # new knot index in mesh 1
    e = jn - 1                                                                                 # old element containing it
    # physical nodal derivatives of all old nodes from the element polynomials
    polys = [op0.poly(op0.local(Y, h0, H0, k)) for k in range(R0.N)]
    def jets(pk, uu, hk):
        out = []; p = pk
        for m in range(5):
            out.append(p(uu) / hk**m); p = p.derivative()
        return out
    D = [jets(polys[k], arb(0), h0[k]) for k in range(R0.N)] + [jets(polys[-1], arb(1), h0[-1])]
    us = (x1[jn] - x0[e]) / h0[e]
    Dn = jets(polys[e], us, h0[e])
    D1 = D[:e + 1] + [Dn] + D[e + 1:]
    # bubbles: old elements keep theirs; the two new sub-elements from the old polynomial
    from flint import arb_mat as AM
    Cinv = AB.Cmat.inv() if not hasattr(AB, 'Cinv') else AB.Cinv; AB.Cinv = Cinv
    bub = {}
    for k_new, (a0, a1) in ((e, (arb(0), us)), (e + 1, (us, arb(1)))):
        sub = polys[e](arb_poly([a0, a1 - a0]))                  # p(a0 + (a1 - a0) v)
        c = sub.coeffs(); c = c + [arb(0)] * (AB.P + 1 - len(c))
        dl = (Cinv * AM(AB.P + 1, 1, c[:AB.P + 1])).entries()
        bub[k_new] = dl[10:]
    y1 = [arb(0)] * (R1.n + 4)
    for j in range(R1.N + 1):
        for m in range(5):
            g = R1.nidx[j, m]
            if g >= 0: y1[g] = D1[j][m] * H1[j]**m * A(AB.nrm[m])
    for k in range(R1.N):
        if AB.nbub == 0: break
        if k in bub: vals = bub[k]
        else:
            kold = k if k < e else k - 1
            vals = [Y[g] for g in R0.bidx[kold]]
        for q, g in enumerate(R1.bidx[k]): y1[g] = vals[q]
    y1[R1.n:] = Y[n0:]
    return y1
