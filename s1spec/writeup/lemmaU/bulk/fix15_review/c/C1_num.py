#!/usr/bin/env python3
"""C1_num.py -- [num] float sanity check of the reparametrisation identity of C1_identification.tex (fix round 14).  Not a proof step.

Independent of the derivation in C1_sym.py: a closed symmetric curve with tangent angle theta = -pi/2 + a + th(a), th odd with even
harmonics only (so the curve is closed), is moved by one explicit Euler step of the purely normal velocity -k_ss N (spectral derivatives
in the curve parameter p), then re-parametrised by normalised arclength from the axis point p = 0, and its new th is compared with
  th + d_sigma * RHS,  RHS = -th'''' - th'' + N(th),  N = -(v^2 + v^3/3)' + (1 + v) P(w^2 - <w^2>)   (step 4, P1 equation),
d_sigma = (2 pi/L)^4 dt, P = zero-mean antiderivative.  Expected: relative defect O(dt) (Euler), checked for two dt.
Output C1_num.json.  Run: nice -n 19 python3 -B C1_num.py"""
import sys, os, json, hashlib, platform
sys.dont_write_bytecode = True
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
M = 1025                       # odd: no Nyquist mode (odd spectral derivatives stay real)
pgrid = 2 * np.pi * np.arange(M) / M
kw = np.fft.fftfreq(M, d=1.0 / M)


def dp(f, m=1):
    return np.real(np.fft.ifft((1j * kw) ** m * np.fft.fft(f)))


def antider_zero_mean(f):
    F = np.fft.fft(f); F[0] = 0; G = np.zeros_like(F); nz = kw != 0
    G[nz] = F[nz] / (1j * kw[nz]); return np.real(np.fft.ifft(G))


def th0(a):
    return 0.05 * np.sin(2 * a) + 0.02 * np.sin(4 * a) - 0.01 * np.sin(6 * a)


def rhs(th):
    v, w = dp(th), dp(th, 2)
    N = -dp(v ** 2 + v ** 3 / 3) + (1 + v) * antider_zero_mean(w ** 2 - np.mean(w ** 2))
    return -dp(th, 4) - dp(th, 2) + N


def curve_from_th(th, L):
    theta = -np.pi / 2 + pgrid + th
    T = np.exp(1j * theta)
    # X(p) = X0 + (L/2pi) int_0^p e^{i theta}; spectral antiderivative of the periodic part
    m = np.mean(T)
    assert abs(m) < 1e-13, ('not closed', m)
    X = (L / (2 * np.pi)) * (antider_zero_mean(T.real) + 1j * antider_zero_mean(T.imag))
    return X - X[0] + 1.0


def step_and_reparam(X, dt):
    Xp = dp(X.real) + 1j * dp(X.imag); r = np.abs(Xp); T = Xp / r; N = 1j * T
    theta = np.unwrap(np.angle(T))
    k = dp(theta - pgrid) / r + 1 / r            # theta - p is periodic
    ks = dp(k) / r; kss = dp(ks) / r
    Xn = X + dt * (-kss) * N
    Xnp = dp(Xn.real) + 1j * dp(Xn.imag); rn = np.abs(Xnp)
    Ln = np.mean(rn) * 2 * np.pi
    s = antider_zero_mean(rn - np.mean(rn)) + np.mean(rn) * pgrid
    s = s - s[0]                                   # arclength from p = 0
    alpha_of_p = 2 * np.pi * s / Ln
    thetan = np.unwrap(np.angle(Xnp))
    thn_p = thetan + np.pi / 2 - alpha_of_p        # th as a function of p (periodic)
    # evaluate th at the grid points a_j: solve alpha(p) = a_j (Newton on the spectral interpolant), then interpolate thn_p
    c_al = np.fft.fft(alpha_of_p - pgrid) / M; c_th = np.fft.fft(thn_p) / M

    def ev(c, z, der=0):
        return np.real(np.sum(c[None, :] * (1j * kw[None, :]) ** der * np.exp(1j * kw[None, :] * z[:, None]), axis=1))
    z = pgrid.copy()
    for _ in range(30):
        f = z + ev(c_al, z) - pgrid; fp = 1 + ev(c_al, z, 1); z = z - f / fp
    return ev(c_th, z), Ln


def run(dt, L=2 * np.pi):
    th = th0(pgrid); X = curve_from_th(th, L)
    thn, Ln = step_and_reparam(X, dt)
    dsig = (2 * np.pi / L) ** 4 * dt
    pred = th + dsig * rhs(th)
    # compare on the modes |n| <= 64 (the data have n <= 6; the modes ~ M/2 carry roundoff amplified by the spectral k_ss, O(1e-3) of the
    # rate per unit dt independent of dt, i.e. an aliasing artefact of the explicit step, not part of the identity)
    D = np.fft.fft(thn - pred); D[np.abs(kw) > 64] = 0
    defect = np.max(np.abs(np.real(np.fft.ifft(D)))) / (dsig * np.max(np.abs(rhs(th))))

    dlnL_num = (np.log(Ln) - np.log(L)) / dsig
    dlnL_exact = -np.mean(dp(th, 2) ** 2)          # -(1/2pi) int_0^{2pi} w^2
    return dict(dt=dt, rel_defect=float(defect), dlnL_num=float(dlnL_num), dlnL_formula=float(dlnL_exact),
                rel_defect_L=float(abs(dlnL_num - dlnL_exact) / abs(dlnL_exact)))


rows = [run(1e-6), run(1e-7)]
for r_ in rows: print(r_)
ok = rows[1]['rel_defect'] < 1e-4 and rows[1]['rel_defect'] < rows[0]['rel_defect'] and rows[1]['rel_defect_L'] < 1e-4
me = os.path.abspath(__file__)
out = dict(label='[num] float sanity check, not a proof step', ok=bool(ok), rows=rows,
           args=dict(argv=sys.argv, python=platform.python_version(), numpy=np.__version__, grid=M,
                     script=dict(path=me, sha256=hashlib.sha256(open(me, 'rb').read()).hexdigest()), inputs={}))
json.dump(out, open(os.path.join(HERE, 'C1_num.json'), 'w'), indent=1)
print('ok', ok)
