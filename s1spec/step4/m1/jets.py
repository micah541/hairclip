"""Stage 3: Lemma O jets (cone cylinders X=3, half-width 0.75 r, L=0.25 r; far cylinders X=2, L=113, unit lambda_1) from the tangent angle theta (spectral)."""
import sys, numpy as np
from geo import *
THINF_ = THINF

def curve_arrays(path, lam, xc):
    S, z = CH.H.load(path)
    N = S.N; L = S.L
    al, th, Z, _ = S.curve(1)
    thp = th - al                                      # periodic part (omega = 1)
    th_hat = S.th_hat
    d1_ = S.d(th_hat, 1); d2_ = S.d(th_hat, 2); d3_ = S.d(th_hat, 3)
    c = 2 * np.pi / L
    idx = (N // 2 - np.arange(N)) % N                  # j=0 tip, j increasing = alpha decreasing (upper strand), wraps through alpha=0 into lower strand
    alj = np.pi - np.arange(N) * (2 * np.pi / N)
    phi = alj + thp[idx] + np.pi                       # outward tangent angle, continuous
    D1 = -lam * c * (1 + d1_[idx]); D2 = lam**2 * c**2 * d2_[idx]; D3 = -lam**3 * c**3 * d3_[idx]
    G = (Z[idx] - xc) / lam
    ds = (L / N) / lam
    return dict(phi=phi, D1=D1, D2=D2, D3=D3, G=G, ds=ds, N=N)

def cyl_jets(A, jc, Lc, X, nwin=None):
    """graph jets of the curve over its tangent line at sample jc, within |x| <= X Lc (units lambda_1). returns dict"""
    N = A['N']; ds = A['ds']
    w = int(1.25 * X * Lc / ds) + 2
    sl = (jc + np.arange(-w, w + 1)) % N
    ph = A['phi'][sl]; ph = ph - ph[w]              # phi is continuous mod wrap: fix 2pi jumps
    ph = (ph + np.pi) % (2 * np.pi) - np.pi
    k0, k1, k2 = A['D1'][sl], A['D2'][sl], A['D3'][sl]
    cs = np.cos(ph)
    xx = np.concatenate([[0], np.cumsum(0.5 * (cs[1:] + cs[:-1]) * ds)]); xx -= xx[w]
    m = np.abs(xx) <= X * Lc
    sec = 1 / cs; t = np.tan(ph)
    h2 = k0 * sec**3
    h3 = sec**4 * (k1 + 3 * k0**2 * t)
    h4 = sec * (4 * sec**4 * t * k0 * (k1 + 3 * k0**2 * t) + sec**4 * (k2 + 6 * k0 * k1 * t + 3 * k0**3 * sec**2))
    x = xx[m]
    n3 = float(np.sqrt(abs(np.trapezoid(h3[m]**2, x))))
    return dict(hx=float(np.max(np.abs(t[m]))), h2=float(np.max(np.abs(h2[m]))), h3=float(np.max(np.abs(h3[m]))), h4=float(np.max(np.abs(h4[m]))), n3=n3,
                h2c=float(h2[w]), h3c=float(h3[w]), h4c=float(h4[w]))

def run(path, lam, xc, cone_ranges=((200, 226), (200, 400), (200, 700)), far_r0=376.0, stride_cone=0.01):
    A = curve_arrays(path, lam, xc)
    N = A['N']; r = np.abs(A['G'])
    # centres on the upper strand j in [0, N/2]; r monotone increasing initially
    out = {'cone': {}, 'far': {}}
    jj = np.arange(N // 2 + 1)
    # first index with r >= value
    def j_of(rv): return int(np.argmax(r[:N // 2 + 1] >= rv))
    # cone family (X=3, L=0.25 r_c)
    cone_rows = []
    rv = 200.0; rmax = max(c[1] for c in cone_ranges) * 1.0
    while rv <= rmax:
        jc = j_of(rv); Lc = 0.25 * r[jc]
        d = cyl_jets(A, jc, Lc, 3.0)
        cone_rows.append((r[jc], Lc, d['hx'], Lc * d['h2'], Lc**2 * d['h3'], Lc**3 * d['h4'], Lc**1.5 * d['n3']))
        rv *= 1 + stride_cone
    out['cone_rows'] = np.array(cone_rows)    # r_c, L, slope, h2, h3, h4, N3 (cylinder units)
    # far family: X=2, L=113; centres r_c>=far_r0 along the whole strand (upper), jets in lambda_1 units
    Lf = 113.0
    far = []
    j = j_of(far_r0)
    while j < N // 2 - 5:
        d = cyl_jets(A, j, Lf, 2.0)
        far.append((r[j], d['hx'], d['h2'], d['h3'], d['h4'], Lf**1.5 * d['n3']))
        j += 20 if r[j] < 2000 else 200
    out['far_rows'] = np.array(far)
    # pointwise centre jets (kappa, kappa_s, kappa_ss in lambda_1 units) for r >= 226 along the whole closed curve (both strands)
    rr = r
    m = rr >= 226.0
    out['pt'] = np.array([np.max(np.abs(A['D1'][m])), np.max(np.abs(A['D2'][m])), np.max(np.abs(A['D3'][m]))])
    j0 = int(np.argmax(np.abs(A['D1'][m]))); out['pt_at'] = float(rr[m][j0])
    out['pt_r'] = [float(rr[m][int(np.argmax(np.abs(A[k][m])))]) for k in ('D1', 'D2', 'D3')]
    return out

if __name__ == '__main__':
    d = dict(np.load(sys.argv[1])); path = sys.argv[2]
    o = run(path, float(d['lam']), complex(d['xc']))
    cr = o['cone_rows']; fr = o['far_rows']
    np.set_printoptions(linewidth=200, precision=4)
    print("cone rows (r_c, L, slope, h2, h3, h4, N3):"); print(cr[::max(1, len(cr)//12)])
    print("far rows max (slope,h2,h3,h4,N3):", fr[:, 1:].max(0), " r of max:", [fr[np.argmax(fr[:, i]), 0] for i in range(1, 6)])
    print("pointwise kappa, kappa_s, kappa_ss max for r>=226:", o['pt'], o['pt_r'])
    print("thresholds far: slope .21 h2 8.9e-3 h3 7.0e-5 h4 6.2e-7 N3 1.0; cone: 0.12 each, N3 0.14")
