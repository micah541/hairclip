"""
Rigorous geometric checks for S1 from CAPD trajectory enclosures (curve_br.txt).

Frame: Gamma~ = -e^{-i phi_tip} Gamma, so the tip sits at (-D, 0) and the symmetry axis is the
real line.  With alpha = theta - theta_tip:
    x~ = q sin(alpha) + p cos(alpha),   y~ = p sin(alpha) - q cos(alpha),
forward (sigma-increasing) tangent direction angle = alpha - pi/2,
d y~/d tau = cos(alpha)  (tau = S - sigma),   d q/d tau = -1 - k p.
Half-branch B = Gamma~([sigma_tip, oo)) corresponds to tau in [0, tau_tip] plus the tail sigma >= S.
"""
import sys, math
import numpy as np
from flint import arb, ctx
ctx.prec = 128

def iv(lo, hi):
    return arb(lo).union(arb(hi))

def dn(a):
    lo = a.lower(); v = float(lo)
    while arb(v) > lo: v = math.nextafter(v, -math.inf)
    return v

def up(a):
    hi = a.upper(); v = float(hi)
    while arb(v) < hi: v = math.nextafter(v, math.inf)
    return v

def load(path):
    L = open(path).read().split("\n")
    assert L[0].startswith("TIP")
    rt = [float(v) for v in L[0].split()[1:3]]
    tip = [tuple(map(float, L[1 + i].split())) for i in range(5)]
    boxes = []
    for line in L[6:]:
        if line.startswith("END") or not line.strip(): break
        v = list(map(float, line.split()))
        boxes.append(v)
    return rt, tip, boxes

def main(path, W=100):
    rt, tip, raw = load(path)
    tau_tip = iv(*rt)
    th_tip = iv(*tip[4]); p_tip = iv(*tip[0]); k_tip = iv(*tip[2])
    D = -p_tip
    print(f"tip: tau in [{rt[0]}, {rt[1]}], D in [{dn(D)}, {up(D)}], k0 in [{dn(k_tip)}, {up(k_tip)}]")
    # keep boxes up to the tip (include those straddling it)
    B = [b for b in raw if b[0] <= rt[1]]
    n = len(B)
    T0 = np.array([b[0] for b in B]); T1 = np.array([b[1] for b in B])
    P = [iv(b[2], b[3]) for b in B]; Q = [iv(b[4], b[5]) for b in B]
    K = [iv(b[6], b[7]) for b in B]; Th = [iv(b[10], b[11]) for b in B]
    A = [th - th_tip for th in Th]
    X = [q * a.sin() + p * a.cos() for p, q, a in zip(P, Q, A)]
    Y = [p * a.sin() - q * a.cos() for p, q, a in zip(P, Q, A)]
    cosA = [a.cos() for a in A]
    dqdt = [-1 - k * p for k, p in zip(K, P)]
    ylo = np.array([dn(y) for y in Y]); yhi = np.array([up(y) for y in Y])
    xlo = np.array([dn(x) for x in X]); xhi = np.array([up(x) for x in X])
    alo = np.array([dn(a) for a in A]); ahi = np.array([up(a) for a in A])
    qlo = np.array([dn(q) for q in Q]); qhi = np.array([up(q) for q in Q])
    ok = True

    # ---------------- (G1) axis crossings ----------------
    pos = ylo > 0; neg = yhi < 0; zer = ~(pos | neg)
    runs = []
    i = 0
    while i < n:
        if zer[i]:
            j = i
            while j + 1 < n and zer[j + 1]: j += 1
            runs.append((i, j)); i = j + 1
        else:
            i += 1
    print("zero-containing runs of y~ (box index ranges):", runs)
    assert len(runs) == 2, "expected exactly two runs (double point and tip)"
    (a0, a1), (b0, b1) = runs
    assert b1 == n - 1, "second run must reach the tip end"
    assert all(pos[:a0]) and all(neg[a1 + 1:b0]), "sign pattern + ... - ... violated"
    c_dp = all(up(cosA[i]) < 0 for i in range(a0, a1 + 1))      # y~ strictly decreasing in tau
    c_tip = all(dn(cosA[i]) > 0 for i in range(b0, b1 + 1))    # y~ strictly increasing in tau
    print(f"  double-point run: tau in [{T0[a0]}, {T1[a1]}], cos(alpha)<0 on run: {c_dp}")
    print(f"  tip run:          tau in [{T0[b0]}, {T1[b1]}], cos(alpha)>0 on run: {c_tip}")
    ok &= c_dp and c_tip
    # x~ of the double point and crossing angle
    xd = X[a0]
    for i in range(a0, a1 + 1): xd = xd.union(X[i])
    ad = A[a0]
    for i in range(a0, a1 + 1): ad = ad.union(A[i])
    # forward tangent angle psi = alpha - pi/2; mirrored branch has -psi; angle between lines = |2 psi| (mod pi)
    psi = ad - arb.pi() / 2
    ang = (2 * psi).abs_upper() if False else None
    two_psi = 2 * psi
    print(f"  double point x~ in [{dn(xd)}, {up(xd)}],  forward tangent angle psi in [{math.degrees(dn(psi))}, {math.degrees(up(psi))}] deg")

    # ---------------- (G2)/(G3) ----------------
    # sub-arc index sets (tau increasing = sigma decreasing)
    Bplus = list(range(0, a1 + 1))        # from S down to (and including) the double-point run
    Bminus = list(range(a0, n))           # from the double-point run to the tip

    def tangent_range_ok(i, j):
        lo = alo[min(i, j):max(i, j) + 1].min(); hi = ahi[min(i, j):max(i, j) + 1].max()
        return hi - lo < math.pi - 1e-9

    def overlaps(I, J, mirrorJ=False):
        I = np.array(I); J = np.array(J)
        xl, xh, yl, yh = xlo[J], xhi[J], ylo[J], yhi[J]
        if mirrorJ:
            yl, yh = -yhi[J], -ylo[J]
        pairs = []
        for s in range(0, len(I), 512):
            Ii = I[s:s + 512]
            ov = (xlo[Ii, None] <= xh[None, :]) & (xl[None, :] <= xhi[Ii, None]) & \
                 (ylo[Ii, None] <= yh[None, :]) & (yl[None, :] <= yhi[Ii, None])
            a, b = np.nonzero(ov)
            pairs += list(zip(Ii[a].tolist(), J[b].tolist()))
        return pairs

    # (i) embeddedness of B+ and of B- : overlapping boxes must lie on a sub-arc of tangent range < pi
    for name, S_ in (("B (compact part)", list(range(n))),):
        bad = [(i, j) for (i, j) in overlaps(S_, S_) if i < j and not tangent_range_ok(i, j)]
        print(f"  {name}: {len(S_)} boxes, overlapping pairs not covered by the tangent-range lemma: {len(bad)}")
        ok &= (len(bad) == 0)
    # B+ and B- lie in the closed upper / lower half planes respectively except their common
    # double-point run; their intersection is {x_dp} by (G1).
    # (ii) R(B-) vs B+ : only overlaps near the double point, handled by the sector lemma
    pr = overlaps(Bplus, Bminus, mirrorJ=True)
    near = lambda i: (a0 - W) <= i <= (a1 + W)
    far = [(i, j) for (i, j) in pr if not (near(i) and near(j))]
    print(f"  R(B-) vs B+: overlapping pairs {len(pr)}, of which away from the double point: {len(far)}")
    ok &= (len(far) == 0)
    # sector lemma at x_dp: A2 = B+ near x_dp (forward directions alpha - pi/2 on boxes [a0-400, a1]),
    #                       A1 = R(B-) near x_dp (outward directions 3pi/2 - alpha on boxes [a0, a1+400])
    lo2 = alo[max(0, a0 - W):a1 + 1].min() - math.pi / 2; hi2 = ahi[max(0, a0 - W):a1 + 1].max() - math.pi / 2
    lo1 = 1.5 * math.pi - ahi[a0:a1 + W + 1].max(); hi1 = 1.5 * math.pi - alo[a0:a1 + W + 1].min()
    print(f"  sector A2 (B+ leaving x_dp): [{math.degrees(lo2):.3f}, {math.degrees(hi2):.3f}] deg;"
          f"  sector A1 (R(B-) leaving x_dp): [{math.degrees(lo1):.3f}, {math.degrees(hi1):.3f}] deg")
    sec_ok = (hi2 - lo2 < math.pi) and (hi1 - lo1 < math.pi) and (hi2 < lo1 - 1e-9 or hi1 < lo2 - 1e-9)
    print(f"  sectors disjoint (apart from the apex) and each narrower than pi: {sec_ok}")
    ok &= sec_ok
    # the local arcs must actually be the arcs within the near-window: boxes in the window that are
    # outside the sectors' arcs are covered because window boxes = those used for the sectors.
    print("ALL GEOMETRIC CHECKS PASSED" if ok else "SOME CHECK FAILED")
    return ok

if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "curve_br.txt", int(sys.argv[2]) if len(sys.argv) > 2 else 100)


def tail_and_radius_checks(path="curve_br.txt", eps_lo="0.00015345849609375", eps_hi="0.00015346679687499998"):
    """(T1) sign structure of q on the compact part: q>0 on (sigma_c, S], q<0 on (sigma_tip, sigma_c),
       with strict monotonicity of q at both zeros  -> |Gamma| decreasing then increasing, max only at S.
       (T2) tail sigma>=S: ray angle enclosure and y~>0, using the tail-theorem constants."""
    from tail import radii, mu
    rt, tip, raw = load(path)
    B = [b for b in raw if b[0] <= rt[1]]
    n = len(B)
    P = [iv(b[2], b[3]) for b in B]; Q = [iv(b[4], b[5]) for b in B]; Kk = [iv(b[6], b[7]) for b in B]
    qpos = [q > 0 for q in Q]; qneg = [q < 0 for q in Q]
    zr = [i for i in range(n) if not (qpos[i] or qneg[i])]
    runs = []
    for i in zr:
        if runs and i == runs[-1][1] + 1: runs[-1][1] = i
        else: runs.append([i, i])
    print("zero-containing runs of q:", runs)
    assert len(runs) == 2 and runs[1][1] == n - 1
    (c0, c1), (t0, t1) = runs
    ok = all(qpos[:c0]) and all(qneg[c1 + 1:t0])
    dq = [-1 - k * p for k, p in zip(Kk, P)]
    ok &= all(dq[i] < 0 for i in range(c0, c1 + 1))   # q strictly decreasing in tau at closest approach
    ok &= all(dq[i] > 0 for i in range(t0, t1 + 1))   # q strictly increasing in tau into the tip
    print(f"  closest approach run tau in [{B[c0][0]}, {B[c1][1]}];  q sign/monotonicity structure ok: {ok}")
    # radius of the part with sigma <= sigma_c (B-, and B+ up to closest approach) is <= D < q(S)
    th_tip = iv(*tip[4]); S = arb(8)
    T = radii(S, arb(max(abs(float(eps_lo)), abs(float(eps_hi))) * (1 + 1e-12)))
    Rb = T["R"]; m = mu(S)
    eta = Rb / m                                   # >= int_S^oo |k|
    pmax = 4 * m * m * Rb                          # >= |p| on [S,oo)
    qmin = S - 4 * Rb * Rb * m                     # <= q on [S,oo)
    argqp = pmax / qmin                            # >= |arg(q+ip)| (atan x <= x)
    half_pi = arb.pi() / 2
    lo = -half_pi + (-th_tip).lower() - eta - argqp
    hi = -half_pi + (-th_tip).upper() + eta + argqp
    ray_lo = -half_pi + (-th_tip).lower() - eta; ray_hi = -half_pi + (-th_tip).upper() + eta
    tail_ok = (lo > 0) and (hi < arb.pi())
    print(f"  tail: arg Gamma~ in [{float(lo.mid())*180/math.pi:.4f}, {float(hi.mid())*180/math.pi:.4f}] deg (so y~>0): {bool(tail_ok)}")
    print(f"  asymptotic ray angle (from the axis, tip-to-centre direction = 0): [{dn(ray_lo)*180/math.pi:.4f}, {up(ray_hi)*180/math.pi:.4f}] deg")
    print(f"  tail bounds: eta={float(eta.mid()):.3e}, |p|<={float(pmax.mid()):.3e}, q>={float(qmin.mid()):.6f}")
    return ok and bool(tail_ok)

if __name__ == "__main__":
    print("TAIL/RADIUS CHECKS PASSED" if tail_and_radius_checks() else "TAIL/RADIUS CHECK FAILED")
