"""
Certified winding numbers from the segment enclosures in contour_result.json.

Each segment i carries, for each class, an affine enclosure  c_i D(lam) in X + V s + E  (s in [-1/2,1/2], E a box),
where c_i > 0 is an (unknown) positive scale used inside the CAPD run (powers of two).  Let S_i = X + V[-1/2,1/2] + E
(convex, 0 notin S_i, checked rigorously in contour.py) and K_i = {r z : r > 0, z in S_i}, a convex cone with 0 notin K_i
containing D(segment i).  We pick unit vectors alpha_i in K_{i-1} cap K_i (both cones contain the ray of D(lam_i)):
the argument interval of the end box of S_{i-1} is intersected with that of the start box of S_i.  Then
      sum_i pv arg(alpha_{i+1}/alpha_i) = change of arg D along the path,
because consecutive alphas lie in a common convex cone avoiding 0.  At the real endpoints alpha = +-1 (D real there).
All arithmetic on the box corners is done in Arb; the representatives alpha_i are exact angles chosen inside the
certified intersections.
"""
import json, math, sys
from flint import arb, acb


def box_of(S, s):
    X = complex(*S["X"]); V = complex(*S["V"]); Er, Ei = S["E"]
    c = X + V * s
    return (c.real - Er, c.real + Er, c.imag - Ei, c.imag + Ei)


def contains_zero(b):
    return b[0] <= 0 <= b[1] and b[2] <= 0 <= b[3]


def arg_interval(b, ref):
    """rigorous enclosure [lo, hi] of arg over the box b (0 notin b), on the branch within pi of ref (float)"""
    assert not contains_zero(b)
    args = []
    for x in (b[0], b[1]):
        for y in (b[2], b[3]):
            a = acb(arb(x), arb(y)).arg()     # principal arg in (-pi, pi], as a ball
            # move to the branch near ref
            k = round((ref - float(a.mid())) / (2 * math.pi))
            args.append(a + 2 * k * arb.pi())
    lo = min(float(a.lower()) for a in args); hi = max(float(a.upper()) for a in args)
    assert hi - lo < math.pi, "box subtends >= pi"
    return lo, hi


def winding_from(done, cls):
    # start: real endpoint
    b0 = box_of(done[0][cls], -0.5)
    assert b0[2] <= 0 <= b0[3], "start box does not meet the real axis"
    assert b0[0] > 0 or b0[1] < 0, "sign at the start not determined"
    theta = [0.0 if b0[0] > 0 else math.pi]
    for r1, r2 in zip(done[:-1], done[1:]):
        ref = theta[-1]
        b1 = box_of(r1[cls], 0.5); b2 = box_of(r2[cls], -0.5)
        I1 = arg_interval(b1, ref); I2 = arg_interval(b2, 0.5 * (I1[0] + I1[1]))
        lo, hi = max(I1[0], I2[0]), min(I1[1], I2[1])
        assert lo <= hi, f"arg intervals do not meet at lam = {r1['lb']}"
        theta.append(0.5 * (lo + hi))
    bN = box_of(done[-1][cls], 0.5)
    assert bN[2] <= 0 <= bN[3], "end box does not meet the real axis"
    assert bN[0] > 0 or bN[1] < 0, "sign at the end not determined"
    thN = 0.0 if bN[0] > 0 else math.pi
    # end representative: the real ray, on the branch of the last theta
    k = round((theta[-1] - thN) / (2 * math.pi))
    thN += 2 * math.pi * k
    I = arg_interval(bN, theta[-1])
    assert I[0] <= thN <= I[1], "end ray not in the last box's argument interval"
    theta.append(thN)
    total = 0.0; maxstep = 0.0
    for a, b in zip(theta[:-1], theta[1:]):
        d = b - a
        d = (d + math.pi) % (2 * math.pi) - math.pi      # principal value
        maxstep = max(maxstep, abs(d))
        total += d
    return dict(arg_change_over_pi=total / math.pi, winding=round(total / math.pi), max_step=maxstep,
                segments=len(done), sign_start=1 if b0[0] > 0 else -1, sign_end=1 if bN[0] > 0 else -1)


if __name__ == "__main__":
    R = json.load(open(sys.argv[1] if len(sys.argv) > 1 else "contour_result.json"))
    done = R["segments"]
    for r in done:
        r["la"] = complex(*r["la"]); r["lb"] = complex(*r["lb"])
        for cls in ("even", "odd"):
            assert r[cls]["ok"], "a segment without a certified enclosure"
    for r1, r2 in zip(done[:-1], done[1:]):
        assert abs(r1["lb"] - r2["la"]) < 1e-12
    out = {cls: winding_from(done, cls) for cls in ("even", "odd")}
    print(json.dumps(out, indent=1))
    json.dump(out, open("winding_report.json", "w"), indent=1)
