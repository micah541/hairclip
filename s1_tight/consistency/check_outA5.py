"""Do the step-3 hulls of the wide CAPD run (s1spec/step3/outA5.txt, as read by step3nl/agent_mod/ibp/ibp_boxes.py)
contain the tight S1 enclosures (ext_tip/prof_tip.txt, prof/prof_tight.txt, ext813/ext813_prof.txt) wherever they overlap?
Ball arithmetic (python-flint / arb) throughout; file numbers are doubles read back exactly.

outA5 rows, in ibp_boxes' conventions (gen_hs3.build_A2):
  MC (1025 rows, stage C, sigma in [sigma_tip, 8]): rescaled time t in [ta, tb], sigma = 8 - c t, c = (8 - st)/T0,
     T0 = 10.249073;  h = (p, q, k, m, st, C, S, Ct2, St2, ws0..ws3, integrals),  (Ct2, St2) = (cos th, -sin th),
     th as in prof_tight (th(8) = 0, dth/dtau = -k).  For S1 itself st = sigma_tip, so the row covers tau = 8 - sigma
     in (tau_tip/T0) [ta, tb].
  MB (2400 rows, stage B, sigma in [8, 32]): t in [ta, tb], sigma = 32 - t, tau = t - 24;  h = (st, C, S, ...),
     (C, S) = (cos K, sin K), K = int_sigma^32 k;  profile = chunk envelope gen_hs2.chunks()[int(ta)] = (k, m, p, q - sigma),
     and q in sigma + dq, exactly as ibp_boxes builds kv, mv, pv_, qv.
Test: on each overlap, the tight enclosure restricted to the overlap (2nd-order Taylor form f(c) + f'(c) h + f''(hull) h^2/2
on every tight row meeting it, then the union) lies inside the outA5 / ibp box, component by component:
  prof_tip, prof_tight vs MC: p, q, k, m;  cos th in Ct2, -sin th in St2;  sigma_tip in st.
  ext813 vs MB:  k, m, p, q - sigma in the chunk envelope, q in sigma + dq;  cos K in C, sin K in S with
     K = th(13) - th(sigma) + int_13^32 k,  |int_13^32 k| <= sum kb(s_i) ds  (kb = hstail.profile_bounds, decreasing);
     sigma_tip in st.
Usage: python3 check_outA5.py      -> check_outA5.json"""
import sys, os, json, time, bisect, math
from flint import arb, ctx
HERE = os.path.dirname(os.path.abspath(__file__)); TIGHT = os.path.dirname(HERE)
S1SPEC = "/home/claude/s1spec"
OUTA5 = os.path.join(S1SPEC, "step3", "outA5.txt")
ctx.prec = 128
NAMES = ("p", "q", "k", "m", "th")


def I(lo, hi): return arb(lo).union(arb(hi))


def f1(x):          # d/dtau, field of s1prof
    p, q, k, m, th = x
    return [k * q, -1 - k * p, -m, -p / 4, -k]


def f2(x):          # d^2/dtau^2
    p, q, k, m, th = x
    return [-m * q - k - k * k * p, m * p - k * k * q, p / 4, -k * q / 4, m]


def load(path):
    """tight rows sorted by ta: (ta, tb, tc, centre, f'(centre), f''(hull), (q-sigma)' (centre), (q-sigma)'' (hull))"""
    out = []
    for line in open(path):
        if line.startswith("END"): break
        v = [float(t) for t in line.split()]
        c = [I(v[3 + 2 * j], v[4 + 2 * j]) for j in range(5)]
        hl = [I(v[13 + 2 * j], v[14 + 2 * j]) for j in range(5)]
        p, q, k, m, th = c; P, Q, K, Mm, T = hl
        out.append((v[0], v[1], v[2], c, f1(c), f2(hl), -k * p, Mm * P - K * K * Q))
    out.sort(key=lambda r: r[0])
    return out


def quad_range(c, a, b, s0, s1):
    """[lo, hi] (exact arb) enclosing {c' + a' h + b' h^2/2 : c' in c, a' in a, b' in b, h in [s0, s1]}: the Lagrange form
    f(tc + h) = f(tc) + f'(tc) h + f''(xi) h^2/2 with f(tc) in c, f'(tc) in a, f''(xi) in b, evaluated without decoupling h
    and h^2.  For each h the value lies in [L(h), U(h)]; L, U are quadratics, minimised / maximised exactly on [s0, s1]
    (split at h = 0, where the sign of a h changes; vertex included only when it can lie in the piece)."""
    pieces = [(s0, arb(0)), (arb(0), s1)] if (s0 < 0 and s1 > 0) else [(s0, s1)]
    cl, cu, bl, bu = c.lower(), c.upper(), b.lower(), b.upper()
    lo = hi = None
    for u0, u1 in pieces:
        if u0 >= 0: aL, aU = a.lower(), a.upper()
        else: aL, aU = a.upper(), a.lower()
        L = lambda x: cl + aL * x + bl * x * x / 2
        U = lambda x: cu + aU * x + bu * x * x / 2
        cl_ = [L(u0).lower(), L(u1).lower()]; cu_ = [U(u0).upper(), U(u1).upper()]
        if bl > 0:
            hs = -aL / bl
            if not (hs.upper() < u0 or hs.lower() > u1): cl_.append((cl - aL * aL / (2 * bl)).lower())
        if bu < 0:
            hs = -aU / bu
            if not (hs.upper() < u0 or hs.lower() > u1): cu_.append((cu - aU * aU / (2 * bu)).upper())
        for v in cl_:
            if lo is None or v < lo: lo = v
        for v in cu_:
            if hi is None or v > hi: hi = v
    return [lo, hi]


def enclose(R, tas, Tlo, Thi):
    """enclosures over tau in [Tlo, Thi] (exact arb): [p, q, k, m, th] and q - sigma (sigma = 8 - tau)"""
    i = max(bisect.bisect_right(tas, float(Tlo)) - 2, 0)
    acc = None
    while i < len(R) and arb(R[i][0]) <= Thi:
        ta, tb, tc, c, d1, d2, g1, g2 = R[i]; i += 1
        ja = Tlo if Tlo > arb(ta) else arb(ta)
        jb = Thi if Thi < arb(tb) else arb(tb)
        if ja > jb: continue
        h0, h1 = ja - arb(tc), jb - arb(tc)
        vals = [quad_range(c[j], d1[j], d2[j], h0, h1) for j in range(5)]
        vals.append(quad_range(c[1] - (8 - arb(tc)), g1, g2, h0, h1))      # q - sigma
        if acc is None: acc = vals
        else:
            for a, v in zip(acc, vals):
                if v[0] < a[0]: a[0] = v[0]
                if v[1] > a[1]: a[1] = v[1]
    if acc is None: raise RuntimeError("no tight row covers [%r, %r]" % (float(Tlo), float(Thi)))
    return acc


def inside(lohi, box):
    return bool(box.lower() <= lohi[0] and lohi[1] <= box.upper())


def ball(lohi): return lohi[0].union(lohi[1])


def trig_range(lohi, fn):
    """tight enclosure of {cos t} (fn = 'cos') or {sin t} (fn = 'sin') for t in [a, b]: endpoint values at full precision
    plus the critical points k pi (cos) or pi/2 + k pi (sin) that may lie in [a, b].  (arb's cos/sin of a ball with poor
    relative accuracy is evaluated at low precision, ~2^-30 radius, which is too coarse for the outA5 boxes.)"""
    a, b = lohi
    f = (lambda x: x.cos()) if fn == "cos" else (lambda x: x.sin())
    vals = [f(a), f(b)]
    pi = arb.pi(); off = arb(0) if fn == "cos" else pi / 2
    k0 = math.floor(float(((a - off) / pi).mid())) - 2; k1 = math.ceil(float(((b - off) / pi).mid())) + 2
    for kk in range(k0, k1 + 1):
        cpt = off + kk * pi
        if not (cpt.upper() < a or cpt.lower() > b):
            vals.append(arb(1 if kk % 2 == 0 else -1))
    lo = vals[0].lower(); hi = vals[0].upper()
    for v in vals[1:]:
        if v.lower() < lo: lo = v.lower()
        if v.upper() > hi: hi = v.upper()
    return [lo, hi]


def margin(lohi, box):
    """smallest distance from the tight range to a face of the box, relative to the box width"""
    w = box.upper() - box.lower()
    if not w > 0: return 0.0
    return float((min(lohi[0] - box.lower(), box.upper() - lohi[1]) / w).mid())


def excess(lohi, box):
    """how far the tight range sticks out of the box: (absolute, relative to the box width)"""
    out = max(box.lower() - lohi[0], lohi[1] - box.upper())
    w = box.upper() - box.lower()
    return float(out.mid()), (float((out / w).mid()) if w > 0 else float("inf"))


def summ(fails):
    return {n: dict(rows=len(f), t_range=[min(z[0] for z in f), max(z[1] for z in f)],
                    max_excess_abs=max(z[2][0] for z in f), max_excess_rel=max(z[2][1] for z in f)) for n, f in fails.items()}


def parse_outA5():
    MB, MC = [], []
    for line in open(OUTA5):
        tag = line[:2]
        if tag not in ("MB", "MC"): continue
        v = [float(t) for t in line.split()[1:]]
        h = [I(v[2 + 2 * i], v[3 + 2 * i]) for i in range((len(v) - 2) // 2)]
        (MB if tag == "MB" else MC).append((v[0], v[1], h))
    return MB, MC


def check_MC(name, path, tau_tip, MC, T0):
    R = load(path); tas = [r[0] for r in R]
    cst = tau_tip / T0
    stip = 8 - tau_tip
    ok = True; nfail = 0; worst = {n: 1.0 for n in ("p", "q", "k", "m", "Ct2", "St2")}; st_ok = True; fails = {}
    for ta, tb, h in MC:
        Tlo = (cst * arb(ta)).lower(); Thi = (cst * arb(tb)).upper()
        if Thi > arb(R[-1][1]): Thi = arb(R[-1][1])           # tight rows end at Tend >= tau_tip
        e = enclose(R, tas, Tlo, Thi)
        tests = [("p", e[0], h[0]), ("q", e[1], h[1]), ("k", e[2], h[2]), ("m", e[3], h[3])]
        cth = trig_range(e[4], "cos"); sth = trig_range(e[4], "sin")
        tests += [("Ct2", cth, h[7]), ("St2", [-sth[1], -sth[0]], h[8])]
        row_ok = True
        for n, lohi, box in tests:
            if not inside(lohi, box):
                row_ok = False; fails.setdefault(n, []).append((ta, tb, excess(lohi, box)))
            else: worst[n] = min(worst[n], margin(lohi, box))
        st_ok &= bool(h[4].lower() <= stip.lower() and stip.upper() <= h[4].upper())
        if not row_ok: nfail += 1
        ok &= row_ok
    return dict(file=name, rows_compared=len(MC), all_inside=bool(ok and st_ok), rows_failing=nfail,
                sigma_tip_inside_st=bool(st_ok), min_relative_margin=worst, failures=summ(fails))


def check_MB(MB, tau_tips):
    sys.path.insert(0, os.path.join(S1SPEC, "step3"))
    import gen_hs2 as G2          # chunks() as used by ibp_boxes (mod_gram imports it from gen_hs2)
    import hstail as Hs
    ctx.prec = 128
    CH = G2.chunks()
    R = load(os.path.join(TIGHT, "ext813", "ext813_prof.txt")); tas = [r[0] for r in R]
    tlo_all, thi_all = arb(R[0][0]), arb(R[-1][1])            # tau in [-5, 0]
    # th(13): ext813 point enclosure at s = 5 (sigma = 13)
    pts = [l.split() for l in open(os.path.join(TIGHT, "ext813", "ext813_points.txt"))]
    last = [float(t) for t in pts[-1]]; assert last[0] == 5.0
    th13 = I(last[9], last[10])
    # |int_13^32 k| <= sum_i kb(s_i) * ds, s_i = 13 + i ds  (kb decreasing in sigma)
    ds = arb(1) / 20; Bt = arb(0); s = arb(13)
    while s < 32:
        Bt += Hs.profile_bounds(s)[0] * ds; s += ds
    Itail = arb(0, Bt.upper())
    ok = True; nfail = 0; n = 0; worst = {k: 1.0 for k in ("k", "m", "p", "q-sigma", "q", "C", "S")}; st_ok = True; fails = {}
    for ta, tb, h in MB:
        Tlo, Thi = arb(ta) - 24, arb(tb) - 24                    # tau = 8 - sigma = t - 24
        if Thi < tlo_all or Tlo > thi_all: continue
        Tlo = Tlo if Tlo > tlo_all else tlo_all; Thi = Thi if Thi < thi_all else thi_all
        n += 1
        e = enclose(R, tas, Tlo, Thi)
        c = min(int(ta + 1e-9), len(CH) - 1)
        kv, mv, pv, dq = [I(lo, hi) for lo, hi in CH[c][2]]
        sig = I(32 - tb, 32 - ta); qv = sig + dq                  # the ibp_boxes q box of this row
        K = th13 - ball(e[4]) + Itail
        Kr = [K.lower(), K.upper()]
        tests = [("k", e[2], kv), ("m", e[3], mv), ("p", e[0], pv), ("q-sigma", e[5], dq), ("q", e[1], qv),
                 ("C", trig_range(Kr, "cos"), h[1]), ("S", trig_range(Kr, "sin"), h[2])]
        row_ok = True
        for nm, lohi, box in tests:
            if not inside(lohi, box):
                row_ok = False; fails.setdefault(nm, []).append((ta, tb, excess(lohi, box)))
            else: worst[nm] = min(worst[nm], margin(lohi, box))
        for tt in tau_tips:
            st = 8 - tt
            st_ok &= bool(h[0].lower() <= st.lower() and st.upper() <= h[0].upper())
        if not row_ok: nfail += 1
        ok &= row_ok
    return dict(file="ext813/ext813_prof.txt", rows_compared=n, all_inside=bool(ok and st_ok), rows_failing=nfail,
                sigma_tip_inside_st=bool(st_ok), int_13_32_k_bound=float(Bt.upper()), min_relative_margin=worst,
                failures=summ(fails))


def main():
    t0 = time.time()
    MB, MC = parse_outA5()
    t0par = [l.split() for l in open(os.path.join(S1SPEC, "step3", "inA5.txt")) if l.startswith("T0 ")][0]
    assert float(t0par[1]) == float(t0par[2]) == MC[-1][1] == 10.249073 and MC[0][0] == 0.0   # stage C: t in [0, T0]
    T0 = arb(10.249073)
    tip = json.load(open(os.path.join(TIGHT, "ext_tip", "ext_tip_check.json")))["tip"]["tau_tip"]
    tt_tip = I(*tip)                                              # conditioned tau_tip (ext_tip)
    tt_prof = I(10.249077650575472, 10.249077651519443)           # box tau_tip (s1_tight, s1curve)
    res = dict(outA5=OUTA5, MC_rows=len(MC), MB_rows=len(MB))
    res["prof_tip_vs_MC"] = check_MC("ext_tip/prof_tip.txt", os.path.join(TIGHT, "ext_tip", "prof_tip.txt"), tt_tip, MC, T0)
    print("prof_tip   vs MC:", res["prof_tip_vs_MC"], flush=True)
    res["prof_tight_vs_MC"] = check_MC("prof/prof_tight.txt", os.path.join(TIGHT, "prof", "prof_tight.txt"), tt_prof, MC, T0)
    print("prof_tight vs MC:", res["prof_tight_vs_MC"], flush=True)
    res["ext813_vs_MB"] = check_MB(MB, [tt_tip, tt_prof])
    print("ext813     vs MB:", res["ext813_vs_MB"], flush=True)
    res["all_inside"] = bool(res["prof_tip_vs_MC"]["all_inside"] and res["prof_tight_vs_MC"]["all_inside"]
                             and res["ext813_vs_MB"]["all_inside"])
    res["seconds"] = time.time() - t0
    json.dump(res, open(os.path.join(HERE, "check_outA5.json"), "w"), indent=1)
    print(f"ALL TIGHT ENCLOSURES INSIDE THE outA5 / ibp_boxes HULLS: {res['all_inside']}   [{res['seconds']:.0f}s]")
    return 0 if res["all_inside"] else 1


if __name__ == "__main__":
    sys.exit(main())
