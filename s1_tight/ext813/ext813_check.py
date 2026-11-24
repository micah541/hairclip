"""Checks of the outward S1 enclosure ext813_prof.txt (sigma in [8, 13]).  Ball arithmetic (python-flint / arb)
throughout; file values are doubles printed with 17 significant digits, read back exactly.

Reader API (for consumers):
  rows(path)        -> list of (ta, tb, tc, centre, hull); ta < tb <= 0 are tau = 8 - sigma (floats, exact),
                       centre/hull = [p, q, k, m, th] as arb balls (centre at tau = tc, hull over [ta, tb]).
  taylor(row, h)    -> (p, q, k, m, th) at sigma = sigma_c + h, h an arb ball with |h| <= half the row length:
                       f(c) + f'(c) h + f''(hull) h^2 / 2,  derivatives from the ODE (' = d/dsigma):
                       k' = m, m' = p/4, p' = -k q, q' = 1 + k p, th' = k.
Checks (written to ext813_check.json):
  0. provenance of the seed: ../in_br2.txt is rebuilt byte-for-byte by the certified generator gen_input2.build,
     with an int K1^2 ball that contains the jhat_rig enclosure of jhat8.json;
  1. tiling of tau in [-Send, 0], tc in [ta, tb], centre inside hull;
  2. widths: centre, pointwise Taylor form, hull; growth of the centre widths from sigma = 8 to 8 + Send;
  3. sigma = 8: the in_br2 parallelepiped, the ext813 point enclosure at s = 0, prof_tight row 0 and ext813 row 0
     (both transported to sigma = 8 by the Taylor form) have a common point in every component;
  4. tail envelopes: hstail.profile_bounds(sigma_a) contains every row's range of k, m, p, q - sigma;
  5. independent of CAPD: (k, k', k'') = eps*Yhat(sigma) + P(sigma) Z_nl(sigma) with Yhat = cn (K0, K1, K2)(sigma)
     exact (laplace.py), eps in the tight bracket, |Z_nl,0| <= rn0/w, |Z_nl,1|, |Z_nl,2| <= rn1/(w g) (tail2.py,
     step-1 weighted norms) -> k, m, p at the points sigma = 8 + s_j must meet the CAPD point enclosures.
Usage: python3 ext813_check.py [prof_file] [points_file]
"""
import sys, os, json, time, math
HERE = os.path.dirname(os.path.abspath(__file__))
TIGHT = os.path.dirname(HERE)
from flint import arb, acb, ctx

NAMES = ("p", "q", "k", "m", "th")


def I(lo, hi):
    """ball containing [lo, hi] (lo, hi exact doubles)"""
    return arb(lo).union(arb(hi))


def rows(path):
    out = []
    for line in open(path):
        if line.startswith("END"):
            break
        v = [float(t) for t in line.split()]
        assert len(v) == 23, "bad row"
        c = [I(v[3 + 2 * j], v[4 + 2 * j]) for j in range(5)]
        h = [I(v[13 + 2 * j], v[14 + 2 * j]) for j in range(5)]
        out.append((v[0], v[1], v[2], c, h, v))
    return out


def d1(x):
    p, q, k, m, th = x
    return [-k * q, 1 + k * p, m, p / 4, k]


def d2(x):
    p, q, k, m, th = x
    return [-(m * q + k * (1 + k * p)), m * p - k * k * q, p / 4, -k * q / 4, m]


def sq_range(h):
    """ball containing {x^2 : x in h} with nonnegative lower end"""
    a = abs(h).upper()
    return arb(0).union(arb(a) * arb(a))


def taylor(row, h):
    """enclosure of (p,q,k,m,th) at sigma = sigma_c + h (h: arb ball inside the row)"""
    c, hl = row[3], row[4]
    f1, f2 = d1(c), d2(hl)
    h2 = sq_range(h)
    return [c[i] + f1[i] * h + f2[i] * h2 / 2 for i in range(5)]


def wid(x):
    """upper bound (double) for the width of the ball x"""
    w = (x.upper() - x.lower()).upper(); f = float(w)
    while arb(f) < w:
        f = math.nextafter(f, math.inf)
    return f


def main():
    ctx.prec = 256
    prof = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "ext813_prof.txt")
    pts = sys.argv[2] if len(sys.argv) > 2 else os.path.join(HERE, "ext813_points.txt")
    res = {}
    t0 = time.time()
    # ---------------------------------------------------------------- 0. provenance of the sigma = 8 seed in_br2.txt
    # The stored in_br2.txt (seed of prof_tight.txt and of this extension) predates the final jhat8.json.  It is
    # rebuilt here byte-for-byte by the certified generator gen_input2.build (bracket2.json, current code) with
    # int K1^2 in I1old = double(1.4404409694e-9) +- 1e-18, and I1old is checked to contain the jhat_rig result.
    sys.path.insert(0, "/home/claude/s1_cap"); sys.path.insert(0, TIGHT)
    import tempfile
    import gen_input2 as G
    ctx.prec = 256
    jreal = json.load(open(os.path.join(TIGHT, "jhat8.json")))
    br = json.load(open(os.path.join(TIGHT, "bracket2.json")))
    stored = open(os.path.join(TIGHT, "in_br2.txt")).read()

    class _J:
        def __init__(self, d): self.d = d
        def load(self, f): return self.d
    def regen(i1):
        jd = dict(jreal); jd["int_K1sq"] = list(i1)
        G._cache.clear(); old = G.json; G.json = _J(jd)
        try:
            with tempfile.TemporaryDirectory() as td:
                fn = os.path.join(td, "in.txt")
                G.write(fn, *G.build("8", repr(br["eps_a"]), repr(br["eps_b"]), 8.5, emax=1.62e-4 * (1 + 1e-12))[:6])
                return open(fn).read()
        finally:
            G.json = old; G._cache.clear()
    I1OLD = (1.4404409694e-9, 1e-18)
    seed_same = regen(I1OLD) == stored
    cur = regen(jreal["int_K1sq"])
    ndiff = sum(a != b for a, b in zip(cur.split(), stored.split()))
    ball_old = arb(I1OLD[0], I1OLD[1])
    cert_str = arb(jreal["int_K1sq_str"]); cert_pair = arb(jreal["int_K1sq"][0], jreal["int_K1sq"][1])
    seed_ok = bool(seed_same and ball_old.contains(cert_str) and ball_old.contains(cert_pair))
    res["seed"] = dict(stored_in_br2_rebuilt_with_I1old=bool(seed_same), I1old=list(I1OLD),
                       I1old_contains_jhat_rig=bool(ball_old.contains(cert_str) and ball_old.contains(cert_pair)),
                       current_code_seed_differs_in_numbers=ndiff, valid=seed_ok)
    print(f"0. seed in_br2.txt rebuilt byte-for-byte by gen_input2 with I1 = double(1.4404409694e-9) +- 1e-18: {seed_same};"
          f" this ball contains the jhat_rig enclosure {jreal['int_K1sq_str']}: {seed_ok}"
          f" (current jhat8.json midpoint gives a seed differing in {ndiff} numbers, q-direction only)", flush=True)
    R = rows(prof)
    N = len(R)
    # ---------------------------------------------------------------- 1. tiling / sanity
    ok_tile = R[0][1] == 0.0 and all(R[i][0] == R[i + 1][1] for i in range(N - 1))
    Send = -R[-1][0]
    ok_tc = all(r[0] < r[1] and r[0] <= r[2] <= r[1] for r in R)
    ok_in = all(r[4][j].lower() <= r[3][j].lower() and r[3][j].upper() <= r[4][j].upper() for r in R for j in range(5))
    res["rows"] = N; res["sigma_range"] = [8.0, 8.0 + Send]
    res["tiling_ok"] = bool(ok_tile); res["tc_in_row_ok"] = bool(ok_tc); res["centre_in_hull_ok"] = bool(ok_in)
    print(f"1. {N} rows, tau in [{-Send}, 0] (sigma in [8, {8 + Send}]): tiling {ok_tile}, tc in row {ok_tc}, "
          f"centre in hull {ok_in}", flush=True)
    # ---------------------------------------------------------------- 2. widths
    wc = [0.0] * 5; wp = [0.0] * 5; wh = [0.0] * 5; wq_row = 0.0
    for r in R:
        ta, tb, tc, c, hl, _ = r
        Hb = arb(max(abs(arb(tc) - arb(ta)).upper(), abs(arb(tb) - arb(tc)).upper()))
        f1, f2 = d1(c), d2(hl)
        for j in range(5):
            wc[j] = max(wc[j], wid(c[j]))
            pw = arb(wid(c[j])) + Hb * arb(wid(f1[j])) + Hb * Hb / 2 * arb(wid(f2[j]))
            wp[j] = max(wp[j], float(pw.upper()))
            wh[j] = max(wh[j], wid(hl[j]))
    P = [l.split() for l in open(pts)]
    Pt = [(float(v[0]), [I(float(v[1 + 2 * j]), float(v[2 + 2 * j])) for j in range(5)]) for v in P]
    w0 = [wid(x) for x in Pt[0][1]]; w1 = [wid(x) for x in Pt[-1][1]]
    res["max_centre_width"] = dict(zip(NAMES, wc))
    res["max_pointwise_width"] = dict(zip(NAMES, wp))
    res["max_hull_width"] = dict(zip(NAMES, wh))
    res["point_widths"] = {f"{8 + s:g}": dict(zip(NAMES, [wid(x) for x in v])) for s, v in Pt}
    res["growth_8_to_end"] = {n: (w1[j] / w0[j] if w0[j] > 0 else None) for j, n in enumerate(NAMES)}
    A_E = arb(3) / 4 * arb(4) ** (-arb(1) / 3)
    dE = A_E * (arb(8 + Send) ** (arb(4) / 3) - arb(8) ** (arb(4) / 3))
    mu_ratio = ((arb(8 + Send) / 8) ** (arb(1) / 3))
    res["wkb_growth"] = dict(eE2=float((dE / 2).exp().mid()), k=float(((dE / 2).exp() / mu_ratio).mid()),
                             m=float((dE / 2).exp().mid()), p=float(((dE / 2).exp() * mu_ratio).mid()))
    print("2. max centre width   ", {n: f"{w:.2e}" for n, w in zip(NAMES, wc)})
    print("   max pointwise width", {n: f"{w:.2e}" for n, w in zip(NAMES, wp)}, "(2nd-order Taylor form, any sigma)")
    print("   max hull width     ", {n: f"{w:.2e}" for n, w in zip(NAMES, wh)})
    print("   growth of centre width sigma=8 ->", 8 + Send, {n: (f"{g:.1f}" if g else "-") for n, g in res["growth_8_to_end"].items()},
          " WKB e^{dE/2} (x mu^-1, 1, mu for k, m, p):", {k: f"{v:.1f}" for k, v in res["wkb_growth"].items()}, flush=True)
    # ---------------------------------------------------------------- 3. sigma = 8 against in_br2 and prof_tight
    L = open(os.path.join(TIGHT, "in_br2.txt")).read().split()
    vals = [float(t) for t in L[1:]]
    xs = [I(vals[2 * i], vals[2 * i + 1]) for i in range(5)]
    Cm = [[arb(vals[10 + 5 * i + j]) for j in range(5)] for i in range(5)]
    r0 = [I(vals[35 + 2 * j], vals[36 + 2 * j]) for j in range(5)]
    par = [xs[i] + sum((Cm[i][j] * r0[j] for j in range(5)), arb(0)) for i in range(5)]
    pt8 = Pt[0][1]; assert Pt[0][0] == 0.0
    first = open(os.path.join(TIGHT, "prof", "prof_tight.txt")).readline().split()
    vP = [float(t) for t in first]
    rowP = (vP[0], vP[1], vP[2], [I(vP[3 + 2 * j], vP[4 + 2 * j]) for j in range(5)],
            [I(vP[13 + 2 * j], vP[14 + 2 * j]) for j in range(5)], vP)
    # prof_tight row 0: centre at sigma_c = 8 - tc;  sigma = 8  <->  h = 8 - sigma_c = tc
    trP = taylor(rowP, arb(vP[2]))
    trE = taylor(R[0], arb(R[0][2]))           # ext813 row 0: h = 8 - (8 - tc) = tc (negative)
    hullmeet = all(max(R[0][4][j].lower(), rowP[4][j].lower()) <= min(R[0][4][j].upper(), rowP[4][j].upper())
                   for j in range(5))
    encl = {"in_br2": par, "ext_s0": pt8, "prof_row0_to_8": trP, "ext_row0_to_8": trE}
    common = {}
    ok8 = True
    for j, n in enumerate(NAMES):
        lo = max(e[j].lower() for e in encl.values()); hi = min(e[j].upper() for e in encl.values())
        meet = bool(lo <= hi); ok8 = ok8 and meet
        common[n] = dict(meet=meet, common_width=float((hi - lo).upper()) if meet else None,
                         widths={k: wid(e[j]) for k, e in encl.items()})
    par_in_pt = all(pt8[j].lower() <= par[j].lower() and par[j].upper() <= pt8[j].upper() for j in range(5))
    res["sigma8"] = dict(all_four_meet=ok8, hull_rows_meet=bool(hullmeet), in_br2_inside_ext_s0=bool(par_in_pt), detail=common)
    print(f"3. sigma=8: in_br2 / ext(s=0) / prof_tight row0 -> 8 / ext row0 -> 8 have a common point in every component: {ok8};"
          f" hull rows meet: {hullmeet}; in_br2 hull inside ext(s=0): {par_in_pt}")
    for n in NAMES:
        d = common[n]
        print(f"   {n:2s} common width {d['common_width']:.2e}  widths " +
              " ".join(f"{k}={w:.2e}" for k, w in d["widths"].items()))
    sys.stdout.flush()
    # ---------------------------------------------------------------- 4. tail envelopes
    sys.path.insert(0, "/home/claude/s1spec/step3")
    import hstail as Hs
    ctx.prec = 256
    worst = {n: (0.0, None) for n in ("k", "m", "p", "q-sigma")}
    ok_env = True
    for r in R:
        ta, tb, tc = r[0], r[1], r[2]
        sa = arb(8) - arb(tb); sc = arb(8) - arb(tc)
        hb = (arb(8) - arb(ta) - sc).union(arb(8) - arb(tb) - sc)   # h = sigma - sigma_c over the row (sigma = 8 - tau)
        kb, mb, pb, qs = Hs.profile_bounds(sa)
        T = taylor(r, hb)
        c, hl = r[3], r[4]
        qms = (c[1] - sc) + (c[2] * c[0]) * hb + (hl[3] * hl[0] - hl[2] * hl[2] * hl[1]) * sq_range(hb) / 2
        for n, val, B in (("k", T[2], kb), ("m", T[3], mb), ("p", T[0], pb), ("q-sigma", qms, qs)):
            ratio = float((abs(val).upper() / B.lower()).upper())
            if ratio > worst[n][0]:
                worst[n] = (ratio, float(sa.mid()))
            if not (abs(val).upper() < B.lower()):
                ok_env = False
        for n, j, B in (("k", 2, kb), ("m", 3, mb), ("p", 0, pb)):
            if not abs(hl[j]).upper() < B.lower():
                ok_env = False
    res["envelopes"] = dict(inside=ok_env, worst_ratio={n: dict(ratio=v[0], sigma_a=v[1]) for n, v in worst.items()})
    print(f"4. tail envelopes hstail.profile_bounds(sigma_a): every row inside: {ok_env}; worst |value|/bound: " +
          ", ".join(f"{n} {v[0]:.4f} (sigma {v[1]:.3f})" for n, v in worst.items()), flush=True)
    # ---------------------------------------------------------------- 5. Laplace / Lyapunov-Perron representation
    sys.path.insert(0, "/home/claude/s1_cap"); sys.path.insert(0, TIGHT)
    import tail as T1
    import tail2
    import laplace as LA
    ctx.prec = 256
    br = json.load(open(os.path.join(TIGHT, "bracket2.json")))
    EPS = arb(repr(br["eps_a"])).union(arb(repr(br["eps_b"])))
    NL = tail2.nl_radii(arb(8), arb(1.62e-4 * (1 + 1e-12)))
    assert NL["verified"]
    rn0, rn1 = NL["rn0"], NL["rn1"]
    Ks8 = [acb(LA.K(8, j)) for j in range(3)]
    Z8 = tail2.solve3(tail2.Pmat(arb(8)), Ks8)
    cnc = 1 / Z8[0]
    cn = cnc.real + arb(0, abs(cnc.imag).upper())
    lap = {}
    ok_lap = True
    for s, v in Pt:
        sig = arb(8) + arb(s)
        Kv = [LA.K(sig, j) for j in range(3)]
        mu_ = T1.mu(sig); g = (sig / 8) ** (arb(4) / 3)
        w = (T1.A_E * (sig ** (arb(4) / 3) - arb(8) ** (arb(4) / 3))).exp() * (sig / 8) ** (arb(1) / 3)
        rem = (rn0 + 2 * rn1 / g) / w
        kL = EPS * cn * Kv[0] + arb(0, rem.upper())
        mL = EPS * cn * Kv[1] + arb(0, (mu_ * rem).upper())
        pL = 4 * (EPS * cn * Kv[2] + arb(0, (mu_ * mu_ * rem).upper()))
        ent = {}
        for n, xl, xc in (("k", kL, v[2]), ("m", mL, v[3]), ("p", pL, v[0])):
            meet = bool(max(xl.lower(), xc.lower()) <= min(xl.upper(), xc.upper()))
            inside = bool(xc.lower() <= xl.lower() and xl.upper() <= xc.upper())
            ok_lap = ok_lap and meet
            ent[n] = dict(meet=meet, laplace_inside_capd=inside, w_laplace=wid(xl), w_capd=wid(xc),
                          offset_mid=float(abs(xl.mid() - xc.mid())))
        lap[f"{float(sig.mid()):g}"] = ent
    keys = list(lap); inside_from = None
    for i in range(len(keys)):
        if all(lap[kk][n]["laplace_inside_capd"] for kk in keys[i:] for n in ("k", "m", "p")):
            inside_from = float(keys[i]); break
    res["laplace"] = dict(all_meet=ok_lap, laplace_inside_capd_for_sigma_ge=inside_from, rn0=float(rn0.upper()),
                          rn1=float(rn1.upper()), points=lap)
    print(f"5. Laplace+LP (k, m, p) at {len(Pt)} points sigma in [8, {8 + Send}]: all meet the CAPD enclosure: {ok_lap}; "
          f"Laplace ball inside the CAPD ball for sigma >= {inside_from} (below, the LP bound rn0 + 2 rn1 is the wider one)")
    for key in list(lap)[::4] + ([list(lap)[-1]] if (len(lap) - 1) % 4 else []):
        e = lap[key]
        print(f"   sigma={key:>5}: " + "  ".join(f"{n}: w_lap {e[n]['w_laplace']:.1e} w_capd {e[n]['w_capd']:.1e} |dmid| {e[n]['offset_mid']:.1e}"
                                          for n in ("k", "m", "p")))
    res["all_checks_pass"] = bool(seed_ok and ok_tile and ok_tc and ok_in and ok8 and hullmeet and ok_env and ok_lap)
    res["seconds"] = time.time() - t0
    json.dump(res, open(os.path.join(HERE, "ext813_check.json"), "w"), indent=1)
    print(f"ALL CHECKS PASS: {res['all_checks_pass']}   [{res['seconds']:.0f}s]  -> ext813_check.json")
    return 0 if res["all_checks_pass"] else 1


if __name__ == "__main__":
    sys.exit(main())
