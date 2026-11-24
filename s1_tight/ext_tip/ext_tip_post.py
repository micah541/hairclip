"""Post-processing and checks for the conditioned S1 profile (ball arithmetic, python-flint/arb; outward rounding).
Input  prof_cond.txt (s1cond: conditioned centres, s1prof hulls) and tip_cond.txt.
Output prof_tip.txt: same rows, with each hull replaced by
       hull_new = (c + f'(c) h + f''(hull_old) [0, H^2]/2)  intersected with  hull_old,   h = tau - tc in [ta-tc, tb-tc],
       ' = d/dtau (field of s1prof), c = conditioned centre; valid because hull_old is valid on the row.
Checks (ext_tip_check.json):
  1. grid identical to prof_tight; centre inside prof_tight centre; hull_new inside hull_old; centre inside hull_new;
  2. widths (centre, pointwise Taylor form) overall and on the tip segment, against prof_tight;
  3. tip: the last row transported to the conditioned tau_tip meets the conditioned tip state and contains q = m = 0;
     the tip data lie inside the s1_tight intervals (tau_tip, D = -p, k(0));
  4. float sanity check (NOT rigorous): a float trajectory from a zero of M in Y lies within the enclosure.
"""
import sys, os, json, math, time
from flint import arb, ctx
HERE = os.path.dirname(os.path.abspath(__file__)); TIGHT = os.path.dirname(HERE)
ctx.prec = 128
NAMES = ("p", "q", "k", "m", "th")


def I(lo, hi): return arb(lo).union(arb(hi))


def dn(a):
    lo = a.lower(); v = float(lo)
    while arb(v) > lo: v = math.nextafter(v, -math.inf)
    return v


def up(a):
    hi = a.upper(); v = float(hi)
    while arb(v) < hi: v = math.nextafter(v, math.inf)
    return v


def f1(x):          # d/dtau
    p, q, k, m, th = x
    return [k * q, -1 - k * p, -m, -p / 4, -k]


def f2(x):          # d^2/dtau^2
    p, q, k, m, th = x
    return [-m * q - k - k * k * p, m * p - k * k * q, p / 4, -k * q / 4, m]


def wid(x): return float((x.upper() - x.lower()).upper())


def main():
    t0 = time.time()
    tip = {}
    for line in open(os.path.join(HERE, "tip_cond.txt")):
        w = line.split()
        if len(w) == 3 and w[0] != "intersections_nonempty":
            try: tip[w[0]] = (float(w[1]), float(w[2]))
            except ValueError: pass
    fc = open(os.path.join(HERE, "prof_cond.txt")); fp = open(os.path.join(TIGHT, "prof", "prof_tight.txt"))
    out = open(os.path.join(HERE, "prof_tip.txt"), "w")
    TT = I(*tip["tau_tip"])
    res = dict(rows=0)
    ok_grid = ok_in = ok_hull = ok_c_in_h = ok_boxhull = True
    wc = [0.0] * 5; wpt = [0.0] * 5; wprof = [0.0] * 5; wseg = {}; whull = [0.0] * 5; whull_old = [0.0] * 5
    SEG = [(0.0, 2.0), (2.0, 5.0), (5.0, 10.3)]          # x = tau_tip - tau ranges (tip segment first)
    for s in SEG: wseg[s] = dict(cond=[0.0] * 5, prof=[0.0] * 5)
    last = None
    for la, lb in zip(fc, fp):
        if la.startswith("END"):
            assert lb.startswith("END"); break
        a = la.split(); b = lb.split()
        ok_grid &= a[:3] == b[:3]
        ok_boxhull &= a[13:] == b[13:]          # s1cond's box-run hulls = prof_tight hulls, byte for byte
        ta, tb, tc = float(a[0]), float(a[1]), float(a[2])
        c = [I(float(a[3 + 2 * j]), float(a[4 + 2 * j])) for j in range(5)]
        ho = [I(float(a[13 + 2 * j]), float(a[14 + 2 * j])) for j in range(5)]
        cp = [I(float(b[3 + 2 * j]), float(b[4 + 2 * j])) for j in range(5)]
        h = (arb(ta) - arb(tc)).union(arb(tb) - arb(tc))
        H = arb(max(abs(arb(ta) - arb(tc)).upper(), abs(arb(tb) - arb(tc)).upper()))
        H2 = arb(0).union(H * H)
        d1, d2 = f1(c), f2(ho)
        row = [a[0], a[1], a[2]]
        for j in range(5):
            row += [a[3 + 2 * j], a[4 + 2 * j]]
            ok_in &= bool(cp[j].lower() <= c[j].lower() and c[j].upper() <= cp[j].upper())
        hn = []
        for j in range(5):
            t = c[j] + d1[j] * h + d2[j] * H2 / 2
            lo = max(dn(t), float(a[13 + 2 * j])); hi = min(up(t), float(a[14 + 2 * j]))
            if not lo <= hi: raise SystemExit("empty hull intersection")
            hn.append((lo, hi)); row += [repr(lo), repr(hi)]
            ok_hull &= bool(float(a[13 + 2 * j]) <= lo and hi <= float(a[14 + 2 * j]))
            ok_c_in_h &= bool(arb(lo) <= c[j].lower() and c[j].upper() <= arb(hi))
            whull[j] = max(whull[j], hi - lo); whull_old[j] = max(whull_old[j], wid(ho[j]))
        out.write(" ".join(row) + "\n")
        x = float((TT - arb(tc)).mid())
        for j in range(5):
            wc[j] = max(wc[j], wid(c[j])); wprof[j] = max(wprof[j], wid(cp[j]))
            hj = [I(*v) for v in hn]
            pw = arb(wid(c[j])) + H * arb(wid(d1[j])) + H * H / 2 * arb(wid(f2(hj)[j]))
            wpt[j] = max(wpt[j], float(pw.upper()))
            for s in SEG:
                if s[0] <= x < s[1]:
                    wseg[s]["cond"][j] = max(wseg[s]["cond"][j], float(pw.upper()))
                    wseg[s]["prof"][j] = max(wseg[s]["prof"][j], wid(cp[j]))
        last = (ta, tb, tc, c, [I(*v) for v in hn])
        res["rows"] += 1
    out.write("END\n"); out.close()
    res.update(grid_identical_to_prof_tight=bool(ok_grid), box_hulls_identical_to_prof_tight=bool(ok_boxhull),
               centre_inside_prof_tight=bool(ok_in),
               hull_inside_old_hull=bool(ok_hull), centre_inside_new_hull=bool(ok_c_in_h))
    res["max_centre_width"] = dict(zip(NAMES, wc)); res["max_pointwise_width"] = dict(zip(NAMES, wpt))
    res["prof_tight_max_centre_width"] = dict(zip(NAMES, wprof))
    res["max_hull_width_new"] = dict(zip(NAMES, whull)); res["max_hull_width_old"] = dict(zip(NAMES, whull_old))
    res["by_x_segment"] = {f"x in [{s[0]}, {s[1]})": {"cond_pointwise": dict(zip(NAMES, v["cond"])),
                                                       "prof_tight_centre": dict(zip(NAMES, v["prof"]))} for s, v in wseg.items()}
    print(f"1. {res['rows']} rows; grid = prof_tight: {ok_grid}; box hulls = prof_tight hulls: {ok_boxhull}; centre in prof_tight centre: {ok_in}; "
          f"new hull in old hull: {ok_hull}; centre in new hull: {ok_c_in_h}")
    print("2. max centre width (all rows)   ", {n: f"{w:.2e}" for n, w in zip(NAMES, wc)})
    print("   max pointwise width (all rows)", {n: f"{w:.2e}" for n, w in zip(NAMES, wpt)})
    print("   prof_tight max centre width   ", {n: f"{w:.2e}" for n, w in zip(NAMES, wprof)})
    for s, v in wseg.items():
        print(f"   x = tau_tip - tau in [{s[0]}, {s[1]}): pointwise", {n: f"{w:.1e}" for n, w in zip(NAMES, v["cond"])},
              " prof_tight", {n: f"{w:.1e}" for n, w in zip(NAMES, v["prof"])})
    print("   hull width max new / old:", {n: f"{a_:.1e}/{b_:.1e}" for n, a_, b_ in zip(NAMES, whull, whull_old)})
    # ------------------------------------------------------------ 3. tip consistency
    ta, tb, tc, c, hl = last
    h = TT - arb(tc)
    assert arb(ta) <= TT.lower() and TT.upper() <= arb(tb)
    tr = [c[j] + f1(c)[j] * h + f2(hl)[j] * (h * h) / 2 for j in range(5)]
    tipst = [I(*tip["tip_" + n]) for n in NAMES]
    meet = all(max(tr[j].lower(), tipst[j].lower()) <= min(tr[j].upper(), tipst[j].upper()) for j in range(5))
    zeros = bool(tr[1].contains(0) and tr[3].contains(0) and tipst[1].contains(0) and tipst[3].contains(0))
    old_tt = I(10.249077650575472, 10.249077651519443); D_old = I(2.50080627654, 2.50080627668); k0_old = I(1.53982063762, 1.53982063822)
    inside_old = bool(old_tt.contains(TT) and D_old.contains(-tipst[0]) and k0_old.contains(tipst[2]))
    res["tip"] = dict(tau_tip=tip["tau_tip"], sigma_tip=[dn(8 - TT), up(8 - TT)], D=[dn(-tipst[0]), up(-tipst[0])],
                      k0=[dn(tipst[2]), up(tipst[2])], theta_tip=list(tip["tip_th"]),
                      transported_meets_tip_state=bool(meet), q_m_zero_contained=zeros, inside_s1_tight_intervals=inside_old,
                      widths=dict(tau_tip=wid(TT), **{n: wid(t) for n, t in zip(NAMES, tipst)}))
    print(f"3. tau_tip in [{tip['tau_tip'][0]!r}, {tip['tau_tip'][1]!r}] (width {wid(TT):.1e}); last row -> tau_tip meets the "
          f"tip state: {meet}; q = m = 0 contained: {zeros}; inside the s1_tight intervals: {inside_old}")
    # ------------------------------------------------------------ 3b. eps of the zeros of M in Y (rigorous by-product)
    # y* = xm + C r* + e, r*_0 = eps* - ec exactly (gen_input2.build); r*_0 = -(M(xm) + sum_{j>=1} A_j r*_j + DM e)/A_0
    vv = [float(t) for t in open(os.path.join(TIGHT, "in_br2.txt")).read().split()]
    r0b = [I(vv[36 + 2 * j], vv[37 + 2 * j]) for j in range(5)]
    br = json.load(open(os.path.join(TIGHT, "bracket2.json")))
    ea, eb = arb(repr(br["eps_a"])), arb(repr(br["eps_b"]))
    ec = arb(float(((ea + eb) / 2).mid()))
    Aj = [I(*tip[f"A{j}"]) for j in range(5)]
    d0 = -(I(*tip["M_xm"]) + sum((Aj[j] * r0b[j] for j in range(1, 5)), arb(0)) + I(*tip["DMe"])) / Aj[0]
    eps_cond = ec + d0
    eps_lo = max(dn(eps_cond), float(ea.lower())); eps_hi = min(up(eps_cond), float(eb.upper()))
    res["eps_star_conditioned"] = dict(interval=[eps_lo, eps_hi], width=eps_hi - eps_lo,
                                       bracket2_width=float((eb - ea).mid()))
    print(f"   eps* (zeros of M in Y) in [{eps_lo!r}, {eps_hi!r}], width {eps_hi - eps_lo:.2e} "
          f"(bracket2: {float((eb - ea).mid()):.2e}); its width is set by the d box, |A_d| rn1/|A_0|")
    # ------------------------------------------------------------ 4. float sanity check (not rigorous)
    import numpy as np
    from scipy.integrate import solve_ivp
    v = [float(t) for t in open(os.path.join(TIGHT, "in_br2.txt")).read().split()]
    xm = np.array([float(arb(v[1 + 2 * i]).union(arb(v[2 + 2 * i])).mid()) for i in range(5)])
    Cm = np.array(v[11:36]).reshape(5, 5)
    M0 = 0.5 * sum(tip["M_xm"]); A0 = 0.5 * sum(tip["A0"])
    y0 = xm + Cm[:, 0] * (-M0 / A0)
    def fv(t, y):
        p, q, k, m, th = y
        return [k * q, -1 - k * p, -m, -p / 4, -k]
    rows = [l.split() for l in open(os.path.join(HERE, "prof_tip.txt")) if not l.startswith("END")]
    pick = [0, 5000, 15000, 25000, 35000, 45000, 50000, len(rows) - 1]
    ts = sorted(float(rows[i][2]) for i in pick)
    sol = solve_ivp(fv, [0, ts[-1]], y0, method="DOP853", rtol=2.3e-14, atol=1e-22, t_eval=ts)
    worst = 0.0
    for i, t in zip(sorted(pick), sol.t):
        r = rows[i]
        for j in range(5):
            lo, hi = float(r[3 + 2 * j]), float(r[4 + 2 * j]); mid, rad = (lo + hi) / 2, (hi - lo) / 2
            dev = abs(sol.y[j][list(sol.t).index(t)] - mid)
            worst = max(worst, dev / max(rad, 1e-300) if rad > 0 else 0.0)
    res["float_check_worst_offset_in_radii"] = worst
    print(f"4. float trajectory (DOP853) from a zero of M in Y: worst |float - centre|/radius = {worst:.2f} over 8 rows x 5 components"
          " (float error ~1e-13 near the tip; not part of the certificate)")
    res["all_rigorous_checks_pass"] = bool(ok_grid and ok_boxhull and ok_in and ok_hull and ok_c_in_h and meet and zeros and inside_old)
    res["seconds"] = time.time() - t0
    json.dump(res, open(os.path.join(HERE, "ext_tip_check.json"), "w"), indent=1)
    print(f"ALL RIGOROUS CHECKS PASS: {res['all_rigorous_checks_pass']}  [{res['seconds']:.0f}s]")
    return 0 if res["all_rigorous_checks_pass"] else 1


if __name__ == "__main__":
    sys.exit(main())
