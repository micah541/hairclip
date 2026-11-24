"""
Rigorous winding numbers of the Evans functions D_even = W13(tip), D_odd = W02(tip) along the upper half of the
contour  1.85 -> 1.85+0.3i -> -0.125+4.1i -> -0.125  (D(conj lam) = conj D(lam), so the full winding is
(arg change along the upper half)/pi).

For each segment [lam_a, lam_b] the CAPD program evans_seg returns an affine enclosure
      D(lam_a + t (lam_b - lam_a))  in  X + V (t - 1/2) + E,   t in [0, 1]      (E a complex box),
we check 0 notin S = X + V[-1/2,1/2] + E (convex), and pick alpha_i in S_{i-1} cap S_i (intersection of the end box
of segment i-1 and the start box of segment i, both of which contain D(lam_i)).  Then
      sum_i pv arg(alpha_{i+1}/alpha_i) = arg change of D along the path   (exactly; each S_i is convex, 0 notin S_i).
At the two real endpoints alpha is taken real (D is real there).
"""
import sys, os, json, math, time, subprocess, cmath
from multiprocessing import Pool
sys.path.insert(0, '/home/claude/s1spec/lin'); sys.path.insert(0, '/home/claude/s1_cap')
from flint import arb, acb
import lintail2 as LT
import seg2

PATH = [complex(1.85, 0.0), complex(1.85, 0.3), complex(-0.125, 4.1), complex(-0.125, 0.0)]
CELL = 0.05
ORDER = "18"
CKPT = "contour_ckpt.json"
ETA_FILE = "eta_cache.json"


def cell_of(z):
    return (math.floor(z.real / CELL), math.floor(z.imag / CELL))


def cell_box(c, pad=0.012):
    i, j = c
    return (i * CELL - pad, (i + 1) * CELL + pad, j * CELL - pad, (j + 1) * CELL + pad)


def eta_for_cell(c):
    box = cell_box(c)
    lam = LT.cbox(*box)
    lc = complex((box[0] + box[1]) / 2, (box[2] + box[3]) / 2)
    eta, Lip, gmin = LT.eta_bound(0.25 - lam, 0.25 - lc)
    return c, float(eta.upper()), float(Lip.upper()), gmin


def parse(out):
    aff = {}
    for l in out.splitlines():
        if l.startswith("AFF"):
            p = l.split()
            aff[p[1]] = tuple(float(v) for v in p[2:])
    return aff


def class_set(aff, re_name, im_name):
    """returns (Xm, Vm, Er, Ei) with Xm, Vm complex doubles and S subset Xm + Vm s + ([-Er,Er] + i[-Ei,Ei])"""
    xr0, xr1, vr0, vr1, er0, er1 = aff[re_name]
    xi0, xi1, vi0, vi1, ei0, ei1 = aff[im_name]
    Xm = complex((xr0 + xr1) / 2, (xi0 + xi1) / 2)
    Vm = complex((vr0 + vr1) / 2, (vi0 + vi1) / 2)
    # absorb widths of X and V (|s| <= 1/2), outward
    Er = max(abs(er0), abs(er1)) + (xr1 - xr0) / 2 + (vr1 - vr0) / 4
    Ei = max(abs(ei0), abs(ei1)) + (xi1 - xi0) / 2 + (vi1 - vi0) / 4
    Er *= 1 + 1e-9; Ei *= 1 + 1e-9
    return Xm, Vm, Er, Ei


def excludes_zero(Xm, Vm, Er, Ei):
    """rigorous: min_{|s|<=1/2} |Xm + Vm s| > |E| (sufficient for 0 notin S); returns margin ratio"""
    X = acb(Xm.real, Xm.imag); V = acb(Vm.real, Vm.imag)
    nV2 = abs(V) ** 2
    cand = [arb(-0.5), arb(0.5)]
    if Vm != 0:
        sstar = -(X.real * V.real + X.imag * V.imag) / nV2
        if sstar > -0.5 and sstar < 0.5:
            cand.append(sstar)
        elif not (sstar < -0.5 or sstar > 0.5):
            cand.append(sstar)   # straddles: include (min over a superset of candidates is still a lower bound)
    # min over s of |X+Vs|^2 is attained at s* (if inside) or an endpoint; evaluate all candidates and take min lower
    d2 = None
    for s in cand:
        val = abs(X + V * acb(s)) ** 2
        d2 = val if d2 is None else (val if val.lower() < d2.lower() else d2)
    E2 = arb(Er) ** 2 + arb(Ei) ** 2
    ok = bool(d2 > E2)
    ratio = float((E2 / d2).upper()) ** 0.5 if d2 > 0 else float('inf')
    return ok, ratio


def run_segment(args):
    la, lb, eta = args
    path = f"/tmp/seg_{os.getpid()}.txt"
    try:
        seg2.make(la, lb, path, eta=arb(eta))
        out = subprocess.run(["./evans_seg", path, ORDER], capture_output=True, text=True).stdout
    except Exception as ex:
        return dict(la=la, lb=lb, ok=False, err=repr(ex))
    if "RESULT OK" not in out:
        return dict(la=la, lb=lb, ok=False, err=out[-400:])
    aff = parse(out)
    if len(aff) < 4:
        return dict(la=la, lb=lb, ok=False, err="no AFF")
    rec = dict(la=la, lb=lb, ok=True)
    for cls, (rn, inn) in (("even", ("W13r", "W13i")), ("odd", ("W02r", "W02i"))):
        Xm, Vm, Er, Ei = class_set(aff, rn, inn)
        ok, ratio = excludes_zero(Xm, Vm, Er, Ei)
        rec[cls] = dict(X=[Xm.real, Xm.imag], V=[Vm.real, Vm.imag], E=[Er, Ei], ok=ok, ratio=ratio)
        rec["ok"] = rec["ok"] and ok
    return rec


def Lrule(eta):
    return 0.004 if eta < 0.03 else (0.002 if eta < 0.05 else (0.001 if eta < 0.07 else 0.0005))


def segments_of_path(cache):
    segs = []
    for a, b in zip(PATH[:-1], PATH[1:]):
        Ltot = abs(b - a); d = (b - a) / Ltot
        s0 = 0.0
        while s0 < Ltot - 1e-12:
            z = a + d * s0
            L = Lrule(cache[cell_of(z)]["eta"]) if cell_of(z) in cache else 0.004
            s1 = min(s0 + L, Ltot)
            if Ltot - s1 < 1e-9:
                s1 = Ltot
            segs.append((a + d * s0 if s0 > 0 else a, a + d * s1 if s1 < Ltot else b))
            s0 = s1
    return segs


def cz(z):
    return [z.real, z.imag]


def main(nproc=2, out_json="contour_result.json"):
    t0 = time.time()
    cache = {}
    if os.path.exists(ETA_FILE):
        cache = {tuple(json.loads(k)): v for k, v in json.load(open(ETA_FILE)).items()}
    # cells needed: cover the path finely
    need = set()
    for a, b in zip(PATH[:-1], PATH[1:]):
        for i in range(4001):
            need.add(cell_of(a + (b - a) * i / 4000))
    need = [c for c in need if c not in cache]
    if need:
        with Pool(nproc) as P:
            for c, eta, Lip, g in P.imap_unordered(eta_for_cell, need):
                cache[c] = dict(eta=eta, Lip=Lip, gmin=g)
        json.dump({json.dumps(list(k)): v for k, v in cache.items()}, open(ETA_FILE, "w"))
    print(f"eta cache: {len(cache)} cells, max eta = {max(v['eta'] for v in cache.values()):.4g}, "
          f"max Lip = {max(v['Lip'] for v in cache.values()):.4g} ({time.time()-t0:.0f}s)", flush=True)

    def eta_seg(a, b):
        cs = {cell_of(a), cell_of(b), cell_of((a + b) / 2)}
        assert abs(b - a) < 0.012
        return max(cache[c]["eta"] for c in cs)

    done = []
    todo = segments_of_path(cache)
    if os.path.exists(CKPT):
        ck = json.load(open(CKPT))
        done = ck["done"]
        todo = [(complex(*a), complex(*b)) for a, b in ck["todo"]]
        for r in done:
            r["la"] = complex(*r["la"]); r["lb"] = complex(*r["lb"])
        print(f"resumed: {len(done)} done, {len(todo)} to do", flush=True)
    print(f"initial segments: {len(todo)}", flush=True)
    rounds = 0
    while todo:
        rounds += 1
        jobs = [(a, b, eta_seg(a, b)) for a, b in todo]
        newtodo = []
        k = 0
        with Pool(nproc) as P:
            for r in P.imap_unordered(run_segment, jobs, chunksize=4):
                k += 1
                if r["ok"]:
                    done.append(r)
                else:
                    a, b = r["la"], r["lb"]
                    if abs(b - a) < 1e-5:
                        raise RuntimeError(f"segment {a}->{b} cannot be resolved: {r.get('err', r)}")
                    m = (a + b) / 2
                    newtodo += [(a, m), (m, b)]
                if k % 200 == 0:
                    print(f"  round {rounds}: {k}/{len(jobs)} ({time.time()-t0:.0f}s)", flush=True)
                    json.dump(dict(done=[dict(r_, la=cz(r_["la"]), lb=cz(r_["lb"])) for r_ in done],
                                   todo=[[cz(a_), cz(b_)] for a_, b_ in newtodo] +
                                        [[cz(a_), cz(b_)] for (a_, b_, e_) in jobs[k:]]), open(CKPT, "w"))
        todo = newtodo
        json.dump(dict(done=[dict(r_, la=cz(r_["la"]), lb=cz(r_["lb"])) for r_ in done],
                       todo=[[cz(a_), cz(b_)] for a_, b_ in todo]), open(CKPT, "w"))
        print(f"round {rounds}: {len(done)} segments certified, {len(todo)} to refine ({time.time()-t0:.0f}s)", flush=True)

    def key(r):
        a = r["la"]
        for e, (p, q) in enumerate(zip(PATH[:-1], PATH[1:])):
            d = q - p
            t = ((a - p) * d.conjugate()).real / abs(d) ** 2
            if -1e-12 <= t <= 1 + 1e-12 and abs((a - p) - t * d) < 1e-9:
                return (e, t)
        raise RuntimeError("point not on path")
    done.sort(key=key)
    for r1, r2 in zip(done[:-1], done[1:]):
        assert abs(r1["lb"] - r2["la"]) < 1e-12, (r1["lb"], r2["la"])
    assert abs(done[0]["la"] - PATH[0]) < 1e-12 and abs(done[-1]["lb"] - PATH[-1]) < 1e-12
    json.dump(dict(segments=[dict(r_, la=cz(r_["la"]), lb=cz(r_["lb"])) for r_ in done]), open(out_json, "w"))
    import winding_cert
    report = {}
    for cls in ("even", "odd"):
        report[cls] = winding_cert.winding_from(done, cls)    # argument cones (per-segment positive scales)
    print(json.dumps(report, indent=1), flush=True)
    json.dump(report, open("contour_report.json", "w"), indent=1)
    return report


def box_of(S, s):
    """complex box (re_lo, re_hi, im_lo, im_hi) = X + V s + E"""
    X = complex(*S["X"]); V = complex(*S["V"]); Er, Ei = S["E"]
    c = X + V * s
    return (c.real - Er, c.real + Er, c.imag - Ei, c.imag + Ei)


def inter(b1, b2):
    lo_r, hi_r = max(b1[0], b2[0]), min(b1[1], b2[1])
    lo_i, hi_i = max(b1[2], b2[2]), min(b1[3], b2[3])
    if lo_r > hi_r or lo_i > hi_i:
        return None
    return (lo_r, hi_r, lo_i, hi_i)


def winding_from(done, cls):
    alphas = []
    # start (real endpoint): real point of the start box of the first segment
    b0 = box_of(done[0][cls], -0.5)
    assert b0[2] <= 0 <= b0[3], "start box does not meet the real axis"
    assert b0[0] > 0 or b0[1] < 0, "sign at the start not determined"
    alphas.append(complex((b0[0] + b0[1]) / 2, 0.0))
    for r1, r2 in zip(done[:-1], done[1:]):
        I = inter(box_of(r1[cls], 0.5), box_of(r2[cls], -0.5))
        assert I is not None, f"empty intersection at {r1['lb']}"
        alphas.append(complex((I[0] + I[1]) / 2, (I[2] + I[3]) / 2))
    bN = box_of(done[-1][cls], 0.5)
    assert bN[2] <= 0 <= bN[3], "end box does not meet the real axis"
    assert bN[0] > 0 or bN[1] < 0, "sign at the end not determined"
    alphas.append(complex((bN[0] + bN[1]) / 2, 0.0))
    total = 0.0; maxstep = 0.0
    for a1, a2 in zip(alphas[:-1], alphas[1:]):
        d = cmath.phase(a2 / a1)
        maxstep = max(maxstep, abs(d))
        total += d
    w = total / math.pi
    return dict(arg_change_over_pi=w, winding=round(w), max_step=maxstep, n=len(done),
                sign_start=1 if alphas[0].real > 0 else -1, sign_end=1 if alphas[-1].real > 0 else -1,
                worst_ratio=max(r[cls]["ratio"] for r in done))


if __name__ == "__main__":
    main()
