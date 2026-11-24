"""T7 far constants on the (SE)-corrected base (STEP4_RESULT §2(ii), p2/se_corr.py).

The corrected base is k^ = k~ + dk, dk = (a1/ell) B(u) + (a2/ell^2) B'(u) on u = x/ell in [0.6, 0.8], so it equals the
polished base for x < 0.6 ell(tau).  For x >= x_c := 0.6 ell_lo (ell_lo = min of ell~ over the slab) every base quantity
q^ lies in q~ + [-w_q, w_q] with the widths below.  Since review B1 the inputs are direct: a1, a2, a1', a2' are the
certified enclosures of p2/t6c/slabc.py (T6's engine run on the corrected base) and the residual widths dz come from the
corrected and polished T6 jets (p2/t6c/jets, t6/jets); p2/out/se_corr.json is no longer read (--source se_corr keeps it):
  k^(m): d_m = a1 Bs_m / L^(m+1) + a2 Bs_(m+1) / L^(m+2)        (Bs_m = sup |B^(m)| in u, L = ell_lo)
  Q = int k'^2:  DQ = 2 K1 dl_1 + d_1 dl_1                      (dl_m = int_window |dk^(m)| dx)
  theta: Dth = a1 + a2 Bs_0 / L;   p, q: dp = (2 Lh + 3) Dth;   n, t1: dn = Dth
  r, A: dz = max(sup|r^| + sup|r~|, sup|A^| + sup|A~|) over the far cells (direct T6 hulls);  V* = -k'' + c p + b n: d_2 + |c| dp + |b| dn;  V = V* + r: + dz
  varpi = -k k' + Q + c (q - x) + b t1 + A:  6 (K0 d1 + K1 d0 + d0 d1) + DQ + |c| dp + |b| dn + dz
  theta_tau at fixed x: dtheta = a1 I_B(u) + (a2/ell) B(u), so
     |d_tau dtheta| <= |a1'| + |a1| Bs_0 0.8 |ell_t| / L + (|a2'| / L + |a2| |ell_t| / L^2) Bs_0 + |a2| Bs_1 0.8 |ell_t| / L^2.
far7.py and nl7.py (--se) add these balls to the base jets on every cell that meets x >= x_c; the farg Gram matrices get
a certified increment (farg_se below).  All widths are rounded up.
  python se7.py widths [--p2 ../p2]                     -> se_widths.json (portable cache)
  python se7.py farg --in run_v2/out_farg --out run_v2se/out_farg"""
import os, sys, json, math, argparse
HERE = os.path.dirname(os.path.abspath(__file__))
STEP4 = os.path.dirname(HERE)
CACHE = os.path.join(HERE, 'se_widths.json')
_W = None


def ub(x):
    from flint import arb
    return float(arb(x).upper()) * (1 + 1e-12) + 1e-300


def sha256(path):
    import hashlib
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for b in iter(lambda: f.read(1 << 20), b''): h.update(b)
    return h.hexdigest()


def sha256_files(paths):
    """sha256 of the concatenated per-file sha256 (sorted paths): one fingerprint for a directory of jets"""
    import hashlib
    h = hashlib.sha256()
    for p in sorted(paths): h.update((os.path.basename(p) + sha256(p)).encode())
    return h.hexdigest()


def build(p2=None, source='t6c'):
    """per-slab widths of the (SE)-corrected base against the polished one, on x >= x_c = 0.6 ell_lo.
    source 't6c' (review B1, default): the correction parameters a1, a2, a1', a2' are the certified enclosures of
    p2/t6c/slabc.py (T6's engine on the corrected base; the 'sharp' subdivided sups, max over the slab's sub-slabs), and
    the residual widths are direct: |r^ - r~| <= sup|r^| + sup|r~|, |A^ - A~| <= sup|A^| + sup|A~| over the far cells,
    r^, A^ from p2/t6c/jets and r~, A~ from t6/jets.  source 'se_corr' (superseded): p2/out/se_corr.json."""
    import numpy as np
    from flint import arb, ctx
    ctx.prec = 128
    p2 = os.path.normpath(p2 or os.path.join(STEP4, 'p2'))
    sys.path.insert(0, p2); sys.path.insert(0, os.path.join(STEP4, 't6'))
    import se_corr as SC, reader as RD
    t6res = os.path.join(STEP4, 't6', 'residuals.json'); t6jets = os.path.join(STEP4, 't6', 'jets')
    T6 = {r['slab']: r for r in json.load(open(t6res))['per_slab']}
    if source == 't6c':
        cres = os.path.join(p2, 't6c', 'residuals.json'); cjets = os.path.join(p2, 't6c', 'jets')
        R6C = json.load(open(cres)); assert R6C.get('_ranges_widened'), 'p2/t6c/residuals.json: ranges not widened'
        rows = {}
        for r in R6C['per_slab']:
            ents = r['se_direct']
            rows[r['slab']] = {k: max(e['sharp'][k] for e in ents) for k in ('a1', 'a2', 'a1p', 'a2p')}
        srcs = dict(t6c_residuals=dict(path=os.path.relpath(cres, STEP4), sha256=sha256(cres)),
                    t6c_jets=dict(path=os.path.relpath(cjets, STEP4), sha256_of_files=sha256_files(
                        [os.path.join(cjets, f'jets_{i:03d}.npz') for i in rows])))
    else:
        f = os.path.join(p2, 'out', 'se_corr.json')
        rows = {r['slab']: r for r in json.load(open(f))}
        srcs = dict(se_corr=dict(path=os.path.relpath(f, STEP4), sha256=sha256(f)))
    psup, pl1 = SC.bump_norms()
    Ip = arb(math.factorial(6) ** 2) / math.factorial(13); Wu = arb(SC.WU)
    Bs = [psup[k] / (Wu ** (k + 1) * Ip) for k in range(10)]; Bl = [pl1[k] / (Wu ** k * Ip) for k in range(10)]
    out = {}
    for i, row in sorted(rows.items()):
        sc = T6[i]['scalars']; J = RD.jets(i, d=t6jets)
        ell_lo, ell_hi = sc['range_ell']; L, Lh = arb(ell_lo), arb(ell_hi)
        ellt = arb(max(abs(v) for v in sc['range_ell_tau']))
        cb = arb(max(abs(v) for v in sc['range_c'])); bb = arb(max(abs(v) for v in sc['range_b']))
        x = J['x']; win = (x[..., 1] >= SC.UA * ell_lo) & (x[..., 0] <= SC.UB * ell_hi)
        K = [arb(float(np.max(np.abs(J[f'k{j}'][win])))) for j in range(3)]
        a1, a2, a1p, a2p = (arb(row[k]) for k in ('a1', 'a2', 'a1p', 'a2p'))
        xc = float((arb('0.6') * L).lower())
        if source == 't6c':
            JC = RD.jets(i, d=cjets)
            fp = J['x'][..., 1] >= xc; fc = JC['x'][..., 1] >= xc
            sup = lambda JJ, m, k: float(np.max(np.abs(JJ[k][m]))) if np.any(m) else 0.0
            dz = arb(max(sup(JC, fc, 'r') + sup(J, fp, 'r'), sup(JC, fc, 'A') + sup(J, fp, 'A')))
        else:
            dz = arb(row['dz'])
        d = [a1 * Bs[k] / L ** (k + 1) + a2 * Bs[k + 1] / L ** (k + 2) for k in range(7)]
        dl = [a1 * Bl[k] / L ** k + a2 * Bl[k + 1] / L ** (k + 1) for k in range(3)]
        DQ = 2 * K[1] * dl[1] + d[1] * dl[1]
        Dth = a1 + a2 * Bs[0] / L; dp = (2 * Lh + 3) * Dth; dn = Dth
        kd = K[0] * d[1] + K[1] * d[0] + d[0] * d[1]
        Vst = d[2] + cb * dp + bb * dn
        tht = a1p + a1 * Bs[0] * arb('0.8') * ellt / L + (a2p / L + a2 * ellt / L ** 2) * Bs[0] + a2 * Bs[1] * arb('0.8') * ellt / L ** 2
        out[str(i)] = dict(xc=xc, k=[ub(v) for v in d], Q=ub(DQ), theta=ub(Dth), p=ub(dp), q=ub(dp),
                           n=ub(dn), t1=ub(dn), r=ub(dz), A=ub(dz), Vstar=ub(Vst), V=ub(Vst + dz),
                           varpi=ub(6 * kd + DQ + cb * dp + bb * dn + dz), tht=ub(tht))
    if source == 't6c':                       # slabs outside p2/t6c (0-1, outside P2: SLAB0 = 2) keep the se_corr widths, flagged
        f = os.path.join(p2, 'out', 'se_corr.json')
        for r in json.load(open(f)):
            i = r['slab']
            if i in rows or str(i) in out: continue
            sc = T6[i]['scalars']; L = arb(sc['range_ell'][0]); Lh = arb(sc['range_ell'][1])
            o = json.load(open(os.path.join(HERE, 'se_widths_secorr.json')))['per_slab'].get(str(i))
            if o is not None: out[str(i)] = dict(o, source='se_corr (slab outside P2 and outside p2/t6c)')
        srcs['outside_P2_se_corr'] = dict(path=os.path.relpath(f, STEP4), sha256=sha256(f), slabs=sorted(int(k) for k, v in out.items() if 'source' in v),
                                          widths=dict(path='t7/se_widths_secorr.json', sha256=sha256(os.path.join(HERE, 'se_widths_secorr.json'))))
    args = dict(source=source, inputs=dict(srcs, t6_residuals=dict(path=os.path.relpath(t6res, STEP4), sha256=sha256(t6res)),
                                           t6_jets=dict(path=os.path.relpath(t6jets, STEP4), sha256_of_files=sha256_files(
                                               [os.path.join(t6jets, f'jets_{i:03d}.npz') for i in rows]))),
                code=dict(se7=sha256(os.path.abspath(__file__)), se_corr_bump_norms=sha256(os.path.join(p2, 'se_corr.py'))),
                note='only se_corr.bump_norms (exact bump sups/L1 norms) is used from se_corr.py when source = t6c')
    json.dump(dict(note=__doc__.split('\n')[0], args=args, per_slab=out), open(CACHE, 'w'), indent=1)
    return out


def widths(slab):
    """the widths of slab `slab` (dict, floats), from the cache se_widths.json"""
    global _W
    if _W is None:
        if not os.path.exists(CACHE): raise SystemExit('se_widths.json missing: run  python se7.py widths  first')
        _W = json.load(open(CACHE))['per_slab']
    return _W[str(slab)]


def farg_se(fg, slab):
    """certified increment of the farg Gram traces for the corrected base: on x >= x_c (zeta = 1 there)
       Gamma^_ii <= Gamma_ii + int_{x_c}^{ell} Wbar (2 |g_i^(m)| Dg_i^(m) + Dg_i^(m)^2) / ||chi g_i||^2,  m = 1 (Gamma0), 3 (Gamma2),
    with sup |g^(m)| from T6's v2 hulls on the far cells (Frenet recursion), Dg from the widths, and
    int Wbar <= 1/(3 x_c^3) + (1 + ell_lo^2/4)^-2 (ell_hi - x_c).  Gamma2E (rho, [0, 34.25]) is unchanged."""
    import numpy as np
    from flint import arb, ctx
    ctx.prec = 128
    sys.path.insert(0, os.path.join(STEP4, 't6')); import reader as RD
    w = widths(slab); sc = RD.residuals(os.path.join(STEP4, 't6', 'residuals.json'))['_by'][slab]['scalars']
    J = RD.jets(slab, d=os.path.join(STEP4, 't6', 'jets')); far = J['x'][..., 1] >= w['xc']
    S = lambda nm: arb(float(np.max(np.abs(J[nm][far])))) if np.any(far) else arb(0)
    K = [S(f'k{j}') for j in range(6)]; P0, Q0, N0, T0 = S('p'), S('q'), S('n'), S('t1')
    cb = arb(max(abs(v) for v in sc['range_c'])); bb = arb(max(abs(v) for v in sc['range_b']))
    d = [arb(v) for v in w['k']]; dp, dn = arb(w['p']), arb(w['n'])
    # sups of the uncorrected jets (Frenet: p' = -k q, q' = 1 + k p, n' = -k t1, t1' = k n)
    P1 = K[0] * Q0; Q1 = 1 + K[0] * P0; P2 = K[1] * Q0 + K[0] * Q1; Q2 = K[1] * P0 + K[0] * P1; P3 = K[2] * Q0 + 2 * K[1] * Q1 + K[0] * Q2
    N1 = K[0] * T0; T1 = K[0] * N0; N2 = K[1] * T0 + K[0] * T1; T2 = K[1] * N0 + K[0] * N1; N3 = K[2] * T0 + 2 * K[1] * T1 + K[0] * T2
    # their changes: delta(ab) <= |da| |b| + |a| |db| + |da| |db|
    pr = lambda da, b, a, db: da * b + a * db + da * db
    Dp1 = pr(d[0], Q0, K[0], dp); Dq1 = pr(d[0], P0, K[0], dp)
    Dp2 = pr(d[1], Q0, K[1], dp) + pr(d[0], Q1, K[0], Dq1); Dq2 = pr(d[1], P0, K[1], dp) + pr(d[0], P1, K[0], Dp1)
    Dp3 = pr(d[2], Q0, K[2], dp) + 2 * pr(d[1], Q1, K[1], Dq1) + pr(d[0], Q2, K[0], Dq2)
    Dn1 = pr(d[0], T0, K[0], dn); Dt1 = pr(d[0], N0, K[0], dn)
    Dn2 = pr(d[1], T0, K[1], dn) + pr(d[0], T1, K[0], Dt1); Dt2 = pr(d[1], N0, K[1], dn) + pr(d[0], N1, K[0], Dn1)
    Dn3 = pr(d[2], T0, K[2], dn) + 2 * pr(d[1], T1, K[1], Dt1) + pr(d[0], T2, K[0], Dt2)
    V1 = K[3] + cb * P1 + bb * N1; V3 = K[5] + cb * P3 + bb * N3
    DV1 = d[3] + cb * Dp1 + bb * Dn1; DV3 = d[5] + cb * Dp3 + bb * Dn3
    G1, G3, DG1, DG3 = (P1, N1, V1), (P3, N3, V3), (Dp1, Dn1, DV1), (Dp3, Dn3, DV3)
    xc = arb(w['xc']); ell_lo, ell_hi = (arb(v) for v in sc['range_ell'])
    IW = 1 / (3 * xc ** 3) + (1 + ell_lo * ell_lo / 4) ** -2 * (ell_hi - xc)
    nrm = [arb(v) for v in fg['norm_chi_g_min']]
    inc0 = [IW * (2 * G1[i] * DG1[i] + DG1[i] ** 2) / nrm[i] ** 2 for i in range(3)]
    inc2 = [IW * (2 * G3[i] * DG3[i] + DG3[i] ** 2) / nrm[i] ** 2 for i in range(3)]
    out = dict(fg)
    sq_up = lambda a, b: ub((arb(a) ** 2 + b).sqrt())
    out['sqrt_lmax_G0'] = sq_up(fg['sqrt_lmax_G0'], sum(inc0, arb(0))); out['sqrt_lmax_G2'] = sq_up(fg['sqrt_lmax_G2'], sum(inc2, arb(0)))
    out['sqrt_diag_G0'] = [sq_up(v, inc0[i]) for i, v in enumerate(fg['sqrt_diag_G0'])]
    out['sqrt_diag_G2'] = [sq_up(v, inc2[i]) for i, v in enumerate(fg['sqrt_diag_G2'])]
    out['se_correction'] = dict(inc_G0=[ub(v) for v in inc0], inc_G2=[ub(v) for v in inc2], xc=w['xc'],
                                note='Gamma0/Gamma2 traces of the (SE)-corrected base: polished-base values + certified far increments (se7.farg_se)')
    return out


if __name__ == '__main__':
    ap = argparse.ArgumentParser(); ap.add_argument('what', choices=['widths', 'farg'])
    ap.add_argument('--p2', default=None); ap.add_argument('--in', dest='inp', default=os.path.join('run_v2', 'out_farg'))
    ap.add_argument('--source', default='t6c', choices=['t6c', 'se_corr'])
    ap.add_argument('--out', default=os.path.join('run_v2se', 'out_farg'))
    A = ap.parse_args()
    if A.what == 'widths':
        o = build(A.p2, A.source)
        print(f'{len(o)} slabs; max widths:', {k: f"{max(max(v[k]) if isinstance(v[k], list) else v[k] for v in o.values()):.3e}"
                                               for k in ('k', 'Q', 'theta', 'p', 'n', 'V', 'varpi', 'tht')})
    else:
        src = A.inp if os.path.isabs(A.inp) else os.path.join(HERE, A.inp); dst = A.out if os.path.isabs(A.out) else os.path.join(HERE, A.out)
        os.makedirs(dst, exist_ok=True); n = 0; worst = 0.0
        for f in sorted(os.listdir(src)):
            if not (f.startswith('farg_') and f.endswith('.json')): continue
            fg = json.load(open(os.path.join(src, f))); o = farg_se(fg, fg['slab'])
            worst = max(worst, max(o['se_correction']['inc_G0'] + o['se_correction']['inc_G2']))
            json.dump(o, open(os.path.join(dst, f), 'w'), indent=1); n += 1
        print(f'{n} farg files corrected -> {dst}; max increment of a normalised Gram diagonal {worst:.3e}')
