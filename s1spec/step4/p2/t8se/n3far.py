"""Lemma N3' (LEMMA_F.md 6.1): certified per-slab increments of the level-j constants, j = 0..3, from the nonlinear terms N
of the exact operator and from the modulation terms Q, with x-dependent majorants on the actual far boxes.

    python3 -I n3far.py --slabs 2-180 --out constants/n3.json
            [--levels 0-3] [--amp constants/amp_fin.json] [--consts constants/constants_fin.json] [--Y 1.7e-12]
            [--fac 1.7] [--b0 1.5e-3] [--b1 3e-4] [--b2 1e-2] [--etaF 1e-2] [--eta1 1e-5] [--noQ]
    then:  camp.py ... --n3 constants/n3.json        (about 30 s for slabs 2-180, one process);
           n3boxes.py constants/n3.json AMP_WITH_N3 CONSTS checks the assumed boxes (review B2).  With --n3 the top term
           2 sup At^+ is in dtheta at every level, so cconst's TOPLOSS may be set to 0 (keeping it counts the top term twice).
Needs n3forms.py (n3sym.py: N) and n3qforms.py (n3qsym.py: Q); n3wcheck.py checks the integration-by-parts identities.

Boxes (bootstrap hypotheses; the defaults hold for every Ybar <= Y = 1.7e-12, checked by n3boxes.py; review B2):
  core y >= 28:  |e| <= fac E(x), fac = 1.7,  E = eta1 + sqrt(F0) y^{3/2}/sqrt(84),  sqrt F0 <= aF0 Y + bF0 (from --amp, per slab)
                 (|e(28)| <= eta1, and H_28(y)^2 = int_28^y T0^-1 <= (y^3 - 28^3)/84);
                 |e'| <= b0 = 1.5e-3,  x|e''| <= b1 = 3e-4,  x^2|e'''| <= b2 = eta_F = 1e-2
                 (certified with N3' at Ybar = 1e-12: 5.43e-4, 7.6e-5, 5.9e-3; margins x2.76, x3.97, x1.69 and x1.70 for |e|).
  y < 28 (x <= 34.25 lies in T5's tip box):  |e|, |e'| <= eta1 = 1e-5 (T5 2.5 hypothesis; review A3),  x|e''|, x^2|e'''| <= eta_F
                 (the e'' part is the assembly's Lemma 4a).
--Y enters the F0 part of the |e| box and the modulation bound; b0, b1, b2 must be certified boxes at the same Ybar (every
smaller Ybar is then covered). --amp: only its level-0 split aF0, bF0 is read (the |e| box). Run once with the camp output
without N3', then once more with the camp output of that n3 file (a fixed point: the sqrt F0 ratio is 1.00000); n3boxes.py
asserts the |e| box (ratio <= fac) and the three derivative boxes for the final camp run.
Modulation (Q): |x| <= C_x Ybar + x_r with x_r = mbar |c_phi|_2 R2 per slab (T5 (M), review B6; C_x from the t7 block of --consts)
and |theta| = |x_3|/n <= (C_x Ybar + x_r)/nV with the per-slab lower bound nV of n = ||chi V~|| from --consts.
Split (LEMMA_F 6.1):  N + Q = (1 - th) N + th e'''' + th k^2 e'' - th c0 e + th P  [weak form, two derivatives onto T_j w_j]
                             + g e' - (th c + x1) e                                    [strong form],
P = (k''' - r') e e'/alpha, g = U/alpha + 5 th k k', U = (x1 + th c) q + (x2 + th b) t1, c0 = -3k'^2 - 4k k''.
Output per slab, per level j (key 'lev'): dtheta, dsigma, dk (to ADD to the theta_j loss, sigma_j and k_jj), d[m], k[m], c[m]
(new couplings d_jm F_m and k_jm V_m^2 in h_j, c_jm in g_j; m < j), hf (absolute forcing in h_j), and diagnostics.
camp.py --n3 FILE applies them: theta_j -= dtheta, lambda_j -= dsigma, k_jj += dk, k_jm += k[m], c_jm += c[m],
h_j += d[m] F_m + hf."""
import os, sys, json, math, time, argparse
from math import comb
HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path: sys.path.insert(0, HERE)
from flint import arb, ctx
ctx.prec = 64
import cconst
from cconst import U, L, iv, H, pos, absb, sqrtp, X0
import n3forms, n3qforms

# Young weights per level (any positive values are valid; D_j-side costs only lower theta_j, which must stay > 0, while F_j,
# N_{j+1}^2 and V_j^2 costs enter sigma_j and k_jj):
#   WD: core pairs with sqrt(D_j); THR: ramp (j, j+2) pair; EA, EB: a_flat with sqrt(D_j) and N_{j+1}; ECC, ECR: a_flat with
#   the T''-weighted w_j on core and ramp.  DE0: level 0, absolute e-terms; DG0: level 0, ramp T'-pair with N_1.
YW = {j: dict(WD=arb(1), THR=arb('1e-3'), EA=arb('1e-3'), EB=arb('1e-4'), ECC=arb(1), ECR=arb('1e-3')) for j in range(3)}
YW[3] = dict(WD=arb(10), THR=arb('0.01'), EA=arb('0.05'), EB=arb('3e-3'), ECC=arb(1), ECR=arb('1e-3'))   # theta_3 has room
# levels <= 2: small D-side weights, because theta_2 enters item 4 through C5 (the coefficient majorants there are <= 1e-5)
DE0 = arb('1e-3')
DG0 = arb('1e-4')     # level 0: Young weight on N_1 for the ramp (T' w0, gamma_1 w1) pair
YW_REPORT = dict({str(j): {k: float(v.mid()) for k, v in w.items()} for j, w in YW.items()}, DE0=float(DE0.mid()), DG0=float(DG0.mid()))
PAIRS = n3forms.PAIRS
assert PAIRS == n3qforms.PAIRS
PI = {pr: n for n, pr in enumerate(PAIRS)}
ZERO = arb(0)


def hull_union(vals):
    lo = min(float(v.lower()) for v in vals); hi = max(float(v.upper()) for v in vals)
    return iv(lo, hi)


def up(v):
    """upper bound of |v| as an exact arb"""
    return arb(U(absb(v)))


def hidx(i, k): return PI[(min(i, k), max(i, k))]


def coef_arrays(V, Z, Zp):
    """N on the box: scalar majorants (Atp, A, A1, A2, B0, B1, dA3e, dCe) and |.|-arrays of the derivative sups."""
    out = n3forms.evaluate(V)
    Za = [arb(z) for z in Z]; Zpa = [arb(z) for z in Zp]
    arr = {k: [up(v) for v in out[k]] for k in ('g3', 'Xg3', 'XXg3', 'h3', 'h2', 'Xh2', 'XXh2')}
    S = dict(Atp=arb(max(0.0, U(out['At'][0]))), A=up(out['At'][0]), A1=up(out['DAt'][0]), A2=up(out['DDAt'][0]))
    S['B0'] = sum((Za[j] * arr['g3'][j] for j in range(3)), ZERO)                                    # |Bt| (mean value)
    S['B1'] = sum((Za[j] * arr['Xg3'][j] for j in range(3)), ZERO) + sum((arr['g3'][i] * Zpa[i] for i in range(3)), ZERO)
    S['dA3e'] = arr['g3'][2] * Zpa[2]                                                                # |d_e'' A3| |e'''|
    S['dCe'] = sum((Za[j] * arr['h2'][hidx(j, 2)] for j in range(3)), ZERO)                          # |d_e'' Ct|
    return S, arr


def effective(S, arr, Parr, th, l2):
    """(1 - th) N + th e'''' + th k^2 e'' + th P:  scaled N majorants, top coefficient shifted, P added to the Hessians"""
    f = 1 + th
    E = {k: f * v for k, v in S.items()}
    E['Atp'] = f * S['Atp'] + th * (1 + S['A'])            # sup((1-th)At + th)^+
    E['A'] = f * S['A'] + th
    E['dCe'] = f * S['dCe'] + l2[0]                         # d_e'' of th k^2 e''
    A = {k: [f * v for v in vals] for k, vals in arr.items()}
    for k in ('h2', 'Xh2', 'XXh2'):
        A[k] = [A[k][n] + th * Parr[k + 'P'][n] for n in range(len(PAIRS))]
    return E, A


def betas(j, A, Z, Zp, Zpp, l2, l0):
    """a_flat (level j) <= beta_e |e| + sum_{m<j} beta_m |w_m|: one jet times boxed factors.  j = 0: gamma_e, gamma_0, gamma_1
    with |C| <= gamma_e |e| + gamma_0 |w_0| + gamma_1 |w_1| (C = the whole lower-order part at level 0)."""
    Za = [arb(z) for z in Z]; Zpa = [arb(z) for z in Zp]; Zppa = [arb(z) for z in Zpp]
    h = lambda key, i, k: A[key][hidx(i, k)]
    if j == 3:                     # [C'']_flat + [B'']_flat e'''  (slots e, w0, w1, w2)
        be = [ZERO] * 4
        for i in range(3): be[i] = be[i] + sum((Za[k] * h('XXh2', i, k) for k in range(3)), ZERO) / 2
        for i in range(3): be[i + 1] = be[i + 1] + 2 * sum((Za[k] * h('Xh2', i, k) for k in range(3)), ZERO)
        for i in range(3): be[i + 1] = be[i + 1] + sum((h('h2', i, k) * Zpa[k] for k in range(3)), ZERO)
        for i in range(2): be[i + 2] = be[i + 2] + sum((Za[k] * h('h2', i, k) for k in range(3)), ZERO)
        be[3] = be[3] + (sum((Za[k] * A['XXg3'][k] for k in range(3)), ZERO) + 2 * sum((A['Xg3'][i] * Zpa[i] for i in range(3)), ZERO)
                         + sum((h('h3', i, k) * Zpa[i] * Zpa[k] for i in range(3) for k in range(3)), ZERO)
                         + sum((A['g3'][i] * Zppa[i] for i in range(3)), ZERO))
        be[3] = be[3] + 2 * l2[1]; be[2] = be[2] + l2[2] + l0[0]; be[1] = be[1] + 2 * l0[1]; be[0] = be[0] + l0[2]
    elif j == 2:                   # [C']_flat = d_X C + d_e C e' + d_e' C e''  (slots e, w0, w1)
        be = [sum((Za[k] * h('Xh2', i, k) for k in range(3)), ZERO) / 2 for i in range(3)]
        be[1] = be[1] + sum((Za[k] * h('h2', 0, k) for k in range(3)), ZERO)
        be[2] = be[2] + sum((Za[k] * h('h2', 1, k) for k in range(3)), ZERO)
        be[2] = be[2] + l2[1]; be[1] = be[1] + l0[0]; be[0] = be[0] + l0[1]
    elif j == 1:                   # C(e, e', 0)  (slots e, w0)
        be = [sum((Za[k] * h('h2', i, k) for k in range(2)), ZERO) / 2 for i in range(2)]
        be[0] = be[0] + l0[0]
    else:                          # level 0: C(e, e', e'') + th k^2 e'' - th c0 e
        be = [sum((Za[k] * h('h2', i, k) for k in range(3)), ZERO) / 2 for i in range(3)]
        be[0] = be[0] + l0[0]; be[2] = be[2] + l2[0]
    return be


def psi(j, E, s0, s1, s2, s3):
    """alpha_0, alpha_1 and |Psi_{j,j+2}|, |Psi_{j,j+1}| over T (core: s0 = 1, s_m = |tau_m|) or over sqrt T (ramp: s0 = sqrt T,
    s_m = |T^(m)|/sqrt T); the Psi are those of LEMMA_F 6.1 (checked in n3wcheck.py)."""
    A, A1, A2, B0, B1, d3, dC = E['A'], E['A1'], E['A2'], E['B0'], E['B1'], E['dA3e'], E['dCe']
    if j == 3:
        a0 = A2 + 2 * B1 + d3 + dC; a1 = 2 * A1 + B0
        p2 = 2 * a0 * s0 + 4 * s2 * A + s1 * (6 * A1 + 4 * B0)
        p1 = s1 * (2 * A2 + 4 * B1 + 4 * d3 + 4 * dC) + 2 * s2 * B0 + 2 * s3 * A
    elif j == 2:
        a0 = B1 + dC; a1 = A1 + B0
        p2 = 2 * a0 * s0 + 4 * s2 * A + s1 * (2 * A1 + 4 * B0)
        p1 = s1 * (4 * dC + 2 * A2) + 2 * s2 * (A1 + B0) + 2 * s3 * A
    else:
        a0 = dC; a1 = B0
        p2 = 2 * a0 * s0 + 4 * s2 * A + s1 * (2 * A1 + 4 * B0)
        p1 = s1 * (4 * dC + 4 * B1 + 2 * A2) + s2 * (2 * B0 + 4 * A1) + 2 * s3 * A
    return a0, a1, p2, p1


def frenet(k, p0, q0, n0, t0):
    P = [p0]; Qq = [q0]; N = [n0]; T1 = [t0]
    for m in range(4):
        P.append(-sum((comb(m, l) * k[l] * Qq[m - l] for l in range(m + 1)), ZERO))
        Qq.append((arb(1) if m == 0 else ZERO) + sum((comb(m, l) * k[l] * P[m - l] for l in range(m + 1)), ZERO))
        N.append(-sum((comb(m, l) * k[l] * T1[m - l] for l in range(m + 1)), ZERO))
        T1.append(sum((comb(m, l) * k[l] * N[m - l] for l in range(m + 1)), ZERO))
    return P, Qq, N, T1


def q_strong(Vb, Z, th, xbar, cabs, babs):
    """majorants of g^(k) (k <= 3 with the box; k = 2, 3, 4 also with the top jet e^(k) removed) and of the top coefficient
    |d_e g| = |U k|/alpha^2, for g = U/alpha + 5 th k k'"""
    k = [Vb[f'k{m}'] for m in range(6)]
    P, Qq, N, T1 = frenet(k, Vb['p'], Vb['q'], Vb['n'], Vb['t1'])
    Ve = {f'k{m}': k[m] for m in range(5)}
    for m in range(4): Ve[f'e{m}'] = iv(-Z[m], Z[m])
    Ve['e4'] = arb(0)
    inv = n3qforms.invalpha(Ve)
    flat = {4: inv[4]}
    for m in (2, 3):
        Vf = dict(Ve); Vf[f'e{m}'] = arb(0); flat[m] = n3qforms.invalpha(Vf)[m]
    dc = arb(xbar) + th * cabs; dt = arb(xbar) + th * babs
    def Ek(F, kk, fl):
        return sum((comb(kk, l) * F[l] * (flat[kk] if (fl and l == 0) else inv[kk - l]) for l in range(kk + 1)), ZERO)
    G, Gf = {}, {}
    for kk in range(5):
        Ekap = sum((comb(kk, l) * k[l] * k[kk - l + 1] for l in range(kk + 1)), ZERO)
        if kk <= 3: G[kk] = up(dc * absb(Ek(Qq, kk, False)) + dt * absb(Ek(T1, kk, False)) + 5 * th * absb(Ekap))
        if kk >= 2: Gf[kk] = up(dc * absb(Ek(Qq, kk, True)) + dt * absb(Ek(T1, kk, True)) + 5 * th * absb(Ekap))
    al = 1 - k[0] * iv(-Z[0], Z[0])
    Gtop = up((dc * absb(Qq[0]) + dt * absb(T1[0])) * absb(k[0]) / (al * al))
    return G, Gf, Gtop


T7CONST = os.path.join(HERE, '..', '..', 't7', 'constants.json')


def xr_per_slab(slabs):
    """residual part of (M), x_r = mbar |c_phi|_2 R2 per slab (as camp.xr_per_slab; review B6)"""
    t7c = json.load(open(T7CONST)); R6 = json.load(open(os.path.join(cconst.T6DIR, 'residuals.json')))
    by6 = {r['slab']: r for r in R6['per_slab']}
    cphi = max((arb(repr(sum(v * v for v in c['c_phi']))).sqrt() for c in t7c['certificates']), key=lambda v: float(v.upper()))
    return {s: U(arb(repr(t7c['slabs'][str(s)]['C3']['mbar'])) * cphi * arb(repr(by6[s]['norms']['R2']))) for s in slabs}


def run_slab(i, amp, cons, Y, levels, fac=1.25, b2=1e-2, etaF=1e-2, eta1=1e-5, b0=7e-4, b1=1.65e-4, withQ=True):
    J, rec = cconst.load_slab(i)
    S = rec['scalars']
    ell = iv(*S['range_ell']); yl = ell - X0
    ps = amp['per_slab']; ka = amp['slabs'].index(i)
    sF0 = arb(repr(ps['0']['aF'][ka])) * arb(Y) + arb(repr(ps['0']['bF'][ka]))
    pc = [q for q in cons['per_slab'] if q['slab'] == i][0]
    xbar = float(U(arb(repr(cons['t7']['Cx'])) * arb(Y) + XR[i])) if withQ else 0.0     # (M): C_x Ybar + x_r (review B6)
    th = arb(U(arb(xbar) / arb(repr(pc['norms']['nV'])))) if withQ else ZERO
    cabs = arb(max(abs(v) for v in S['range_c'])); babs = arb(max(abs(v) for v in S['range_b']))
    LN = arb(cconst.YB / cconst.YA).log()
    x_a = arb(cconst.YA) + X0; x_b = arb(cconst.YB) + X0
    xa_s = L(x_a); xb_s = math.nextafter(U(x_b), math.inf); xmax = float(ell.upper())
    pcs = list(cconst.pieces(J, xa_s, xmax + 1.0, [xa_s, xb_s, 30.25 + 1e-9, 34.25], cconst.NSUB_TIP, (xa_s, xb_s)))
    # ---- x-ranges for the coefficient majorants: ramp, [x_b, 30.25], [30.25, 34.25], then dyadic in x up to ell
    edges = [xa_s, xb_s, 30.25 + 1e-9, 34.25]
    while edges[-1] < xmax: edges.append(edges[-1] * 2.0)
    groups = {}
    for e, c, xi, ln, fx, u0, u1 in pcs:
        g = max(q for q in range(len(edges) - 1) if u0 >= edges[q] * (1 - 1e-15)) if u0 >= edges[0] * (1 - 1e-15) else 0
        groups.setdefault(g, []).append((e, c, xi, ln, u0, u1))
    out = dict(slab=i, box=dict(fac=fac, b0=b0, b1=b1, b2=b2, etaF=etaF, eta1=eta1, Y=Y, sqrtF0=U(sF0)),
               mod=dict(xbar=xbar, thbar=U(th), nV=pc['norms']['nV']), ranges=[], lev={})
    zero_acc = lambda: dict(eH=ZERO, eH1=ZERO, e0=ZERO, e01=ZERO, fR=ZERO, f1R=ZERO,
                            Eabs2=ZERO, EF2=ZERO, E1abs2=ZERO, E1F2=ZERO, ER2=ZERO)
    acc = {j: zero_acc() for j in levels}
    sup = {j: {} for j in levels}
    def smax(j, key, v):
        u = U(v) if not isinstance(v, float) else v
        if u > sup[j].get(key, -math.inf): sup[j][key] = u
    diag = dict(Atp=0.0, Atp_x=None, ke=0.0, ke_x=None)
    eta1a = arb(eta1)
    for g in sorted(groups):
        cells = groups[g]
        xlo, xhi = min(cc[4] for cc in cells), max(cc[5] for cc in cells)
        Vb = {f'k{m}': hull_union([H(J[f'k{m}'][e, c]) for e, c, *_ in cells]) for m in range(6)}
        Vb['A'] = hull_union([H(J['A'][e, c]) for e, c, *_ in cells])
        for nm in ('r1', 'r2', 'r3', 'p', 'q', 'n', 't1'): Vb[nm] = hull_union([H(J[nm][e, c]) for e, c, *_ in cells])
        ylo, yhi = xlo - 2.2490776515, xhi - 2.2490776515                 # floats; E0 carries a 1e-12 cushion
        core_box = ylo >= 28.0 - 1e-9
        if core_box:
            E0 = fac * (eta1 + U(sF0) * yhi ** 1.5 / math.sqrt(84.0)) * (1 + 1e-12)      # |e(28)| <= eta1 (tip box)
            Z = (E0, b0, b1 / xlo, b2 / xlo ** 2)
        else:
            Z = (eta1, eta1 if xhi <= 34.25 else etaF, etaF / xlo, etaF / xlo ** 2)
        Zp = (Z[1], Z[2], Z[3]); Zpp = (Z[2], Z[3], 0.0)
        V = dict(Vb)
        for m in range(4): V[f'e{m}'] = iv(-Z[m], Z[m])
        Sn, arr = coef_arrays(V, Z[:3], Zp)
        Parr = {k: [up(v) for v in vals] for k, vals in n3qforms.hessP(V).items()}
        k = [Vb[f'k{m}'] for m in range(6)]
        l2 = [up(th * k[0] * k[0]), up(th * 2 * k[0] * k[1]), up(th * 2 * (k[1] * k[1] + k[0] * k[2]))]
        l0 = [up(th * (3 * k[1] * k[1] + 4 * k[0] * k[2])), up(th * (10 * k[1] * k[2] + 4 * k[0] * k[3])),
              up(th * (10 * k[2] * k[2] + 14 * k[1] * k[3] + 4 * k[0] * k[4]))]
        Ef, Af = effective(Sn, arr, Parr, th, l2)
        if withQ: G, Gf, Gtop = q_strong(Vb, Z, th, xbar, cabs, babs)
        else: G, Gf, Gtop = {m: ZERO for m in range(4)}, {m: ZERO for m in (2, 3, 4)}, ZERO
        ke = U(absb(Vb['k0'])) * Z[0]
        if ke > diag['ke']: diag['ke'], diag['ke_x'] = ke, (xlo, xhi)
        if U(Sn['Atp']) > diag['Atp']: diag['Atp'], diag['Atp_x'] = U(Sn['Atp']), (xlo, xhi)
        rng = dict(x=[xlo, xhi], core_box=core_box, Z=list(Z), ke=ke, **{k_: U(v) for k_, v in Sn.items()},
                   G=[U(G[m]) for m in range(4)], Gtop=U(Gtop))
        bet = {j: betas(j, Af, Z[:3], Zp, Zpp, l2, l0) for j in levels}
        rng['beta'] = {j: [U(v) for v in bet[j]] for j in levels}
        for e, c, xi, ln, u0, u1 in cells:
            y = xi - X0; lnA = arb(ln)
            ramp = u1 <= xb_s
            H2 = pos((y ** 3 - arb(28) ** 3) / 84)
            for j in levels:
                p = 2 * j - 2; be = bet[j]
                smax(j, 'nu', 2 * Ef['Atp'])
                if not ramp:
                    w = cconst.core_weight(y, yl, p); tau = w['tau']; T0 = w['T0']; Yv = w['Y']
                    t1, t2, t3 = absb(tau[1]), absb(tau[2]), absb(tau[3])
                    smax(j, 't2', t2)
                    if j >= 1:
                        a0, a1, m2, m1 = psi(j, Ef, arb(1), t1, t2, t3)
                        smax(j, 'm2', m2); smax(j, 'm1', m1); smax(j, 'm45', 2 * a1); smax(j, 'm33', 2 * t2 * a0)
                        for m in range(j):
                            rat = Yv ** (2 * (j - m))
                            smax(j, f'dC{m}', be[m + 1] ** 2 * rat); smax(j, f'd1C{m}', be[m + 1] ** 2 * t1 * t1 * rat)
                        acc[j]['eH'] += lnA * T0 * be[0] ** 2 * H2; acc[j]['eH1'] += lnA * T0 * t1 * t1 * be[0] ** 2 * H2
                        acc[j]['e0'] += lnA * T0 * be[0] ** 2 * eta1a ** 2; acc[j]['e01'] += lnA * T0 * t1 * t1 * be[0] ** 2 * eta1a ** 2
                    else:
                        A_, A1, A2_, B0 = Ef['A'], Ef['A1'], Ef['A2'], Ef['B0']
                        smax(j, 'm12', 2 * (A1 + B0)); smax(j, 'm02', 4 * t2 * A_ + 2 * t1 * (2 * A1 + B0))
                        smax(j, 'm01', 2 * t3 * A_ + 4 * t2 * A1 + 2 * t1 * A2_)
                        smax(j, 'G0', be[1]); smax(j, 'G1', be[2]); smax(j, 'G0t', t1 * be[1]); smax(j, 'G1t', t1 * be[2])
                        acc[j]['Eabs2'] += lnA * T0 * be[0] ** 2 * 2 * eta1a ** 2; acc[j]['EF2'] += lnA * T0 * be[0] ** 2 * 2 * H2
                        acc[j]['E1abs2'] += lnA * T0 * t1 * t1 * be[0] ** 2 * 2 * eta1a ** 2
                        acc[j]['E1F2'] += lnA * T0 * t1 * t1 * be[0] ** 2 * 2 * H2
                    # Q, strong form (transport g e' and the constant zeroth order; the latter is added in the assembly)
                    qs = (2 * j + 1) * G[1] + G[0] * t1 + (2 * Gtop * arb(Z[1]) if j >= 1 else ZERO)
                    smax(j, 'Qsig', qs)
                    for kk in range(2, j + 2):
                        GG = Gf[kk] if kk == j + 1 else G[kk]
                        smax(j, f'Qc{j + 1 - kk}', 2 * comb(j + 1, kk) * GG * Yv ** (kk - 1))
                else:
                    w = cconst.ramp_weight(y, p, LN); T = w['T']; Tt = w['Tt']; Zw = w['Z']
                    Winv = (1 + xi * xi) ** 2
                    rt = (28 * Zw).sqrt()
                    sqT, r1_, r2_, r3_ = absb(Tt[0]) / rt, absb(Tt[1]) / rt, absb(Tt[2]) / rt, absb(Tt[3]) / rt
                    smax(j, 'MfR', r2_ * r2_ * Winv); smax(j, 'T2pR', absb(T[2]) * Winv)
                    if j >= 1:
                        a0, a1, P2, P1 = psi(j, Ef, sqT, r1_, r2_, r3_)
                        smax(j, 'M2', P2 * P2 * Winv); smax(j, 'M1', P1 * P1 * Winv)
                        smax(j, 'M33', 2 * absb(T[2]) * a0 * Winv); smax(j, 'm45', 2 * a1)
                        for m in range(j):
                            smax(j, f'dR{m}', be[m + 1] ** 2 * y ** (2 * (j - m)))
                            smax(j, f'k1R{m}', be[m + 1] ** 2 * r1_ * r1_ * Winv)
                        acc[j]['fR'] += lnA * absb(T[0]) * be[0] ** 2 * eta1a ** 2; acc[j]['f1R'] += lnA * r1_ * r1_ * be[0] ** 2 * eta1a ** 2
                    else:
                        A_, A1, A2_, B0 = Ef['A'], Ef['A1'], Ef['A2'], Ef['B0']
                        smax(j, 'm12', 2 * (A1 + B0))
                        P02 = 4 * r2_ * A_ + 2 * r1_ * (2 * A1 + B0); P01 = 2 * r3_ * A_ + 4 * r2_ * A1 + 2 * r1_ * A2_
                        smax(j, 'M02', P02 * P02 * Winv); smax(j, 'M01', P01 * P01 * Winv)
                        smax(j, 'G0', be[1]); smax(j, 'G1', be[2])
                        smax(j, 'Mg0', absb(T[1]) * be[1] * Winv); smax(j, 'Mg1', r1_ * r1_ * be[2] ** 2 * Winv)
                        smax(j, 'MT1', r1_ * r1_ * Winv)
                        acc[j]['Eabs2'] += lnA * absb(T[0]) * be[0] ** 2 * eta1a ** 2
                        acc[j]['ER2'] += lnA * absb(T[0]) * be[0] ** 2 * eta1a ** 2
                    qk = ((2 * j + 1) * G[1] + (2 * Gtop * arb(Z[1]) if j >= 1 else ZERO)) * absb(T[0]) + G[0] * absb(T[1])
                    smax(j, 'Qk', qk * Winv)
                    for kk in range(2, j + 2):
                        GG = Gf[kk] if kk == j + 1 else G[kk]
                        smax(j, f'Qc{j + 1 - kk}', 2 * comb(j + 1, kk) * GG * y ** (kk - 1))
        out['ranges'].append(rng)
    for j in levels:
        out['lev'][str(j)] = assemble(j, sup[j], acc[j], xbar, th, cabs) if j >= 1 else assemble0(sup[0], acc[0], xbar, th, cabs)
    l3 = out['lev'].get('3')
    if l3:                                       # level-3 keys of the first version (camp reads 'lev' when present)
        out.update(dtheta3=l3['dtheta'], dsigma3=l3['dsigma'], dk33=l3['dk'], d3=l3['d'], k3=l3['k'], hf3=l3['hf'])
    out.update(sup_Atp=diag['Atp'], sup_Atp_x=diag['Atp_x'], sup_ke=diag['ke'], sup_ke_x=diag['ke_x'])
    return out


def assemble(j, sp_, ac, xbar, th, cabs):
    """level j >= 1 (LEMMA_F 6.1): weak form for (1-th)N + th e'''' + th k^2 e'' - th c0 e + th P, strong form for Q_s"""
    g = lambda k_: arb(sp_.get(k_, 0.0))
    yw = YW[j]; WD, THR, EA, EB, ECC, ECR = yw['WD'], yw['THR'], yw['EA'], yw['EB'], yw['ECC'], yw['ECR']
    nu, m2, m1, m45, m33 = g('nu'), g('m2'), g('m1'), g('m45'), g('m33')
    M2, M1, M33, MfR, T2pR, t2s = g('M2'), g('M1'), g('M33'), g('MfR'), g('T2pR'), g('t2')
    eps1 = sqrtp(4 * M1 / T2pR) if U(T2pR) > 0 and U(M1) > 0 else arb('1e-3')
    cD = nu + (m2 + m45) * WD / 2 + THR + EA
    cN = m1 / 2 + m45 / (2 * WD) + eps1 / 2 + EB
    cF = m2 / (2 * WD) + m1 / 2 + m33 + ECC * t2s * t2s
    cV = M2 / (4 * THR) + M1 / (2 * eps1) + M33 + ECR * MfR
    Qsig = g('Qsig') + 2 * (arb(xbar) + th * cabs)
    dtheta = cD + cN
    dsigma = cF + cN * (arb('0.25') + t2s / 2) + Qsig
    dk = cV + cN * T2pR / 2 + g('Qk')
    nt = j + 1                                    # number of terms in the a_flat bound
    fC = 2 / EA + 1 / ECC; fR = 2 / EA + 1 / ECR; f1 = 8 / EB
    d = {m: U(nt * (fC * g(f'dC{m}') + f1 * g(f'd1C{m}') + fR * g(f'dR{m}'))
              + (2 * nt * (fC * ac['eH'] + f1 * ac['eH1']) if m == 0 else 0)) for m in range(j)}
    k = {m: U(nt * f1 * g(f'k1R{m}')) for m in range(j)}
    c = {m: U(g(f'Qc{m}')) for m in range(j)}
    hf = U(2 * nt * (fC * ac['e0'] + f1 * ac['e01']) + nt * (fR * ac['fR'] + f1 * ac['f1R']))
    return dict(dtheta=U(dtheta), dsigma=U(dsigma), dk=U(dk), d=d, k=k, c=c, hf=hf,
                parts=dict(sp_, cD=U(cD), cN=U(cN), cF=U(cF), cV=U(cV), Qsig_total=U(Qsig)))


def assemble0(sp_, ac, xbar, th, cabs):
    """level 0: S0 = 2 int T At w2^2 + Psi12 w1 w2 + Psi02 w0 w2 + Psi01 w0 w1 - 2 int (T w1 + T' w0) C (LEMMA_F 6.1)"""
    g = lambda k_: arb(sp_.get(k_, 0.0))
    yw = YW[0]; WD, THR = yw['WD'], yw['THR']
    nu, m12, m02, m01, t2s = g('nu'), g('m12'), g('m02'), g('m01'), g('t2')
    M02, M01, T2pR, MT1 = g('M02'), g('M01'), g('T2pR'), g('MT1')
    G0, G1, G0t, G1t, Mg0, Mg1 = g('G0'), g('G1'), g('G0t'), g('G1t'), g('Mg0'), g('Mg1')
    Eabs, EF, E1abs, E1F = sqrtp(ac['Eabs2']), sqrtp(ac['EF2']), sqrtp(ac['E1abs2']), sqrtp(ac['E1F2'])
    eps01 = sqrtp(4 * M01 / T2pR) if U(T2pR) > 0 and U(M01) > 0 else arb('1e-3')
    sMg1 = sqrtp(Mg1)
    cD = nu + (m12 + m02) * WD / 2 + THR
    cN = m12 / (2 * WD) + m01 / 2 + eps01 / 2 + 2 * G1 + G0 + EF + DE0 + G1t + DG0
    cF = m02 / (2 * WD) + m01 / 2 + G0 + EF + 2 * G0t + G1t + 2 * E1F + DE0
    cV = M02 / (4 * THR) + M01 / (2 * eps01) + 2 * Mg0 + Mg1 / DG0 + DE0
    hf = (ac['Eabs2'] + ac['E1abs2'] + MT1 * ac['ER2']) / DE0
    Qsig = g('Qsig') + 2 * (arb(xbar) + th * cabs)
    dtheta = cD + cN
    dsigma = cF + cN * (arb('0.25') + t2s / 2) + Qsig
    dk = cV + cN * T2pR / 2 + g('Qk')
    return dict(dtheta=U(dtheta), dsigma=U(dsigma), dk=U(dk), d={}, k={}, c={}, hf=U(hf),
                parts=dict(sp_, cD=U(cD), cN=U(cN), cF=U(cF), cV=U(cV), Qsig_total=U(Qsig), Eabs=U(Eabs), EF=U(EF)))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--slabs', required=True)
    ap.add_argument('--levels', default='0-3')
    ap.add_argument('--amp', default=os.path.join(HERE, 'constants', 'amp_fin.json'))
    ap.add_argument('--consts', default=os.path.join(HERE, 'constants', 'constants_fin.json'))
    ap.add_argument('--Y', default='1.7e-12')
    ap.add_argument('--fac', type=float, default=1.7)
    ap.add_argument('--b2', type=float, default=1e-2)
    ap.add_argument('--etaF', type=float, default=1e-2)
    ap.add_argument('--eta1', type=float, default=1e-5, help="tip box |e|, |e'| on x <= 34.25 (T5 2.5 box; review A3)")
    ap.add_argument('--b0', type=float, default=1.5e-3, help="core box |e'| for Ybar <= Y")
    ap.add_argument('--b1', type=float, default=3e-4, help="core box x|e''| for Ybar <= Y")
    ap.add_argument('--noQ', action='store_true', help='leave out the modulation terms Q (for comparison only)')
    ap.add_argument('--out', default=os.path.join(HERE, 'constants', 'n3.json'))
    a = ap.parse_args()
    amp = json.load(open(a.amp)); cons = json.load(open(a.consts))
    global XR
    XR = xr_per_slab(cconst.parse_slabs(a.slabs))
    levels = cconst.parse_slabs(a.levels)
    res = {}; t0 = time.time()
    for s in cconst.parse_slabs(a.slabs):
        o = run_slab(s, amp, cons, float(a.Y), levels, a.fac, a.b2, a.etaF, a.eta1, a.b0, a.b1, not a.noQ); res[s] = o
        print(f"slab {s}: " + ' | '.join(f"j{j} dth {o['lev'][str(j)]['dtheta']:.2e} dsig {o['lev'][str(j)]['dsigma']:.2e} "
                                         f"dk {o['lev'][str(j)]['dk']:.2e}" for j in levels)
              + f" | sup At+ {o['sup_Atp']:.3e}, sup|k e| {o['sup_ke']:.2e} ({time.time()-t0:.0f}s)", flush=True)
    summ = {}
    for j in levels:
        lv = [o['lev'][str(j)] for o in res.values()]
        summ[j] = dict(dtheta=max(v['dtheta'] for v in lv), dsigma=max(v['dsigma'] for v in lv), dk=max(v['dk'] for v in lv),
                       hf=max(v['hf'] for v in lv), d={m: max(v['d'][m] for v in lv) for m in range(j)},
                       k={m: max(v['k'][m] for v in lv) for m in range(j)}, c={m: max(v['c'][m] for v in lv) for m in range(j)},
                       Qsig=max(v['parts']['Qsig_total'] for v in lv))
    summ['sup_Atp'] = max(o['sup_Atp'] for o in res.values()); summ['sup_ke'] = max(o['sup_ke'] for o in res.values())
    summ['thbar'] = max(o['mod']['thbar'] for o in res.values()); summ['xbar'] = max(o['mod']['xbar'] for o in res.values())
    args = dict(vars(a), argv=sys.argv[1:], t6dir=cconst.T6DIR, t7const=T7CONST, young=YW_REPORT,
                env={k: os.environ[k] for k in ('SUBW', 'T9JSON', 'T9WIN', 'T6DIR', 'TOPLOSS', 'NLFAC', 'BASE_T1', 'YBAR', 'T6RES', 'SLAB0', 'YOUT', 'OMP_NUM_THREADS') if k in os.environ})
    json.dump(dict(args=args, summary=summ, per_slab=res), open(a.out, 'w'), indent=1)
    print('summary', json.dumps(summ))


if __name__ == '__main__':
    main()
