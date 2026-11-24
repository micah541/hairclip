"""(SE) settlement (task A/B): the corrected base
    kap^ = kap~ + dk,   dk = (a1/ell) b1(u) + (a2/ell^2) b2(u),   u = x/ell(tau) in [0.6, 0.8],
    b1 = B, b2 = B', B(u) = psi((u - 0.6)/0.2)/(0.2 I_psi), psi = t^6 (1-t)^6 (C^5), int B du = 1, int b2 = 0, int u b2 = -1,
with a1 = -delta(tau) (turning defect) and a2 the root of  a2 J2 = -(h + a1 ell J1 + N)  (far-point height h, quadratic N),
is exactly symmetric about the axis at the far end (turning exactly pi, far point on the axis), and unchanged on x < 0.6 ell:
so (SE) holds exactly, the tip and all tip-region quantities (R2, Rs, Rs_chi, T7 certificates) are unchanged.
Certified per slab (sup over the slab): delta from T6 (sup_turn_def); h from T1's 192-bit stage offsets (se_v2) + the rate
|h'| <= |V(ell)| delta + c(1.01 ell delta + h) + b delta + |A(ell)| + |eps_ell| between stages; |delta'| <= k(ell) eps_ell + |r'(ell)|
+ k(ell)|A(ell)| (identity d/dtau theta(ell) = k eps_ell + rho_theta(ell)); then a1, a2, a1', a2', the residual change dr (r, r', r'',
r''', S1 and A) and the increments of every far-field T6 norm (R1, RIs, RIs3, RF/RFa/RFb_j opt4n, far-end scalars).
Writes out/se_corr.json and out/t6se/residuals.json (T6 v2 norms + increments; scalars of the corrected base)."""
import sys, json, os, math, copy
import numpy as np
from flint import arb, arb_poly, ctx
ctx.prec = 128
S4 = '/home/claude/s1spec/step4'
sys.path.insert(0, f'{S4}/t6'); import reader as RD
UA, UB = 0.6, 0.8; WU = UB - UA
DS1 = arb('2.500806276619802'); X0 = arb('2.2490776515')
HGAP = arb('0.0973') * arb('0.15') * arb('1.001')       # half the largest Radau-node gap (0.4528 -> 0.6474) in tau


def ub(x): return math.nextafter(float(arb(x).abs_upper()), math.inf)   # outward


def bump_norms(NS=400):
    """sup and L1 (on [0,1]) of psi^(k), psi = t^6 (1-t)^6, k <= 9: per cell the Taylor-shifted polynomial at the cell
    midpoint (exact coefficients), |p(t_m + s)| <= sum |c_j| r^j (rigorous upper bounds, no ball-evaluation blow-up)"""
    p = arb_poly([0, 0, 0, 0, 0, 0, 1]) * arb_poly([1, -1]) ** 6
    sup, l1 = [], []
    r = arb(1) / (2 * NS)
    for k in range(10):
        m = arb(0); s = arb(0)
        for i in range(NS):
            tm = arb(2 * i + 1) / (2 * NS)
            q = p(arb_poly([tm, 1]))                                   # Taylor shift to the midpoint
            v = sum((abs(c) * r ** j for j, c in enumerate(q.coeffs())), arb(0))
            m = max(m, v, key=ub); s += v / NS
        sup.append(arb(ub(m))); l1.append(arb(ub(s))); p = p.derivative()
    return sup, l1


def main():
    T6 = json.load(open(f'{S4}/t6/residuals.json')); SE = json.load(open(f'{S4}/t1/se_polished_v2.json'))
    psup, pl1 = bump_norms()
    Ip = arb(math.factorial(6) ** 2) / math.factorial(13); W_ = arb(WU)
    # in u: B^(k) sup = psup[k]/(W^(k+1) Ip),  L1 = pl1[k]/(W^k Ip)
    Bs = [psup[k] / (W_ ** (k + 1) * Ip) for k in range(10)]; Bl = [pl1[k] / (W_ ** k * Ip) for k in range(10)]
    rows = []; newT6 = copy.deepcopy(T6)
    for r_ in newT6['per_slab']:            # T6 stores ranges round-to-nearest: widen outward by 2 ulps (sup_ scalars: merge.py)
        for k_, v_ in r_['scalars'].items():
            if k_.startswith('range_') and isinstance(v_, list) and len(v_) == 2:
                r_['scalars'][k_] = [math.nextafter(math.nextafter(v_[0], -math.inf), -math.inf), math.nextafter(math.nextafter(v_[1], math.inf), math.inf)]
    for rec, nrec in zip(T6['per_slab'], newT6['per_slab']):
        i = rec['slab']; sc = rec['scalars']; J = RD.jets(i); se = SE[i]; assert se['slab'] == i
        x = J['x']; ell_lo, ell_hi = sc['range_ell']; ellt = arb(max(abs(v) for v in sc['range_ell_tau']))
        L, Lh = arb(ell_lo), arb(ell_hi)
        cb = arb(max(abs(v) for v in sc['range_c'])); bb = arb(max(abs(v) for v in sc['range_b']))
        H = lambda nm, sel: arb(float(np.max(np.abs(J[nm][sel])))) if np.any(sel) else arb(0)
        xlo, xhi = x[..., 0], x[..., 1]
        win = (xhi >= UA * ell_lo) & (xlo <= UB * ell_hi)
        far = xhi >= UA * ell_lo
        endc = np.zeros(x.shape[:2], bool); endc[-1, :] = True
        K = [H(f'k{j}', win) for j in range(6)]; Kf = [H(f'k{j}', far) for j in range(6)]; Qw = H('Q', win)
        thw = J['theta'][win]; cs = arb(float(np.min(thw[..., 0]))).union(arb(float(np.max(thw[..., 1])))).cos()
        assert float(cs.upper()) < 0, f'slab {i}: cos theta not negative on the window'
        cmin = -arb(float(cs.upper()))
        # ---- far-end defects
        dlt = arb(sc['sup_turn_def']); Aend = arb(sc['sup_A_end']); epsl = arb(sc['sup_eps_ell'])
        kend = H('k0', endc); r1end = H('r1', endc); Vend = H('V', endc)
        hrow = arb(se['offset_max']) + 10 * Lh ** 2 * (rec['n'] + 4) * arb(se['Yn_rad_max']) + arb('1e-30')
        H0 = Vend * dlt + cb * arb('1.01') * Lh * dlt + bb * dlt + Aend + epsl
        h = (hrow + HGAP * H0) / (1 - HGAP * cb); hp = H0 + cb * h
        rend = json.load(open(f'{S4}/p2/out/rhoend/slab_{i:03d}.json'))['sup_rho_theta_end']
        dlt_p = kend * epsl + arb(rend)                                  # d/dtau theta(ell) = k(ell) eps_ell + rho_theta(ell)  (exact)
        # theta_tau bounds (T1 local form + rho_theta = r' + k A)
        def thtau(sel):
            g = lambda nm: np.max(np.abs(J[nm]), axis=-1)[sel]
            v = g('k3') + g('k0') ** 2 * g('k1') + g('k0') * g('Q') + float(cb.upper()) * g('x') * g('k0') + g('r1') + g('k0') * g('A')
            return arb(float(np.max(v)) * (1 + 1e-12))
        Tth, Tthw = thtau(np.ones(x.shape[:2], bool)), thtau(win)
        # ---- a1, a2 and their tau-derivatives
        a1, a1p = dlt, dlt_p
        J1b = arb('0.4')                                                  # |J1| = |int_0.6^1 B1 cos| <= 0.4
        a2 = (h + a1 * Lh * J1b) / cmin
        for _ in range(3):
            Dth = a1 + a2 * Bs[0] / L; a2 = (h + a1 * Lh * J1b + Lh * Dth ** 2) / cmin
        Dth = a1 + a2 * Bs[0] / L
        J2p = Tthw + K[0] * ellt; LJ1p = J1b * ellt + Lh * (Tth + Kf[0] * ellt) * J1b
        a2p = arb(0)
        for _ in range(3):
            Dthp = a1p + a2p * Bs[0] / L + a2 * Bs[0] * ellt / L ** 2 + ellt / L * UB * (a1 * Bs[0] + a2 * Bs[1] / L)
            Np = 2 * Lh * Dth * Dthp + ellt * Dth ** 2
            a2p = (hp + a1p * Lh * J1b + a1 * LJ1p + Np + a2 * J2p) / cmin
        # ---- dk in x: sup and L1 over the window
        d = [a1 * Bs[k] / L ** (k + 1) + a2 * Bs[k + 1] / L ** (k + 2) for k in range(8)]
        dl = [a1 * Bl[k] / L ** k + a2 * Bl[k + 1] / L ** (k + 1) for k in range(8)]      # int_window |dk^(k)| dx
        u = arb(UB)
        dtk = [a1p * Bs[k] / L ** (k + 1) + a2p * Bs[k + 1] / L ** (k + 2)
               + ellt / L * (a1 * ((k + 1) * Bs[k] + u * Bs[k + 1]) / L ** (k + 1) + a2 * ((k + 2) * Bs[k + 1] + u * Bs[k + 2]) / L ** (k + 2))
               for k in range(3)]                                       # |d_tau dk^(k)| at fixed x
        xb = Lh * UB; xa = L * UA; w = (UB * Lh - UA * L)
        DQ = 2 * K[1] * dl[1] + d[1] * dl[1]                            # Q^ - Q~ (constant beyond the window)
        dsum = lambda m: sum((d[a] for a in range(m + 1)), arb(0))
        Ks = K[0] + K[1] + K[2] + K[3] + K[4] + K[5]
        quad = lambda m: 20 * dsum(m + 4) ** 2 * (1 + Ks) ** 3
        # Delta Phi = -dk''' - (k^2 dk' + 2 k k' dk) + dk Q~ + k DQ + dk DQ - c x dk + O(dk^2); sup of d^j/dx^j (j = 0, 1, 2)
        def dPhi(j):                                                    # d^j/dx^j of (Delta Phi + c x dk): the dilation part is removed
            return (d[3 + j] + arb(4) ** j * 3 * Ks ** 2 * dsum(j + 1) + arb(2) ** j * (Qw + 2 * Ks ** 2) * dsum(j)
                    + arb(2) ** j * Ks * (DQ + 2 * Ks * dsum(j + 1) + dsum(j + 1) ** 2)
                    + arb(2) ** j * dsum(j) * (DQ + 2 * Ks * dsum(j + 1)) + quad(j))
        # comoving combination: d_tau dk + c d_x(x dk) = a1'/l B + a2'/l^2 B' + (c - l_tau/l)[...] - (l_tau/l)(a2/l^2) B'
        Qend = arb(max(abs(v) for v in sc['range_Q_end'])); cdef = (Qend + epsl) / L
        def com(m):
            return (a1p * Bs[m] / L ** (m + 1) + a2p * Bs[m + 1] / L ** (m + 2)
                    + cdef * (a1 * ((m + 1) * Bs[m] + UB * Bs[m + 1]) / L ** (m + 1) + a2 * ((m + 1) * Bs[m + 1] + UB * Bs[m + 2]) / L ** (m + 2))
                    + ellt / L * a2 * Bs[m + 1] / L ** (m + 2))
        Phw = dPhi(0)
        IPhi = dl[3] + K[0] ** 2 * dl[1] + 2 * K[0] * K[1] * dl[0] + Qw * dl[0] + w * (K[0] + d[0]) * DQ + w * quad(0)
        Tth_w = a1p + (a2p / L + a2 * ellt / L ** 2) * Bs[0] + cdef * UB * (a1 * Bs[0] + a2 * Bs[1] / L)   # |d_tau Dtheta + c x dk|
        ITth = a1p * w + (a2p + a2 * ellt / L) * (Lh / L) * Bl[0] + cdef * UB * Lh * (a1 * Bl[0] + a2 * Bl[1] / L)
        rth_w = Tth_w + Phw                                             # |d rho_theta| on the window
        rth_a = a1p + Kf[0] * DQ                                        # |d rho_theta| beyond the window (exact identity)
        IW = ITth + IPhi                                                # int_window |d rho_theta|
        rk_w = [com(m) + dPhi(m + 1) for m in range(2)]                 # |d rho_kappa^(m)| on the window (comoving form)
        rk_a = [Kf[1 + m] * DQ for m in range(2)]                       # beyond: d rho_kappa = -k' DQ
        rho0 = H('r1', far) + H('k0', far) * H('A', far); rA0 = H('r', far) + H('A', far)
        lena = Lh - xb; lenf = Lh - xa
        dz = IW + lena * rth_a + Dth * (rA0 + lenf * rho0)              # |dr|, |dA| on [x_a, ell]
        dz_w = IW + Dth * (rA0 + w * rho0)                              # on the window
        def derivs(dzr, rth, rk, Kr):
            r1 = rth + (Kr[0] + d[0]) * dzr + d[0] * H('A', far)
            r2 = rk[0] + Kr[1] * dzr + Kr[0] ** 2 * dzr + d[1] * (H('A', far) + dzr) + 3 * (Kr[0] + d[0]) * d[0] * (H('r', far) + dzr)
            r3 = (rk[1] + Kr[2] * dzr + 3 * Kr[0] * Kr[1] * dzr + Kr[0] ** 2 * r1
                  + 8 * (Kr[0] + Kr[1] + d[0] + d[1] + d[2]) ** 2 * (H('A', far) + H('r', far) + H('r1', far) + dzr + r1))
            return r1, r2, r3
        dr_w = derivs(dz_w, rth_w, rk_w, K); dr_a = derivs(dz, rth_a, rk_a, Kf)
        S_w = w * dz_w; S_a = w * dz_w + lena * dz                      # |dS1| (and |int_x^ell dr|) on window / after
        # ---- norm increments: window [x_a, x_b] and after [x_b, ell] with monotone weight majorants
        ya, yw, yb = xa - X0, xb - X0, Lh - X0
        def Wn(lo): return (1 + lo * lo) ** -2
        inc = {}
        for nm, (vw, va) in dict(R1=(S_w, S_a), RIs=(dr_w[0], dr_a[0]), RIs3=(dr_w[2], dr_a[2])).items():
            inc[nm] = (Wn(xa) * w * vw ** 2 + Wn(xb) * lena * va ** 2).sqrt()
        der_w = [S_a, dz_w, dr_w[0], dr_w[1], dr_w[2]]; der_a = [S_a, dz, dr_a[0], dr_a[1], dr_a[2]]
        for j in range(5):
            def wts(y):
                if j == 0: return 112 / y ** 2, 1792 / y ** 4, 28 * arb('17.2') * 16 / y ** 4
                if j == 1: return arb(28), arb(0), arb(0)
                return (28 * y ** (2 * j - 2), 28 * (2 * j - 2) ** 2 * y ** (2 * j - 4),
                        28 * ((2 * j - 2) * abs(2 * j - 3) + arb('5.6') * (2 * j - 2)) * y ** (2 * j - 4))
            ww = wts(yw if j >= 1 else ya); wa = wts(yb if j >= 1 else yw)
            for t, nm in enumerate(('RF', 'RFa', 'RFb')):
                inc[f'{nm}{j}_opt4n'] = (ww[t] * w * der_w[j] ** 2 + wa[t] * lena * der_a[j] ** 2).sqrt()
        sc_new = dict(sup_A_end=Aend + dz, sup_eps_ell=epsl + DQ, sup_closure_q=arb(0), sup_turn_def=arb(0), sup_turning_t1=arb(0),
                      sup_D_end=Aend + dz + epsl + DQ, sup_r_end=arb(sc['sup_r_end']) + dz, sup_S1_end=arb(sc['sup_S1_end']) + S_a)
        for k, v in inc.items():
            if k in nrec['norms']: nrec['norms'][k] = ub(arb(repr(nrec['norms'][k])) + v)
        for k, v in sc_new.items(): nrec['scalars'][k] = ub(v)
        # ---- widened jets of the corrected base (far region x >= 0.6 ell_lo), for T7/T8 constants and hulls
        Jn = {k: (np.array(v, copy=True) if isinstance(v, np.ndarray) and v.ndim == 3 else v) for k, v in J.items()}
        kk = K[0] + K[1]; kd = K[0] * d[1] + K[1] * d[0] + d[0] * d[1]
        DG = Lh * Dth; dp = DG + (Lh + 3) * Dth; dn = Dth
        dr1m, dr2m, dr3m = (max(a_, b_, key=ub) for a_, b_ in zip(dr_w, dr_a)); dzm = max(dz, dz_w, key=ub)
        dvar = 6 * kd + DQ + cb * dp + bb * dn
        wid = dict(Q=DQ, theta=Dth, p=dp, q=dp, n=dn, t1=dn, r=dzm, A=dzm, r1=dr1m, r2=dr2m, r3=dr3m,
                   A1=(Kf[0] + d[0]) * dzm + d[0] * H('r', far), Vstar=d[2] + cb * dp + bb * dn, V=d[2] + cb * dp + bb * dn + dzm,
                   varpistar=dvar, varpi=dvar + dzm, b1=dvar + 5 * kd + dzm, kr=(Kf[0] + d[0]) * dzm + d[0] * H('r', far))
        for j_ in range(6): wid[f'k{j_}'] = d[j_]
        for nm, dv in wid.items():
            if nm in Jn and isinstance(Jn[nm], np.ndarray) and Jn[nm].ndim == 3:
                v = ub(dv) * (1 + 1e-12)
                Jn[nm][far, 0] = np.nextafter(Jn[nm][far, 0] - v, -np.inf)          # outward (sub-ulp widths kept)
                Jn[nm][far, 1] = np.nextafter(Jn[nm][far, 1] + v, np.inf)
        os.makedirs(f'{S4}/p2/t6/jets', exist_ok=True)
        np.savez(f'{S4}/p2/t6/jets/jets_{i:03d}.npz', **{k: v for k, v in Jn.items() if isinstance(v, np.ndarray)})
        rows.append(dict(slab=i, k=rec['k'], a1=ub(a1), a2=ub(a2), a1p=ub(a1p), a2p=ub(a2p), h=ub(h), hp=ub(hp), cmin=float(cmin.lower()),
                         dk0=ub(d[0]), dz=ub(dz), dr=[ub(v) for v in dr_a], inc={k: ub(v) for k, v in inc.items()},
                         rel={k: ub(v) / max(rec['norms'].get(k, 0), 1e-300) for k, v in inc.items() if k in rec['norms']}))
    summ = {}                                # summary regenerated from the corrected per-slab data (as t6/merge.py)
    for nm in sorted(set(k for r in newT6['per_slab'] for k in r['norms'])):
        mx, at = max((r['norms'].get(nm, 0.0), r['slab']) for r in newT6['per_slab'])
        tg = newT6['targets'].get(nm)
        summ[nm] = dict(max=mx, at_slab=at) if tg is None else dict(max=mx, at_slab=at, target=tg, margin=(tg / mx if mx > 0 else None),
                                                                     fails=[r['slab'] for r in newT6['per_slab'] if r['norms'].get(nm, 0.0) > tg])
    for nm in sorted(set(k for r in newT6['per_slab'] for k in r['scalars'] if k.startswith('sup_'))):
        mx, at = max((r['scalars'][nm], r['slab']) for r in newT6['per_slab']); summ['scalar:' + nm] = dict(max=mx, at_slab=at)
    newT6['summary'] = summ
    os.makedirs(f'{S4}/p2/out/t6se', exist_ok=True); json.dump(newT6, open(f'{S4}/p2/t6/residuals.json', 'w'))
    newT6['se_correction'] = 'T6 v2 norms + certified increments of the (SE)-corrected base (p2/se_corr.py); far-end scalars of the corrected base'
    json.dump(newT6, open(f'{S4}/p2/out/t6se/residuals.json', 'w'))
    json.dump(rows, open(f'{S4}/p2/out/se_corr.json', 'w'), indent=1)
    worst = {}
    for r in rows:
        for k, v in r['rel'].items():
            if v > worst.get(k, (0, 0))[0]: worst[k] = (v, r['slab'])
    print('max relative increments:', {k: f'{v:.2e}@{s}' for k, (v, s) in worst.items()})
    print('max a1 %.2e a2 %.2e h %.2e hp %.2e |dr| %.2e; min cos %.4f' % (max(r['a1'] for r in rows), max(r['a2'] for r in rows),
          max(r['h'] for r in rows), max(r['hp'] for r in rows), max(r['dz'] for r in rows), min(r['cmin'] for r in rows)))


if __name__ == '__main__':
    main()
