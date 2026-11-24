"""Assemble T6 outputs:  python merge.py --in OUTDIR [--first 5 --last 187] [--json residuals.json] [--jets jets]
-> residuals.json (per-slab certified enclosures + summary vs the T5 targets) and the jets directory (copied npz + index)."""
import os, sys, json, glob, shutil, argparse

HERE = os.path.dirname(os.path.abspath(__file__))
TARGETS = {'R2': 2e-14, 'R1': 5e-13, 'Rs': 5e-11, 'RIs': 1e-10}
for j in range(5):
    t = 1e-7 if j <= 3 else 1e-5
    for nm in (f'RF{j}_T5', f'RF{j}_opt4n', f'RFlit{j}'): TARGETS[nm] = t

DEFS = {
    'base': 'slab i: unknowns = degree-8 polynomials in tc in [-1,1] (tau = tau0 + dtau (1+tc)/2) through the stored stage balls '
            'at nodes 0 and the radau.py Radau IIA nodes (60-digit binary values as exact rationals); per element kappa = sum d_q N_q(s), '
            'T1 Hermite-C4 basis and knots (outer knots XT (1 + C ell^3)^(1/3)); see base.py',
    'r': 'r = V~ + k~_xx - (1/4 - a~) p~ - b~ n~ (T5 1.1); computed as r = Re z, A = Im z, z = e^{i theta} int_0^x e^{-i theta} rho_theta, '
         'rho_theta = theta_tau - Phi; A is the residual part of the transport varpi',
    'R2': '||S2||_{L2_rho(0,34.25)}, S2 = int_0^x int_0^x r', 'Rs': '||r||_{L2_rho(0,34.25)}',
    'Rs_chi': '||chi_32 r||_rho (T7 clock-direction correction)',
    'R1': '||zeta_inf S1||_{W_inf}, whole half curve', 'RIs': '||zeta_inf r_x||_{W_inf}', 'RIs3': '||zeta_inf r_xxx||_{W_inf} (T7 gamma_2 correction)',
    'RF{j}_{design}': 'Lemma F R_j = (int T_j (d^{j-1} r)^2)^{1/2}, T_j = 28 eta_j^2 Yl^{2j-2} (flattened), d^{-1} r = S1 - S1(ell); '
                      'designs T5 (ramps 14,28) and opt4n (14,22 for j<=3; 22,30 for j=4)',
    'RFa{j}_{design}': '(int T_j\'^2/T_j (d^{j-1} r)^2)^{1/2}  (residual part of Lemma F R^_j, first term)',
    'RFb{j}_{design}': '(int |T_j\'\'| (d^{j-1} r)^2)^{1/2}  (residual part of Lemma F R^_j, second term)',
    'RFlit{j}': 'T5 (R) literal / T1 resid.py: (28 int eta^2 x^{2j-2} (d^{j-1} r)^2)^{1/2}, eta = S(ln(y/14)/ln 2), y = x + STIP, d^{-1} r = S1',
    'scalars': 'sup over the slab of |defects|: eps_ell = ell_tau + Q(ell) - c ell, eps_lnlam, eps_c, eps_t (scalar ODEs); closure_q = q~(ell) '
               '(far-end height above the axis), turning_t1 = t~1(ell) = cos theta(ell), A_end = A(ell), D_end = c q(ell) + b t1(ell) + A(ell) + eps_ell '
               '(= b1(ell) + ell_tau, the far-end transport defect), S1_end = S1(ell)',
    't7': 'crude certified sups (coefficient sums over element x slab): r, A, r_x, A_x = kappa r on [0,34.25] and on x >= 10.25; '
          'supArho_tip = sup |A rho\'/rho|, supAx_far = sup |A|/x',
}


def _up(v):
    """round a nonnegative bound upward by >= 1 ulp (the slab files hold nearest-rounded floats of exact arb bounds)"""
    import math
    return math.nextafter(math.nextafter(v, math.inf), math.inf) if isinstance(v, float) and v > 0 else v


def _round_up(r):
    for key in ('norms', 'norms_int', 't7'):
        r[key] = {k: _up(v) for k, v in r.get(key, {}).items()}
    r['scalars'] = {k: (_up(v) if k.startswith('sup_') else v) for k, v in r['scalars'].items()}
    return r


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--in', dest='inp', required=True)
    ap.add_argument('--first', type=int, default=5)
    ap.add_argument('--last', type=int, default=187)
    ap.add_argument('--json', default=os.path.join(HERE, 'residuals.json'))
    ap.add_argument('--jets', default=os.path.join(HERE, 'jets'))
    a = ap.parse_args()
    files = sorted(glob.glob(os.path.join(a.inp, 'slab_*.json')))
    recs = [_round_up(json.load(open(f))) for f in files]
    recs = [r for r in recs if a.first <= r['slab'] <= a.last]
    recs.sort(key=lambda r: r['slab'])
    have = [r['slab'] for r in recs]
    missing = [i for i in range(a.first, a.last + 1) if i not in have]
    summ = {}
    for nm, tg in TARGETS.items():
        vals = [(r['norms'].get(nm, 0.0), r['slab']) for r in recs]
        if not vals: continue
        mx, at = max(vals)
        fails = [s for v, s in vals if v > tg]
        summ[nm] = dict(max=mx, at_slab=at, target=tg, margin=(tg / mx if mx > 0 else None), fails=fails)
    other = sorted(set(k for r in recs for k in r['norms']) - set(TARGETS))
    for nm in other:
        vals = [(r['norms'].get(nm, 0.0), r['slab']) for r in recs]
        mx, at = max(vals); summ[nm] = dict(max=mx, at_slab=at)
    for nm in sorted(set(k for r in recs for k in r['scalars'] if k.startswith('sup_'))):
        mx, at = max((r['scalars'][nm], r['slab']) for r in recs); summ['scalar:' + nm] = dict(max=mx, at_slab=at)
    for nm in sorted(set(k for r in recs for k in r['t7'])):
        mx, at = max((r['t7'][nm], r['slab']) for r in recs); summ['t7:' + nm] = dict(max=mx, at_slab=at)
    cpu = [r.get('cpu', r['diag']['cpu_total']) for r in recs]
    out = dict(definitions=DEFS, targets=TARGETS, slabs_range=[a.first, a.last], missing=missing,
               cpu_total_s=sum(cpu), cpu_per_slab_mean_s=(sum(cpu) / len(cpu) if cpu else None),
               summary=summ, per_slab=recs)
    with open(a.json, 'w') as f: json.dump(out, f, indent=1)
    os.makedirs(a.jets, exist_ok=True)
    idx = []
    for r in recs:
        src = os.path.join(a.inp, 'jets', f'jets_{r["slab"]:03d}.npz')
        dst = os.path.join(a.jets, os.path.basename(src))
        if os.path.exists(src): shutil.copy2(src, dst)
        if os.path.exists(dst): idx.append(r['slab'])
    with open(os.path.join(a.jets, 'index.json'), 'w') as f:
        json.dump(dict(slabs=idx, format='jets_XXX.npz: for every name in _names an array (N_el, NC, 2) of [lo, hi] (outward-rounded '
                       'float64) over sub-cell c of element e (s in [-1 + 2c/NC, -1 + 2(c+1)/NC], x = x_e + h_e (1+s)/2) and the whole slab; '
                       '_meta = [slab, tau0, dtau, k, N_el, NC]. Names: x, k0..k5 (kappa^(m) in x), Q, theta, p, q, n, t1, V (true, incl. r), '
                       'varpi (true, incl. A), b1 (true), Vstar, varpistar (residual-free parts), r, A, r1, r2, r3, A1 = kappa r.'), f, indent=1)
    print(f'{len(recs)} slabs merged, missing {missing}; json -> {a.json}; jets -> {a.jets}')
    for nm in list(TARGETS)[:4] + [f'RF{j}_T5' for j in range(5)] + [f'RF{j}_opt4n' for j in range(5)]:
        if nm in summ: print(f"  {nm:12s} max {summ[nm]['max']:.3e} (slab {summ[nm]['at_slab']}) target {summ[nm]['target']:.0e} fails {summ[nm]['fails']}")


if __name__ == '__main__':
    main()
