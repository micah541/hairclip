"""T7: merge the per-unit outputs into t7/constants.json (keyed by slab, plus the certificate intervals).
Inputs (relative to this file): out/cert7_XX.json (C1+C2 local, reserved rank-3 budget), out_c2r3/c2r3_XXX.json (rank-3 variation
V_rank3 and C3 per slab), out_far/far_XXX.json (mu1, a_mu, mu2, Hardy), out_tip/tip_XXX.json (+ weight_consts.json),
out_farg/farg_XXX.json (W-Grams for C_gamma0/2), out_nl/nl_XXX.json (n4, n3, n2).  Missing units are reported as pending.
Checks: per certificate, sum of V_rank3 over the slabs on each side of tau_c <= c_R reserved;  lambda < 2 mu1 margin per slab.
usage: python3 merge.py [--lambda_max 0.2527]"""
import os, json, glob, argparse
HERE = os.path.dirname(os.path.abspath(__file__))


WORK = HERE


def load(pattern):
    out = {}
    for f in sorted(glob.glob(os.path.join(WORK, pattern))):
        d = json.load(open(f)); out[d.get('slab', d.get('id'))] = d
    return out


def slab_list(*units):
    """the slabs of P2: SLABS ('a-b') if given, else the whole store BASE (manifest.json), else every slab with any output"""
    if SLABS:
        a_, b_ = [int(v) for v in SLABS.split('-')]; return list(range(a_, b_ + 1))
    if BASE:
        p = BASE if os.path.isabs(BASE) else os.path.normpath(os.path.join(HERE, BASE))
        return sorted(r['slab'] for r in json.load(open(os.path.join(p, 'manifest.json'))))
    return sorted(set().union(*[set(u) for u in units[1:]], *[set(range(c['slabs'][0], c['slabs'][1] + 1)) for c in units[0].values()]))


def main(lam_max=0.2527):
    certs = load(os.path.join('out', 'cert7_*.json')); c2 = load(os.path.join('out_c2r3', 'c2r3_*.json'))
    far = load(os.path.join('out_far', 'far_*.json')); tip = load(os.path.join('out_tip', 'tip_*.json'))
    farg = load(os.path.join('out_farg_se' if SE else 'out_farg', 'farg_*.json')); nl = load(os.path.join('out_nl', 'nl_*.json'))
    SL = slab_list(certs, c2, far, tip, farg, nl)
    wc = os.path.join(WORK, 'out_tip', 'weight_consts.json'); wcon = json.load(open(wc)) if os.path.exists(wc) else {}
    C = dict(meta=dict(
        note='T7 constants for the P2 driver. All values are certified upper bounds (lower bounds for c_I, c_tau_c, mu1, mu2, nV_min) '
             'unless status says pending. theta_I = 0.025 (closure theta_0). Directions g = (p~, n~, V_s), V_s = -k\'\' + c p + b n; '
             'T6 residual corrections (V~ = V_s + r in g3, b1 += int k r) not included (sensitivity in T7_REPORT).',
        weight='W = (1+x^2)^-2 on [10.25, ell~] with T8 far-end flattening X = ell Phi(x/ell) (s0 = 1/2)',
        constraints='phi_1 = chi_32 p_L, phi_2 = chi_32 n_L (t7/phi_store.json, base7.late_base), phi_3 = chi_32 psi3^ (c4fun.PSI)',
        lambda_max_T8=lam_max, weight_consts=wcon), certificates=[], slabs={})
    for i in sorted(certs):
        d = certs[i]; lo, hi = d['slabs']; sc = d['slab_c']
        right = [s for s in range(sc, hi + 1)]; left = [s for s in range(lo, sc)]
        have = all(s in c2 for s in range(lo, hi + 1))
        Vr = sum(c2[s]['V_rank3'] for s in right if s in c2); Vl = sum(c2[s]['V_rank3'] for s in left if s in c2)
        status = ('verified' if max(Vr, Vl) <= d['c_R_reserved'] else 'FAILED') if have else 'pending (c2r3)'
        C['certificates'].append(dict(id=i, slabs=[lo, hi], tau=d['tau_I'], tau_c=d['tau_c'], slab_c=sc, k_c=d['k_c'],
                                      theta=d['theta'], c_I=d['c_I'], c_tau_c=d['c_tau_c'], c_R_reserved=d['c_R_reserved'],
                                      rank3_V=dict(right=Vr, left=Vl, status=status), mu0_max=d['mu0_max'], mu1_max=d['mu1_max'],
                                      C3_tau_c=d.get('C3_tau_c'), c_phi=d.get('c_phi_normalised'), n=d['n'], h=d['h']))
        for s in range(lo, hi + 1):
            C['slabs'].setdefault(str(s), {})['cert'] = dict(id=i, c_I=d['c_I'], theta=d['theta'], status=status)
    for s in SL:
        e = C['slabs'].setdefault(str(s), {})
        src = c2.get(s) or far.get(s) or tip.get(s)
        if src: e['tau'], e['k'] = src['tau'], src['k']
        if s in c2: e['C3'] = c2[s]['C3']; e['V_rank3'] = c2[s]['V_rank3']
        if s in far:
            f = far[s]; e['E_I'] = dict(mu1=f['mu1'], a_mu=f['a_mu'], mu_f=f['muf'], split_ok=f['split_ok'], mu1_raw=f['mu1_raw'],
                                       hardy=f['hardy'], mu2=f['mu2'], lam_lt_2mu1_margin=2 * f['mu1'] / lam_max - 1)
        if s in tip:
            t = tip[s]; e['E2'] = dict(G22=t['G22'], G12=t['G12'], G02=t['G02'])
            e['strip'] = dict(sB=t['sB'], sA=t['sA'], sA0=t['sA0'], IA=t['IA'], d3=t['d3'], d2=t['d2'], d1=t['d1'], d0=t['d0'])
            e['ramp'] = dict(a_f=t['a_f'] + (far[s]['a_mu'] if s in far else 0.0), b_f=t['b_f'], c3=t['c3'], c2=t['c2'], c1=t['c1'],
                             c0=t['c0'], a_f_without_split=t['a_f'])
        if s in farg and s in c2:
            g = farg[s]; Cx = c2[s]['C3']['Cx']; Cxi = c2[s]['C3'].get('Cx_i')
            comp = lambda key: sum(a_ * b_ for a_, b_ in zip(Cxi, g[key])) if Cxi and key in g else float('inf')
            e['modulation'] = dict(Cgamma0=min(Cx * g['sqrt_lmax_G0'], comp('sqrt_diag_G0')), Cgamma2=min(Cx * g['sqrt_lmax_G2'], comp('sqrt_diag_G2')),
                                   Gx=min(Cx * g['sqrt_lmax_G2E'], comp('sqrt_diag_G2E')),
                                   note='|sum x_i f_i| <= min(C_x sqrt(tr Gamma), sum_i C_x,i ||f_i||), normalised directions')
        if s in farg and s not in c2 and 'cert' in e:
            cert = certs[e['cert']['id']]; Cx = (cert.get('C3_tau_c') or {}).get('Cx'); g = farg[s]
            if Cx: e['modulation'] = dict(Cgamma0=Cx * g['sqrt_lmax_G0'], Cgamma2=Cx * g['sqrt_lmax_G2'], Gx=Cx * g['sqrt_lmax_G2E'],
                                          note='PROVISIONAL: C_x at tau_c of the certificate (c2r3 batch pending), times sqrt(tr Gamma)')
        if s in nl: e['nonlinear'] = dict(tip=nl[s]['tip'], far=nl[s]['far'], far_max=nl[s]['far_max'], farB_max=nl[s].get('farB_max'),
                                          box=nl[s]['box'], se=nl[s].get('se', False))
        if s in farg:                                               # raw Gram values (the P2 driver's (Y) uses them)
            g = farg[s]; e['farg'] = {k: g[k] for k in ('sqrt_lmax_G0', 'sqrt_lmax_G2', 'sqrt_lmax_G2E', 'sqrt_diag_G0', 'sqrt_diag_G2',
                                                         'sqrt_diag_G2E', 'norm_chi_g_min', 'nsub') if k in g}
            if 'se_correction' in g: e['farg']['se_increments'] = {k: g['se_correction'][k] for k in ('inc_G0', 'inc_G2')}
        if s in far: e['E_I']['se'] = far[s].get('se', False)
    # closure (certified rho_c per interval, common mu_e) and the e-bound K_e
    cc = os.path.join(WORK, 'closure_cert.json')
    if os.path.exists(cc):
        R = json.load(open(cc))
        C['meta']['closure'] = dict(mu_e=R['mu_e'], rho_c_min=R['rho_c_min'], all_certified=R['all_certified'],
                                    rho_c={str(o['id']): o['rho_c'] for o in R['per_interval']}, kappa_D=0.5, theta_0=0.025,
                                    note='closure_cert.py: Q >= rho_c Lambda_e on all shapes (AM-GM relaxation + arb LDL), T5 2.4 form, (E2) AB form')
        for cert in C['certificates']: cert['rho_c'] = C['meta']['closure']['rho_c'].get(str(cert['id']))
        import provenance as PV
        C['meta']['closure']['args'] = dict(file=PV.entry(cc), closure_cert_args=R.get('args'))
    NT = None
    ccn, nlt = os.path.join(WORK, 'closure_cert_nl.json'), os.path.join(WORK, 'nl_terms.json')
    if os.path.exists(ccn) and os.path.exists(nlt) and 'closure' in C['meta']:
        nonlinear_closure(C, json.load(open(ccn)), NT := json.load(open(nlt)), ccn, nlt)
    from flint import arb, ctx
    ctx.prec = 128
    eps_, c1_ = arb('1e-5'), arb('0.5804'); xa = arb(14) + arb('2.2490776515'); xr = arb('34.25')
    phi = lambda x: (1 + x * x) ** (arb(2) / 3) / 4
    dphi = xr * (1 + xr * xr) ** (-arb(1) / 3) / 3                     # sup_[0,34.25] phi' (phi' increasing)
    s_ = (1 / (2 * eps_.sqrt()) + c1_).sqrt()
    Ke = (phi(xa) / 2).exp() * (2 * (s_ + dphi / 2)).sqrt()
    C['meta']['K_e'] = dict(value=float(Ke.upper()), x_a=float(xa.mid()),
                            definition='|e(x)| <= K_e Ybar for x <= x_a = 14 + X0 (T8 level-0 ramp start): w = e^{-phi/2} u, w(34.25) = 0, '
                                       'w^2 <= 2|w||w\'|, |w\'| <= W1 + sup|phi\'|/2 A, W1^2 <= AB + c1 A^2 <= Lambda (1/(2 sqrt eps) + c1)')
    # sensitivity of the certificates to T6's corrections (V~ = V_s + r in g_3, b1 += int k r), from t6/residuals.json:
    #   per slab  2 mbar C_x ||chi r||_rho / n_min  +  1/2 sup_[0,34.25] (|k r| + |phi' int k r|),  per interval the max over its slabs
    p6 = os.path.join(HERE, '..', 't6', 'residuals.json')
    if os.path.exists(p6):
        R6 = {r['slab']: r for r in json.load(open(p6))['per_slab']}
        sens = {}
        for s in SL:
            if s in R6 and s in c2:
                r, c3 = R6[s], c2[s]['C3']
                sens[s] = (2 * c3['mbar'] * c3['Cx_full'] * r['norms']['Rs_chi'] / c3['nV_min'] + 0.5 * (r['t7']['supA1_tip'] + r['t7']['supArho_tip'])) * (1 + 1e-12)
        if len(sens) == len(SL):
            per = {}
            for cert in C['certificates']:
                lo, hi = cert['slabs']; v = max(sens[s] for s in range(lo, hi + 1)); per[str(cert['id'])] = v
                cert['T6_sensitivity'] = v; cert['c_I_with_T6'] = cert['c_I'] - v
            C['meta']['T6_sensitivity'] = dict(max=max(per.values()), per_interval=per, source=os.path.relpath(p6, HERE),
                                               formula='2 mbar Cx_full ||chi r||_rho / nV_min + (sup|k r| + sup|phi\' int k r|)/2 on [0, 34.25]')
        else:
            C['meta']['T6_sensitivity'] = dict(status=f'pending: t6/residuals.json or c2r3 lacks {len(SL) - len(sens)} slabs')
    # summary
    m = [e['E_I']['lam_lt_2mu1_margin'] for e in C['slabs'].values() if 'E_I' in e]
    C['meta']['summary'] = dict(n_cert=len(C['certificates']), min_c_I=min((c['c_I'] for c in C['certificates']), default=None),
                                min_c_I_with_T6=min((c['c_I_with_T6'] for c in C['certificates'] if 'c_I_with_T6' in c), default=None),
                                rank3_status={c['id']: c['rank3_V']['status'] for c in C['certificates']},
                                min_lam_margin=min(m) if m else None,
                                slabs=[SL[0], SL[-1]], uncovered=[s for s in SL if 'cert' not in C['slabs'][str(s)]],
                                pending=dict(c2r3=[s for s in SL if s not in c2], far=[s for s in SL if s not in far],
                                             tip=[s for s in SL if s not in tip], farg=[s for s in SL if s not in farg],
                                             nl=[s for s in SL if s not in nl]))
    C['meta']['work_dir'] = os.path.relpath(WORK, HERE); C['meta']['base_dir'] = BASE
    if SE:
        import provenance as PV
        sw = os.path.join(HERE, 'se_widths.json'); SW = json.load(open(sw))
        chk = None; so = os.path.join(HERE, 'se_widths_secorr.json')
        if os.path.exists(so):                 # far7 (out_far) was computed with the se_corr widths: valid if they contain the direct ones
            O = json.load(open(so))['per_slab']; N = SW['per_slab']
            rat = max(max(N[s]['k'][m_] / O[s]['k'][m_] for m_ in range(4)) for s in N if s in O and 'source' not in N[s])
            ratQ = max(N[s]['Q'] / O[s]['Q'] for s in N if s in O and 'source' not in N[s])
            chk = dict(far7_inputs_k0_k3_Q_max_ratio_direct_over_secorr=max(rat, ratQ), far7_valid=max(rat, ratQ) <= 1,
                       note='far7 widens only k^(0..3) and Q; out_far (se_corr widths) contains the direct widths, so it holds for the corrected base',
                       secorr_widths=PV.entry(so))
        C['meta']['se_correction'] = dict(widths=PV.entry(sw), widths_args=SW.get('args'), far7_check=chk,
            direct='review B1: p2/t6c (T6 engine on the corrected base); nl7 rerun and farg increments recomputed with the direct widths',
            note='far constants (far7 mu1/mu2/a_mu/Hardy, nl7 far boxes, farg Gamma0/Gamma2) hold for the (SE)-corrected base '
                 '(STEP4_RESULT 2(ii)): jets widened on x >= 0.6 ell_lo, farg by certified increments; all tip-region constants '
                 '(certificates, C3, E2, strip, ramp, Gamma2E, tip n_j) are unchanged since the corrected base equals the polished one there')
    if nl and all('m23' in v['tip'] for v in nl.values()): C['meta']['nonlinear_tip'] = tip_lemma(C, nl, c_e1=(NT['tip_box']['c_e1'] if NT else None))
    json.dump(C, open(os.path.join(WORK, 'constants.json'), 'w'), indent=1)
    if PUBLISH and WORK != HERE:
        import shutil; shutil.copy(os.path.join(WORK, 'constants.json'), os.path.join(HERE, 'constants.json'))
    S = C['meta']['summary']
    print('certificates', S['n_cert'], 'min c_I', S['min_c_I'], '(with T6 corrections', S['min_c_I_with_T6'], ') min lambda<2mu1 margin', S['min_lam_margin'])
    print('pending counts', {k: len(v) for k, v in S['pending'].items()}, 'slabs', S['slabs'], 'not covered by a certificate', len(S['uncovered']))
    return C


BASE = None
SLABS = None
PUBLISH = False
SE = False


def nonlinear_closure(C, RN, NT, ccn=None, nlt=None):
    """T7 Lemma 4c (nlterms7.py + closure_cert.py --nl): the certified rho_c with the nonlinear terms of (Y) included.  The
    reserve Ds and its split are unchanged (kappa_D = 1/2 for Lemma J, the other half for the residual Young terms of the
    P2 (Y) comparison), and no separate allowance is needed (DN = D4CUB = 0).  Valid on the bootstrap box Lambda_e <= Ybox^2
    with T8's far box (n3.json) and Lemma 4b's tip box.  Lemma 4b's tip box and the (M) bound are published alongside."""
    import provenance as PV
    files = dict(nl_terms=PV.entry(nlt), closure_cert_nl=PV.entry(ccn)) if nlt else {}
    args = dict(files=files, nl_terms_args=NT.get('args'), closure_cert_nl_args=RN.get('args'))
    lin = C['meta']['closure']; rn = {str(o['id']): o for o in RN['per_interval']}
    assert RN['all_certified'] and abs(RN['mu_e'] - lin['mu_e']) < 1e-15
    C['meta']['closure_linear'] = dict(lin)
    rho = {k: (rn[k]['rho_c'] if k in rn else None) for k in lin['rho_c']}
    eps, epsI, mu, th0 = 1e-5, 3e-5, RN['mu_e'], 0.025
    C['meta']['closure'] = dict(
        mu_e=mu, mode='nonlinear (T7 Lemma 4c)', rho_c=rho, rho_c_min=min(v for v in rho.values() if v is not None),
        all_certified=True, kappa_D=0.5, theta_0=th0, Ybox=NT['Ybox'], DN=0.0, D4CUB=0.0,
        reserve='Ds = theta0/4 B^2 + eps/8 D4^2 + mu_e/8 V2^2 + mu_e eps_I/8 V4^2 (unchanged): Q(z) - Ds - NL(z) >= rho_c Lambda_e',
        shares_ycomp=dict(B=th0 / 8, D4=eps / 16, V2=mu / 16, V4=mu * epsI / 16), shares_lemmaJ='the other half (kappa_D = 1/2)',
        not_certified={k: 'outside P2 (SLAB0 = 2): T8 n3.json has no far majorants for these slabs' for k, v in rho.items() if v is None},
        note='closure_cert.py --nl nl_terms.json: Q - Ds - NL >= rho_c Lambda_e on all shapes (AM-GM relaxation + arb LDL); '
             'NL = nonlinear terms N, Q and the nonlinear/E_ij part of (M), reduced to quadratic terms on the box (T7_REPORT 4c)',
        args=args)
    for cert in C['certificates']:
        cert['rho_c_linear'] = lin['rho_c'].get(str(cert['id'])); cert['rho_c'] = rho.get(str(cert['id']))
    per = {}
    for r in NT['intervals']:
        k = str(r['id'])
        per[k] = dict(rho_c=rho[k], rho_c_linear=lin['rho_c'][k], loss=(lin['rho_c'][k] - rho[k]) if rho[k] is not None else None,
                      top_I0_V2sq=r['top']['I0_V2sq'], top_I2_V4sq=r['top']['I2_V4sq'], sup_At_plus=r['top']['sup_Atp'],
                      n_terms=len(r['terms']), max_dissipation_degree=r['max_dissipation_degree'], box=r['box'],
                      parts={p: dict(n_terms=v['n_terms'], max_coef=v['max_coef']) for p, v in r['parts'].items()})
    C['meta']['nonlinear_closure'] = dict(
        Ybox=NT['Ybox'], terms='t7/nl_terms.json (closure units; closure_cert.py --nl)', far_majorants=NT['n3'],
        validity='Lambda_e <= Ybox^2 on the interval, T8 box B\' (n3.json run at Y = 1.7e-12), Lemma 4b tip box (eta1 = c_e1 Ybox <= 1e-5)',
        n_dependence='tip: nl7 n4, n3_0, n2_0, n21, m2, m23 (T5 box, valid for eta1 <= 1e-5) times eta1, S_W (both linear in Ybox); '
                     'A31 = 10 e\'/J^6 and At analytic; ramp/middle: D(N) monomial coefficients (TD2 enclosures, nlsym/nlforms_gen.py) '
                     'with T6 jets; far core x >= 30.25: T8 n3.json per-range majorants (Atp, A, A1, B0, B1, dCe, beta, G) per slab',
        per_interval=per, args=args)
    C['meta']['tip_box'] = dict(NT['tip_box'], args=dict(files=files, nl_terms_args=NT.get('args')), eta1_at_Ybox=NT['tip_box']['c_e1'] * NT['Ybox'], Ybox=NT['Ybox'],
                                asserted_against_T5_eta1=1e-5, holds=NT['tip_box']['c_e1'] * NT['Ybox'] <= 1e-5,
                                Ybar_max_for_T5_box=1e-5 / NT['tip_box']['c_e1'], lemma='T7_REPORT 4b')
    xs = [r['M'] for r in NT['intervals']]
    C['meta']['modulation_nl'] = dict(
        statement='(M): |x^ - x^_lin - x^_res| <= mbar |b_NL| + mbar|E|/(1 - mbar|E|) (C_x A + x_r + mbar |b_NL|) (T7_REPORT 4c(M))',
        xNL_box_max=max(x['xNL_box'] for x in xs), slack_min_at_Ybar_1e12=min(x['slack_1e12'] for x in xs),
        ratio_max=max(x['xNL_box'] / x['slack_1e12'] for x in xs), E_box_max=max(x['E_box'] for x in xs),
        mbarE_max=max(x['E_box'] * x['mbar'] for x in xs), xr_increment_max=max(x['xr_increment'] for x in xs),
        Ycert_used=NT['Ycert'], Ybox=NT['Ybox'], args=dict(files=files, nl_terms_args=NT.get('args')),
        per_interval={str(r['id']): r['M'] for r in NT['intervals']})


def tip_lemma(C, nl, Yref='1e-12', c_e1=None):
    """Gap 9 (T7_REPORT 4a): the tip constants hold with no condition on e''; sup|e''| on [0, 34.25] is controlled by
    S^2 <= S_W^2 + 2 K_reg B W3,  S_W^2 <= (V1^2/18 + 2 V1 V2)/W(34.25)  (E_I where zeta = 1; K_reg = 1/rho(16.25)).
    Cost of the one new term m23|e''|^3 in T5 (E0), (E2) at the verified box for Ybar = Yref (|e|,|e'| <= 2e3 Ybar/sqrt(mu_e),
    Lambda_e <= Ybar^2), split by Young into a rate loss and a fraction of eps D4^2 (formulas in T7_REPORT 4a)."""
    from flint import arb, ctx
    ctx.prec = 128
    tips = [v['tip'] for v in nl.values()]
    mx = lambda k: max(t[k] for t in tips)
    m = arb(repr(mx('m23')))
    eps, epsI, c1 = arb('1e-5'), arb('3e-5'), arb('0.5804')
    c1W = arb(repr(C['meta']['weight_consts'].get('c1W', 0.0943)))
    muE = arb(repr(C['meta'].get('closure', {}).get('mu_e', 0.036)))
    Phi1, Phi2 = arb('1.083'), arb('1.3')
    rho = lambda x: (-(1 + x * x) ** (arb(2) / 3) / 4).exp(); W = lambda x: (1 + x * x) ** -2
    Kreg = 1 / rho(arb('16.25')); Wmin = W(arb('34.25')); sig = (rho(arb('31.249')) / Wmin).sqrt()
    kA = arb(1) / 2 + c1.sqrt(); kV1 = arb(1) / 2 + c1W.sqrt(); c2 = 2 * Phi1 * kA + Phi2
    Y = arb(Yref); eta1 = (arb(repr(c_e1)) * Y) if c_e1 else 2000 * Y / muE.sqrt()      # Lemma 4b (c_e1) when published
    eB = Y / eps.sqrt(); eV0 = Y / muE.sqrt(); eV2 = Y / (muE * epsI).sqrt()
    eV1 = kV1 * eV0 + eV2 / 2; cV1 = 2 * kV1 ** 2 / muE + 1 / (2 * muE * epsI)
    SW = ((eV1 ** 2 / 18 + 2 * eV1 * eV2) / Wmin).sqrt()
    eg = eB + sig * eV1; eQ2 = eB ** 2 + sig ** 2 * eV1 ** 2; cQ = 1 / eps + sig ** 2 * cV1
    bet = arb('1e-9')                                                  # Young weights (fractions of eps D4^2)
    d1 = eta1 * m * (SW + (2 * Kreg * kA).sqrt() * eB) * cQ           # (E0): S <= S_W + sqrt(2 K kA) B + sqrt(K B D4)
    a = eta1 * m * (Kreg * eB).sqrt(); lam = (4 * bet * eps) ** (arb(1) / 4)
    d2 = arb(3) / 4 * (a / lam) ** (arb(4) / 3) * eQ2 ** (arb(1) / 3) * cQ
    P = SW ** 2 + 2 * Kreg * kA * eB ** 2                             # (E2): S^2 <= P + K B D4
    ba = m * Kreg * (1 + Phi1) * eB * eg
    db = eps * (m * P * (1 + Phi1)) ** 2 * (2 / eps + 2 * sig ** 2 * cV1) / (4 * bet)
    dc = eps * m * P * c2 * ((1 + sig / 2) / eps + sig / 2 * cV1)
    dd = (m * Kreg * c2 * eB * eg) ** 2 / (4 * bet)
    U = lambda v: float(v.upper())
    return dict(constants=dict(n4=mx('n4'), n3=mx('n3'), n2=mx('n2'), m23=mx('m23'), m3=mx('m3'), m2=mx('m2'), n21=mx('n21'),
                               n2_0=mx('n2_0'), n3_0=mx('n3_0'), box='|e|, |e\'| <= 1e-5 on [0, 34.25]; e\'\' free'),
                K_reg=U(Kreg), W_34=U(Wmin), sigma_strip=U(sig),
                at_Ybar=dict(Ybar=float(Y.mid()), eta1=U(eta1), B=U(eB), V1=U(eV1), V2=U(eV2), S_W=U(SW),
                             cubic_rate_loss=U(d1 + d2 + db + dc + dd), cubic_D4_fraction=U(3 * bet + ba),
                             parts=dict(E0_rate=U(d1 + d2), E2_rate=U(db + dc + dd), E2_D4_fraction=U(ba))),
                eta1_source=('Lemma 4b: eta1 = c_e1 Ybar' if c_e1 else 'T5 2e3 Ybar/sqrt(mu_e) (superseded)'),
                superseded_by='T7_REPORT 4c: the cubic term is inside the nonlinear closure (DN = D4CUB = 0)' if c_e1 else None,
                lemma='T7_REPORT 4a')

if __name__ == '__main__':
    ap = argparse.ArgumentParser(); ap.add_argument('--lambda_max', type=float, default=0.2527)
    ap.add_argument('--work', default='.', help='directory holding out/, out_c2r3/, out_far/, ... (relative to t7/)')
    ap.add_argument('--base', default=None, help='polished store used (its slabs = the slabs of P2)'); ap.add_argument('--publish', action='store_true')
    ap.add_argument('--slabs', default=None, help='a-b: restrict to slabs a..b (v1 store: 5-187)')
    ap.add_argument('--se', action='store_true', help='far constants of the (SE)-corrected base (out_farg_se, far7/nl7 --se)')
    A = ap.parse_args(); WORK = A.work if os.path.isabs(A.work) else os.path.normpath(os.path.join(HERE, A.work))
    BASE = A.base; PUBLISH = A.publish; SLABS = A.slabs; SE = A.se
    main(A.lambda_max)
