#!/usr/bin/env python3
"""Task 2 (r2), last stage: out/A_interface.json (CONTRACTS.md Task 2 'A_interface.json', amended by CONTRACTS_r2.md A6).
Every number is exported as an arb string and as an outward-rounded double ('up' for upper bounds, 'lo' for lower bounds), with a label.
Reads only outputs of this task (paths below; all pinned in the args block).  Run: PYTHONHASHSEED=0 nice -n 19 python3 -B rigor/export_t2.py"""
import os, sys, json
sys.dont_write_bytecode = True
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import arlib as L
from arlib import A, Q_, arb
O = L.OUT
F = dict(asm_ref=O + '/asm_design_N2ref.json', asm_ref4=O + '/asm_design_N2ref4.json', asm_0=O + '/asm_design_N20.json',
         asm_lo=O + '/asm_lo_N2ref.json', asm_hi=O + '/asm_hi_N2ref.json', asm_probe=O + '/asm_probe_N2ref.json', asm_pm3=O + '/asm_design_pm3ibp.json',
         hjw=O + '/hjw_rep.json', ic_design=O + '/ic_design.json', ic_lo=O + '/ic_lo.json', ic_hi=O + '/ic_hi.json', ic_probe=O + '/ic_probe.json',
         weights=O + '/weights_t2.json', weights_fb=O + '/weights_t2_fb.json', outer=O + '/outer_L2.json', ramp=O + '/ramp_rep.json', signs=O + '/a3_signs.json')
BETA_A = '5/4'
N2REF = ['0.1', '0.2', '5.0']

def num(x, kind='up', label='cert'):
    x = x if isinstance(x, arb) else A(Q_(str(x)))
    return dict(arb=x.str(15, radius=True), double=(L.fup(x) if kind == 'up' else L.flo(x)), bound=kind, label=label)
def vec(xs, kind='up', label='cert'): return [num(x, kind, label) for x in xs]

def main():
    D = {k: L.load(v) for k, v in F.items()}
    out = dict(provisional=False, labels='cert = arb, this task; est = reference value (not a bound); hyp = hypothesis (not certified here)',
               beta_A=BETA_A, N2_ref=dict(values=N2REF, label='est'), res={})
    W, OL, HJ = D['weights'], D['outer'], D['hjw']
    kap = [[OL['res'][str(j)]['kappa'][k] for k in range(3)] for j in range(3)]
    for st in ('tight', 'v31'):
        R, R4, R0 = D['asm_ref']['res'][st], D['asm_ref4']['res'][st], D['asm_0']['res'][st]
        gp = R['gp']; EA0 = R0['E']; EAr = R['E']; EAr4 = R4['E']
        acore = (1 + A(Q_(11, 100))**2)**-2; r = [A(Q_(x)) for x in W['r_export'][:2]]
        rhoJ = R['export']['0']['rhoJ']
        # rule check: EA0 + (o_j + 2 sum_{k<j} Ghat_jk eps28 (E_k - EA0_k))/gp_j reproduces the assembled EA_ref (coupling included)
        GH = {j: {int(k): v for k, v in W['levels'][str(j)]['core']['Ghat'].items()} for j in range(3)}; e28 = arb(28)**-4
        rule = []
        for j in range(3):
            o = sum((kap[j][k]*A(Q_(N2REF[k])) for k in range(3)), arb(0))
            cpl = sum((GH[j][k]*e28*(rule[k] - EA0[k]) for k in range(j)), arb(0))
            rule.append(EA0[j] + (o + 2*cpl)/gp[j])
        rule_ok = all((EAr[j] - rule[j]).upper() <= 1e-15*max(1.0, float(EAr[j].upper())) for j in range(3))
        exp = R['export']; exp0 = R0['export']
        Dsrc = [R['info'][k]['D'] for k in range(2)]; F2 = [R['info'][k]['F2'] for k in range(2)]; coup = [R['info'][k]['coup_rate'] for k in range(2)]
        fl = R['flux']
        res = dict(
            EA0=vec(EA0), kappa_o=[vec(kap[j]) for j in range(3)], gp=vec(gp, 'lo'),
            EA_of_N2='EA_j = EA0_j + (o_j + 2 sum_{k<j} Ghat_core[j][k] 28^-4 (EA_k - EA0_k))/gp_j,  o_j = sum_k kappa_o[j][k] Nbar_k^2 (the coupling term is < 1e-4 of o_j/gp_j)',
            EA_ref=vec(EAr), EA_ref4=vec(EAr4), N2_ref=N2REF, rule_reproduces_EA_ref=rule_ok,
            R5prime_addition=dict(at_N2ref=vec([EAr[j] - EA0[j] for j in range(3)]), at_4N2ref=vec([EAr4[j] - EA0[j] for j in range(3)]),
                                  o0_over_gp0=dict(N2ref=num(sum((kap[0][k]*A(Q_(N2REF[k])) for k in range(3)), arb(0))/gp[0]), N2ref4=num(4*sum((kap[0][k]*A(Q_(N2REF[k])) for k in range(3)), arb(0))/gp[0]))),
            Dsrc=vec(Dsrc), r_export=W['r_export'][:2], amin_c=num(acore, 'lo'), rhoJ=num(rhoJ), F2=vec(F2), coup=vec(coup),
            C0_rule='C0_k = (EA_k + Dsrc_k)/(r_k amin_c)', C1_rule='C1_k = C1_0_k + o_k/(r_k amin_c)',
            C0_ref=vec([exp[str(k)]['C0'] for k in range(2)]), C1_0=vec([exp0[str(k)]['C1'] for k in range(2)]), C1_ref=vec([exp[str(k)]['C1'] for k in range(2)]),
            IC=dict(setting='design', H2=D['ic_design']['E_tau1']['%s/design' % st]['H2'], H2f=D['ic_design']['H2f'], H2f_label='hyp',
                    E_tau1=vec(D['ic_design']['E_tau1']['%s/design' % st]['E']), rule='E^A_j(tau_1) <= H2_j + H2f_j', step4_H2_range=D['ic_design']['step4_H2_range']),
            core=dict(S0v=num(R['S0']), S1v=num(R['S1']), margins={k: num(v, 'lo') for k, v in R['margins'].items()},
                      margin_thresholds={'S0v': 0.02, 'S0v_vs_0.0153': 0.0153, 'M0core': 0.11, 'M1core': 0.52, 'Ghat22': 153664, 'ramp M1': 0.52, 'ramp z1': 0.52, 'ramp M0': 0.05},
                      S0v_margin_vs_0p0153=num(R['margins']['S0v_vs_0.0153'], 'lo'), beta_A=BETA_A,
                      K0_at_E=num(fl['K0_at_E']), K2_at_E=num(fl['K2_at_E']),
                      K0_at_betaE={b: num(v['K0']) for b, v in fl['beta'].items()}, K2_at_betaE={b: num(v['K2']) for b, v in fl['beta'].items()},
                      flux_margin={b: dict(FLUX0_over_K0=num(v['margin_FLUX0'], 'lo'), FLUX2_over_K2=num(v['margin_FLUX2'], 'lo')) for b, v in fl['beta'].items()},
                      FLUX0=1.64e-20, FLUX2=8.0e-20, FLUX0_needed=num(fl['FLUX0_needed']), FLUX2_needed=num(fl['FLUX2_needed']),
                      FLUX_needed_note='4 K(EA): the verifier inputs that would give x2 at beta_A = 2'),
            v31_only=dict(note='(B)-exit claims of v3.1 (lemmaS_F4B); not consumed by STATEMENTS_r2.tex (replaced by Lemma U)',
                          S0_B_yE50=num(R['B50']['S0']), S0V=0.0153, S0_B_le_S0V=bool((R['B50']['S0'] - A(Q_('0.0153'))).upper() <= 0),
                          claims_failing=[list(c) for c in R['claims_failing']]),
            parts={str(j): {k: num(v) for k, v in R['parts'][j].items()} for j in range(3)})
        out['res'][st] = res
    # A4 / A4b: the four (H2f) settings at the tightened (H2)
    hs = {}
    for nm in ('design', 'lo', 'hi', 'probe'):
        R = D['asm_' + ('ref' if nm == 'design' else nm)]['res']['tight']; fl = R['flux']
        hs[nm] = dict(H2f=D['ic_' + nm]['H2f'], E_tau1=vec(D['ic_' + nm]['E_tau1']['tight/' + nm]['E']), EA=vec(R['E']),
                      EA0_parts={k: num(R['parts'][0][k]) for k in ('handoff', 'inflow', 'ramp_comm', 'top', 'cls', 'outer')},
                      K0_at_E=num(fl['K0_at_E']), K2_at_E=num(fl['K2_at_E']),
                      beta={b: dict(K0=num(v['K0']), K2=num(v['K2']), FLUX0_over_K0=num(v['margin_FLUX0'], 'lo'), FLUX2_over_K2=num(v['margin_FLUX2'], 'lo')) for b, v in fl['beta'].items()},
                      FLUX0_needed=num(fl['FLUX0_needed']), FLUX2_needed=num(fl['FLUX2_needed']),
                      S0v_vs_0p0153=num(R['margins']['S0v_vs_0.0153'], 'lo'), M1core=num(R['margins']['M1core'], 'lo'),
                      S0_B_yE50=num(R['B50']['S0']))
    out['H2f_settings'] = hs
    out['HJW'] = dict(W1_ok=HJ['W1_ok'], W2_ok=HJ['W2_ok'], W3_ok=HJ['W3_ok'], W1_note=HJ['W1_note'], W2_note=HJ['W2_note'], W3_note=HJ['W3_note'],
                      kappa_J=num(Q_(repr(HJ['kappa_J']))), kappa_J_W1=num(Q_(repr(HJ['kappa_J_W1']))), kappa_J_W3=num(Q_(repr(HJ['kappa_J_W3']))),
                      F0=HJ['F0'], tau_march=HJ['tau_march'], X_const=HJ['X_const'], worst_cell=HJ['worst_cell'], verifier_file=HJ['verifier_file'], verifier_sha256=HJ['verifier_sha256'],
                      statement='Lemma 4 of STATEMENTS_r2.tex (F0 in c_k): int over any window of V_m^2 <= J_m (1 + kappa_J) + s_m T, m = 0,1,2,4')
    rp = D['ramp']['pointterm']; Rt, Rp = D['asm_ref']['res']['tight'], D['asm_pm3']['res']['tight']
    out['pointterm'] = dict(gc_28=num(rp['gc_28']), routes=rp['routes'], signs=D['signs']['res'],
                            dE=vec([Rt['E'][j] - Rp['E'][j] for j in range(3)]), dE_note='E (direct data route for the m=3 monomial of (j,j+1), levels 0,1) minus E (phase-2a IBP route, point term missing); design, tight')
    lv = W['levels']
    out['interface28'] = dict(
        convention_found='v3.1/phase 2a: interpolation, L4 and {2,3,3} identities taken separately on the ramp [14,28] (box (0.05,0.9)) and the core [28,350] '
                         '(box (0.11,0.52)), per-region constants, boundary terms at y = 28 dropped; not covered by Lemma L(i). L(ii) with Ytilde falling from '
                         'Y_ramp to Y_core on [28,56] fails (core Ghat_22 at 28 would be 2.4e5 > 153664).',
        repair='ramp coefficient box := core box (verifier: ramp |y f_y| <= 0.4396, |z1| <= 0.4328 < 0.52); one C^{1,1} Ytilde on [14,700] (rises on [40,50] and [175,250]); '
               'L4/S33 localised by partitions of unity on [50,60] (S33 weight) and [250,350] (P1|P2), one constant per product per piece (L(iii)); weights_t2.py',
        Ghat_core=[vec([lv[str(j)]['core']['Ghat'][str(k)] for k in range(j + 1)]) for j in range(3)],
        gp_before=dict(values=[0.24292, 0.24139, 0.12087], source='phase2a A_REPAIR.md (lower bounds)', label='cert'), gp_after=vec(D['asm_ref']['res']['tight']['gp'], 'lo'),
        Ghat22_before=79364.6, Ghat22_after=num(lv['2']['core']['Ghat']['2']))
    def segs(Wd, box):
        o = {}
        for j in range(3):
            l = Wd['levels'][str(j)]
            o[str(j)] = dict(box=box, seg_60_175=dict(ok=l['seg']['check']['ok'], rows=l['seg']['check']['rows'], retained=l['seg']['check']['retained_abs_over_amin_c']),
                             rise_40_50=dict(ok=l['transition1']['rise']['ok'], worst=l['transition1']['rise']['worst']),
                             S33rise_50_60=dict(ok=l['transition1']['partition']['ok'], worst=l['transition1']['partition']['worst']),
                             rise_175_250=dict(ok=l['transition2']['rise']['ok'], worst=l['transition2']['rise']['worst']),
                             chi_250_350=dict(ok=l['transition2']['partition']['ok'], worst=l['transition2']['partition']['worst']),
                             O_350_700=dict(ok=l['O']['check']['ok'], rows=l['O']['check']['rows'], Bc=l['O']['Bc'], delta_o=l['O']['delta_o']),
                             Ytilde={k: num(v) for k, v in l['Ytilde'].items()})
        return dict(all_ok=Wd['all_ok'], levels=o, U='<= 1/20 on [28,60], <= 0.08 on (60,350], <= 0.162 on O', z1=Wd['boxes'])
    out['segments'] = {'(0.15,0.52)': segs(W, '(0.15,0.52)'), '(0.20,0.60)': segs(D['weights_fb'], '(0.20,0.60)')}
    out['classes_outer'] = {str(j): [dict(I=e['I'], ft=e['ft'], tool=e['tool'], routes=e['routes'], mq=e['mq'], cost={k: num(v) for k, v in e.get('cost', {}).items()})
                                     for e in OL['res'][str(j)]['classes']] for j in range(3)}
    out['outer_other'] = {str(j): dict(interp_remainder=num(OL['res'][str(j)]['interp_remainder']['cost']), outer_foot_N2=num(OL['res'][str(j)]['outer_foot']['N2']),
                                       chain=({k: num(v) for k, v in OL['res'][str(j)]['chain'].items()} if OL['res'][str(j)]['chain'] else None),
                                       reserve=dict(delta_o=num(OL['res'][str(j)]['reserve']['delta_o']), used=num(OL['res'][str(j)]['reserve']['used'])),
                                       pointterm=OL['res'][str(j)]['pointterm']) for j in range(3)}
    out['lemmaN'] = 'Nbar_k^2 = 2 (sqrt(beta Bhat_k/w_min) + p^k c_dU 700^p sqrt(ln 1.5))^2 + n_III_k (STATEMENTS_r2 Lemma N, r_chi = 253, tau_III = 4.071); this file is affine in Nbar^2'
    L.dump(out, O + '/A_interface.json', list(F.values()) + [os.path.abspath(__file__)])
    t = out['res']['tight']
    print('A_interface.json: EA_ref(tight) =', [x['double'] for x in t['EA_ref']], 'R5 addition at 4 N2ref, level 2:', t['R5prime_addition']['at_4N2ref'][2]['double'],
          'rule ok:', t['rule_reproduces_EA_ref'])

if __name__ == '__main__':
    main()
