# gen_numbers_A_v37.py -- generator of fix12_texts/numbers_A_v37.tex (and numbers_A_v37.json) from the v36 chain on fin13v4b.
# Derived from fix2/n35_gen_numbers.py (v35); that script is not edited.  Same rounding rules:
# bounds up (Decimal ROUND_CEILING), lower bounds and margins down (ROUND_FLOOR), from the arb midpoint +/- radius.
# Usage: PYTHONHASHSEED=0 nice -n 19 python3 -B gen_numbers_A_v37.py <out_dir>
import json, re, sys, hashlib, os
from decimal import Decimal as Dc, getcontext, ROUND_CEILING, ROUND_FLOOR
getcontext().prec = 80
B = '/home/claude/s1spec/writeup/lemmaU/bulk/'
N = B + 'fix6_v36/'
V4B = '/home/claude/s1spec/step3nl_v4b/runs/P/'
SRC = {
    'A': N + 't2_A/out/A_interface.json',
    'A4': N + 't2_A/out/asm_design_N2ref4.json',
    'A1': N + 't2_A/out/asm_design_N2ref.json',
    'HJW': N + 't2_A/out/hjw_rep.json',
    'U': N + 't4_cert/out/U_cert_v36_pass2.json',
    'U34': N + 't4_cert/out/u34_primary.json',
    'R3': N + 'r3/r3_forcing.json',
    'Z': N + 't5_zone3/out/zone3_interface.json',
    'H': '/home/claude/s1spec/step4/p2/out_v33/h2u_v33_arb.json',
    'D': '/home/claude/s1spec/step4/p2/out_v33/b0_dprof.json',
    'RT15': V4B + 'rate_table_fin13v4b_a0.00015.json',
    'RT125': V4B + 'rate_table_fin13v4b_a0.000125.json',
    'MS': V4B + 'modsups_fin13v4b/final_alpha_points_5e-8.json',
    'FX': B + 'fix8_flux2_impl/out/flux2_cert.json',
    'FXN': B + 'fix8_flux2_impl/texts/numbers_flux2.json',
    'G1': B + 'fix10_proofs/G1_check.json',
}
def J(k): return json.load(open(SRC[k]))
def sha(p): return hashlib.sha256(open(p, 'rb').read()).hexdigest()

def val(x, up):
    """x: {'arb':...}, str '[m +/- r]', [mant, exp, rad] (exact mant*2^exp), or float"""
    if isinstance(x, dict) and 'arb' in x: x = x['arb']
    if isinstance(x, list):
        m = Dc(int(x[0])) * (Dc(2) ** int(x[1])); r = Dc(repr(float(x[2])))
    elif isinstance(x, str):
        mm = re.match(r'\[(\S+)(?: \+/- (\S+))?\]$', x)
        if mm: m = Dc(mm.group(1)); r = Dc(mm.group(2)) if mm.group(2) else Dc(0)
        else:
            mm = re.match(r'\[\+/- (\S+)\]$', x)
            if mm: m = Dc(0); r = Dc(mm.group(1))
            elif '/' in x: a, b = x.split('/'); m = Dc(a) / Dc(b); r = Dc(0)
            else: m = Dc(x); r = Dc(0)
    else: m = Dc(repr(float(x))); r = Dc(0)
    return m + r if up else m - r

def rnd(x, n, up):
    x = Dc(x); e = x.adjusted(); q = Dc(1).scaleb(e - n + 1)
    r = x.quantize(q, rounding=ROUND_CEILING if up else ROUND_FLOOR)
    if len(r.as_tuple().digits) > n: r = rnd(r, n, up)
    return r

def fmt(x, n, up):
    r = rnd(x, n, up); e = r.adjusted()
    digs = ''.join(map(str, r.as_tuple().digits)).ljust(n, '0')[:n]
    if -2 <= e <= 5: return format(r, 'f')
    return f"{digs[0]}.{digs[1:]}\\cdot10^{{{e}}}" if n > 1 else f"{digs[0]}\\cdot10^{{{e}}}"

A = J('A'); T = A['res']['tight']; A4 = J('A4')['res']['tight']; A1 = J('A1')['res']['tight']
U = J('U'); H = J('H'); D = J('D'); HJ = J('HJW'); R3 = J('R3'); U34 = J('U34')
RT = J('RT15'); RT2 = J('RT125'); MS = J('MS'); FX = J('FX'); FXN = J('FXN')
def tup(vals, n, up): return '(' + ',\\,'.join(fmt(v, n, up) for v in vals) + ')'
M = {}; RAW = {}
def put(name, tex, raw, src, key, direction):
    M[name] = tex; RAW[name] = {'tex': tex, 'raw': [str(r) for r in raw] if isinstance(raw, list) else str(raw), 'source': src, 'key': key, 'dir': direction}
# ---- datum (unchanged step-4 sources) ----
hr = [val(r['bound'], True) for r in H['rows']]; put('HtwoCertV', tup(hr, 4, True), hr, 'H', 'rows[].bound', 'up')
hm = [val(r['margin'], False) for r in H['rows']]; put('HtwoMarV', '\\times' + ',\\ \\times'.join(fmt(v, 4, False) for v in hm), hm, 'H', 'rows[].margin', 'down')
de = [val(x, True) for x in D['at_alpha_t1']['E_14_60']]; put('dEtOneV', tup(de, 4, True), de, 'D', 'at_alpha_t1.E_14_60', 'up')
for nm, k in [('dsupRampV', 'sup_ramp'), ('dsupCoreV', 'sup_core'), ('dWV', 'W_ramp')]:
    v = [val(x, True) for x in T['B0'][k]]; put(nm, tup(v, 3, True), v, 'A', 'res.tight.B0.' + k, 'up')
# ---- derived energies ----
v = [val(x, True) for x in U['supersolution']['EA']]; put('EhV', tup(v, 5, True), v, 'U', 'supersolution.EA', 'up')
v = [val(x, True) for x in T['EA0']]; put('EAzeroV', tup(v, 5, True), v, 'A', 'res.tight.EA0', 'up')
v = [val(x, False) for x in T['gp']]; put('gpV', tup(v, 5, False), v, 'A', 'res.tight.gp', 'down')
gc = A['interface28']['Ghat_core']
v = [val(gc[0][0], True), val(gc[1][1], True), val(gc[2][2], True)]; put('GcoreV', tup(v, 6, True), v, 'A', 'interface28.Ghat_core diag', 'up')
v = [val(gc[1][0], True), val(gc[2][0], True), val(gc[2][1], True)]
put('GcoreOffV', '\\hat G_{10}\\le%s,\\ \\hat G_{20}\\le%s,\\ \\hat G_{21}\\le%s' % tuple(fmt(x, 5, True) for x in v), v, 'A', 'interface28.Ghat_core off-diag', 'up')
v = val(HJ['kappa_J'], True); put('kappaJV', fmt(v, 5, True), v, 'HJW', 'kappa_J (certified + pad 1e-9)', 'up')
v = [val(x, True) for x in T['Dsrc'][:2]]; put('DsrcV', tup(v, 5, True), v, 'A', 'res.tight.Dsrc[0:2]', 'up')
v = val(T['rhoJ'], True); put('rhoJV', fmt(v, 5, True), v, 'A', 'res.tight.rhoJ', 'up')
put('rexpV', '(\\tfrac{1249}{1000},\\,\\tfrac{1249}{1000},\\,\\tfrac1{1000})', T['r_export'], 'A', 'res.tight.r_export', 'exact')
v = val(T['amin_c'], False); put('acV', fmt(v, 5, False), v, 'A', 'res.tight.amin_c', 'down')
v = [val(x, True) for x in T['F2'][:3]]; put('FtwoV', tup(v, 5, True), v, 'A', 'res.tight.F2[0:3]', 'up')
v = [val(T['coup'][1], True), val(T['coup'][2], True)]; put('coupV', '(0,\\,%s,\\,%s)' % (fmt(v[0], 5, True), fmt(v[1], 5, True)), v, 'A', 'res.tight.coup[1:3]', 'up')
v = [val(x, True) for x in T['C1_0'][:2]]; put('ConeZeroV', tup(v, 5, True), v, 'A', 'res.tight.C1_0[0:2]', 'up')
v = [val(x, True) for x in U['supersolution']['C0']]; put('CzeroV', tup(v, 5, True), v, 'U', 'supersolution.C0', 'up')
v = [val(x, True) for x in U['supersolution']['C1']]; put('ConeV', tup(v, 5, True), v, 'U', 'supersolution.C1', 'up')
v = [val(x, True) for x in U['supersolution']['N2']]; put('NbarV', tup(v, 5, True), v, 'U', 'supersolution.N2', 'up')
v = val(T['kappa_o'][2][2], True); put('kapOtwotwoV', fmt(v, 5, True), v, 'A', 'res.tight.kappa_o[2][2]', 'up')
# ---- level-2 export (r_2 = 1/1000) ----
r2 = T['r2_export']
for nm, k, n in [('CzeroTwoRefV', 'C0_ref', 6), ('ConeTwoRefV', 'C1_ref', 6), ('CzeroTwoRefFourV', 'C0_ref4', 7), ('ConeTwoRefFourV', 'C1_ref4', 7), ('ConeTwoZeroV', 'C1_0', 6)]:
    v = val(r2[k], True); put(nm, fmt(v, n, True), v, 'A', 'res.tight.r2_export.' + k, 'up')
wd = r2['window_data']
for nm, k in [('EAtwoV', 'EA_2'), ('DtwoV', 'D_2'), ('DtwoSlopeV', 'D_window_slope')]:
    v = val(wd[k], True); put(nm, fmt(v, 6, True), v, 'A', 'res.tight.r2_export.window_data.' + k, 'up')
put('BcCoreTwoV', '\\tfrac{1688999}{10^6}', '1688999/1000000', 'G1', 'B1 (weights_t2 core B_c at level 2)', 'exact')
# ---- R3: honest forcing constant ----
v = val(R3['r0_new'], True); put('rzeroV', fmt(v, 6, True), v, 'R3', 'r0_new', 'up')
v = val(R3['DL'], True); put('DLV', '3.2224577463\\cdot10^{-6}', v, 'R3', 'DL (exact)', 'exact')
for j in range(3):
    v = val(R3['levels'][str(j)]['Feps2'] if isinstance(R3['levels'], dict) else R3['levels'][j]['Feps2'], True)
    put('FepsV%s' % 'abc'[j], fmt(v, 4, True), v, 'R3', 'levels.%d.Feps2' % j, 'up')
v = val(R3['c_dU']['new'], True); put('cdUV', fmt(v, 6, True), v, 'R3', 'c_dU.new', 'up')
v = val(R3['c_dU']['betabar_prime'], True); put('betabarPrimeV', fmt(v, 6, True), v, 'R3', 'c_dU.betabar_prime', 'up')
v = Dc('2e-8') / val(R3['c_dU']['new'], True); put('cdUMarV', fmt(v, 4, False), v, 'R3', '2e-8 / c_dU.new', 'down')
# ---- core / ramp conclusions and margins at 4 N2_ref ----
v = val(A4['S0'], True); put('SzerovV', fmt(v, 5, True), v, 'A4', 'res.tight.S0', 'up')
v = val(A4['S1'], True); put('SonevV', fmt(v, 5, True), v, 'A4', 'res.tight.S1', 'up')
fl4 = A4['flux']['beta']['2']
v = val(fl4['margin_FLUX0'], False); put('MarFluxZeroV', fmt(v, 4, False), v, 'A4', 'flux.beta.2.margin_FLUX0 (RECORD, old flux block)', 'down')
v = val(fl4['margin_FLUX2'], False); put('MarFluxTwoV', fmt(v, 4, False), v, 'A4', 'flux.beta.2.margin_FLUX2 (RECORD, old flux block)', 'down')
mg = A4['margins']
for nm, k in [('MarSzeroV', 'S0v_vs_0.0153'), ('MarSoneV', 'S1v_vs_0.6'), ('MarMoneV', 'M1core'), ('MarMzeroV', 'M0core')]:
    v = val(mg[k], False); put(nm, fmt(v, 4, False), v, 'A4', 'res.tight.margins.' + k, 'down')
cm = T['core']['margins']
for nm, k in [('MarGV', 'Ghat22'), ('MarRampMoneV', 'ramp M1'), ('MarRampZoneV', 'ramp z1'), ('MarRampMzeroV', 'ramp M0')]:
    v = val(cm[k], False); put(nm, fmt(v, 4, False), v, 'A', 'res.tight.core.margins.' + k, 'down')
fr = T['core']['flux_margin']['2']
v = [val(fr['FLUX0_over_K0'], False), val(fr['FLUX2_over_K2'], False)]
put('MarFluxRefV', '%s,\\ %s' % (fmt(v[0], 4, False), fmt(v[1], 4, False)), v, 'A', 'res.tight.core.flux_margin.2 (RECORD, old flux block)', 'down')
v = val(A['pointterm']['gc_28'], True); put('gcV', fmt(v, 5, True), v, 'A', 'pointterm.gc_28', 'up')
v = max(val(T_i['Phi'], True) for T_i in A4['info']); put('PhiV', fmt(v, 5, True), v, 'A4', 'max_j res.tight.info[j].Phi', 'up')
# ---- verifier fin13v4b values used by (A) ----
fm = RT['far_max']
v = Dc(repr(fm['M1_ramp'])); put('MoneRampV', fmt(v, 7, True), v, 'RT15', 'far_max.M1_ramp', 'up')
v = Dc(repr(fm['sI1'])); put('sIoneV', fmt(v, 7, True), v, 'RT15', 'far_max.sI1', 'up')
v = Dc(repr(fm['sI2'])); put('sItwoV', fmt(v, 5, True), v, 'RT15', 'far_max.sI2', 'up')
v = Dc(repr(RT['S_ramp_max_along_trajectory'])); put('SrampV', fmt(v, 4, True), v, 'RT15', 'S_ramp_max_along_trajectory', 'up')
v = Dc(repr(RT['sqrtL_max'])); put('peakV', fmt(v, 6, True), v, 'RT15', 'sqrtL_max', 'up')
v = Dc(repr(RT['tail']['rho'])); put('rhoTailV', fmt(v, 5, False), v, 'RT15', 'tail.rho', 'down')
v = Dc(repr(RT['tail']['sqrt_LF'])); put('sqrtLFV', fmt(v, 5, True), v, 'RT15', 'tail.sqrt_LF', 'up')
v = Dc(repr(RT['tail']['alpha_F'])); put('alphaFV', fmt(v, 5, True), v, 'RT15', 'tail.alpha_F', 'up')
v = Dc(repr(RT['tail']['alpha_condition_whole_tail']['margin_whole_tail'])); put('tailAlphaMarV', fmt(v, 5, False), v, 'RT15', 'tail.alpha_condition_whole_tail.margin_whole_tail', 'down')
v = Dc(repr(RT2['tail']['alpha_condition_whole_tail']['margin_whole_tail'])); put('tailAlphaMarBV', fmt(v, 5, False), v, 'RT125', 'tail.alpha_condition_whole_tail.margin_whole_tail', 'down')
v = Dc(str(MS['1.5e-4']['sup_ab'])); put('abarV', fmt(v, 6, True), v, 'MS', '1.5e-4.sup_ab', 'up')
v = Dc(str(MS['1.5e-4']['sup_beta'])); put('betabarV', fmt(v, 6, True), v, 'MS', '1.5e-4.sup_beta', 'up')
for m, nm in [('V0sq', 'JzeroV'), ('V1sq', 'JoneV'), ('V2sq', 'JtwoV'), ('V4sq', 'JfourV')]:
    c = Dc(repr(RT['unc'][m]['const'])); s = Dc(repr(RT['unc'][m]['slope_per_tau']))
    put(nm, '%s,\\ %s' % (fmt(c, 5, True), fmt(s, 5, True)), [c, s], 'RT15', 'unc.%s.{const,slope_per_tau}' % m, 'up')
# ---- Lemma U, v36 ----
v = [val(x, True) for x in U['supersolution']['Bhat_num']]; put('BhV', tup(v, 5, True), v, 'U', 'supersolution.Bhat_num', 'up')
um = U['margins']
v = [val(um['late_M1'], False), val(um['late_M0'], False), Dc(repr(um['present_M1_min']['value'])), Dc(repr(um['present_M0_min']['value']))]
put('UmarV', 'late M1 $\\times%s$, late M0 $\\times%s$, present M1 $\\ge\\times%s$, present M0 $\\ge\\times%s$' % tuple(fmt(x, 4, False) for x in v), v, 'U', 'margins', 'down')
v = [Dc(repr(x)) for x in U['ufine']['Bhat']]; put('BhFineV', tup(v, 5, True), v, 'U', 'ufine.Bhat', 'up')
v = val(U['supersolution']['EA'][0], False) / Dc('1.01e-6'); put('KsixV', fmt(v, 4, False), v, 'U', 'supersolution.EA[0] / 1.01e-6', 'down')
# ---- flux lemma v2 (fix8) ----
res = FX['res']
for nm, k, n in [('KzeroBarV', 'K0bar', 5), ('KtwoPtV', 'K2pt', 5), ('CprimeV', 'Cp', 5), ('CppV', 'Cpp', 5)]:
    v = val(res[k]['arb'], True); put(nm, fmt(v, n, True), v, 'FX', 'res.%s' % k, 'up')
# combined margin (FLUX0 + eps FLUX2)/(K0bar + eps K2pt), lower bound
eps = Dc('3e-5'); comb = (Dc('1.64e-20') + eps * Dc('8.0e-20')) / (val(res['K0bar']['arb'], True) + eps * val(res['K2pt']['arb'], True))
put('FluxCombMarV', fmt(comb, 4, False), comb, 'FX', '(FLUX0+eps_I FLUX2)/(K0bar+eps_I K2pt)', 'down')
if __name__ == '__main__':
    out = sys.argv[1] if len(sys.argv) > 1 else '.'
    hdr = ['% numbers_A_v37.tex -- THE NUMERICAL VALUES of Lemma S(A), Lemma U and the verifier inputs, v37 (chain v36 on verifier fin13v4b).',
           '% GENERATED by fix12_texts/tools/gen_numbers_A_v37.py; do not edit by hand.  Rounding: bounds up, lower bounds and margins down.',
           '% Supersedes numbers_A_v33.tex (v35 values; that file is unchanged).  Macro names of v33 are kept; new macros are added.',
           '% Sources (sha256 at generation time):']
    for k, p in SRC.items(): hdr.append('%%   %-5s = %s  (%s)' % (k, p, sha(p)[:16]))
    hdr.append('% RECORD ONLY: MarFluxZeroV, MarFluxTwoV, MarFluxRefV come from the (A)-chain flux block (old kappa, old W(57));')
    hdr.append('% they are superseded by flux lemma v2 (FX: KzeroBarV, KtwoPtV, FluxCombMarV) and are not used by any v37 argument.')
    lines = hdr + ['\\providecommand{\\%s}{}\\renewcommand{\\%s}{%s}' % (k, k, v) for k, v in M.items()]
    open(os.path.join(out, 'numbers_A_v37.tex'), 'w').write('\n'.join(lines) + '\n')
    js = {'args': {'argv': sys.argv, 'env': {k: os.environ.get(k) for k in ('PYTHONHASHSEED',)}, 'inputs': {k: {'path': p, 'sha256': sha(p)} for k, p in SRC.items()},
                   'script_sha256': sha(os.path.abspath(__file__))}, 'macros': RAW}
    json.dump(js, open(os.path.join(out, 'numbers_A_v37.json'), 'w'), indent=1, sort_keys=True)
    for k, v in M.items(): print(k, '=>', v)
