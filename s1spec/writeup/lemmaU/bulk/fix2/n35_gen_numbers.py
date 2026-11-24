# N35: generator of the macro values of fix1/b12_text/numbers_A_v33.tex from the v35 sources (python3 -B n35_gen_numbers.py prints name => value).
# Two tweaks were applied when patching the file: dEtOneV at 4 significant digits, PhiV = max over the three info[j].Phi entries (5 digits, up).
# Outward rounding: bounds up (Decimal ROUND_CEILING), margins down (ROUND_FLOOR).
# regenerate the numeric macros of numbers_A_v33.tex from the v35 sources (round UP for bounds, DOWN for margins)
import json, re, sys
from decimal import Decimal as Dc, getcontext, ROUND_CEILING, ROUND_FLOOR
getcontext().prec = 60
B='/home/claude/s1spec/writeup/lemmaU/bulk/'
def J(p): return json.load(open(p))
def val(x, up):
    """x: arb-dict {'arb': str}, str '[m +/- r]', {'arb':[mant,exp,rad]}, or float"""
    if isinstance(x, dict) and 'arb' in x: x = x['arb']
    if isinstance(x, list):
        m = Dc(int(x[0])) * (Dc(2) ** int(x[1])); r = Dc(repr(x[2])); 
    elif isinstance(x, str):
        mm = re.match(r'\[(\S+)(?: \+/- (\S+))?\]$', x)
        if mm: m = Dc(mm.group(1)); r = Dc(mm.group(2)) if mm.group(2) else Dc(0)
        else: m = Dc(x); r = Dc(0)
    else: m = Dc(repr(float(x))); r = Dc(0)
    return m + r if up else m - r
def rnd(x, n, up):
    x = Dc(x); e = x.adjusted(); q = Dc(1).scaleb(e - n + 1)
    return x.quantize(q, rounding=ROUND_CEILING if up else ROUND_FLOOR)
def fmt(x, n, up):
    r = rnd(x, n, up); s = format(r.normalize() if False else r, 'f') if False else None
    e = r.adjusted(); digs = str(abs(r.as_tuple().digits and int(''.join(map(str, r.as_tuple().digits))))).rjust(n, '0')[:n]
    # r has exactly n significant digits after quantize unless rounding carried (e.g. 9.99->10.0): re-quantize
    if len(r.as_tuple().digits) > n: r = rnd(r, n, up); e = r.adjusted(); 
    digs = ''.join(map(str, r.as_tuple().digits)).ljust(n, '0')[:n]
    if e >= -2 and e <= 5:
        return format(r, 'f')
    return f"{digs[0]}.{digs[1:]}\\cdot10^{{{e}}}" if n > 1 else f"{digs[0]}\\cdot10^{{{e}}}"
A=J(B+'fix2/n12/t2_A_v35/out/A_interface.json'); T=A['res']['tight']
U=J(B+'t4_cert/out/U_cert_v35.json'); A4=J(B+'fix2/n12/t2_A_v35/out/asm_design_N2ref4.json')['res']['tight']
H=J('/home/claude/s1spec/step4/p2/out_v33/h2u_v33_arb.json'); D=J('/home/claude/s1spec/step4/p2/out_v33/b0_dprof.json')
def tup(vals, n, up, wrap='({})'): return wrap.format(',\\,'.join(fmt(v, n, up) for v in vals))
M = {}
M['HtwoCertV'] = tup([val(r['bound'], True) for r in H['rows']], 4, True)
M['HtwoMarV'] = '\\times' + ',\\ \\times'.join(fmt(val(r['margin'], False), 4, False) for r in H['rows'])
M['dEtOneV'] = tup([val(x, True) for x in D['at_alpha_t1']['E_14_60']], 3, True)
M['dsupRampV'] = tup([val(x, True) for x in T['B0']['sup_ramp']], 3, True)
M['dsupCoreV'] = tup([val(x, True) for x in T['B0']['sup_core']], 3, True)
M['dWV'] = tup([val(x, True) for x in T['B0']['W_ramp']], 3, True)
M['EhV'] = tup([val(x, True) for x in U['supersolution']['EA']], 5, True)
M['EAzeroV'] = tup([val(x, True) for x in T['EA0']], 5, True)
M['gpV'] = tup([val(x, False) for x in T['gp']], 5, False)
gc = A['interface28']['Ghat_core']
M['GcoreV'] = tup([val(gc[0][0], True), val(gc[1][1], True), val(gc[2][2], True)], 6, True)
M['GcoreOffV'] = '\\hat G_{10}\\le%s,\\ \\hat G_{20}\\le%s,\\ \\hat G_{21}\\le%s' % (fmt(val(gc[1][0], True), 5, True), fmt(val(gc[2][0], True), 5, True), fmt(val(gc[2][1], True), 5, True))
M['kappaJV'] = fmt(val(A['HJW']['kappa_J'], True), 5, True)
M['DsrcV'] = tup([val(x, True) for x in T['Dsrc'][:2]], 5, True)
M['rhoJV'] = fmt(val(T['rhoJ'], True), 5, True)
M['acV'] = fmt(val(T['amin_c'], False), 5, False)
M['FtwoV'] = tup([val(x, True) for x in T['F2'][:2]], 5, True)
M['coupV'] = '(0,\\,%s)' % fmt(val(T['coup'][1], True), 5, True)
M['ConeZeroV'] = tup([val(x, True) for x in T['C1_0'][:2]], 5, True)
M['CzeroV'] = tup([val(x, True) for x in U['supersolution']['C0']], 5, True)
M['ConeV'] = tup([val(x, True) for x in U['supersolution']['C1']], 5, True)
M['NbarV'] = tup([val(x, True) for x in U['supersolution']['N2']], 5, True)
M['kapOtwotwoV'] = fmt(val(T['kappa_o'][2][2], True), 5, True)
M['SzerovV'] = fmt(val(A4['S0'], True), 5, True)
M['SonevV'] = fmt(val(A4['S1'], True), 5, True)
fl4 = A4['flux']['beta']['2'] if 'beta' in A4['flux'] and '2' in A4['flux']['beta'] else None
M['_fl4'] = json.dumps(list(A4['flux']['beta']['2'].keys()) if fl4 else None)
for k, v in [('MarFluxZeroV', fl4 and fl4['margin_FLUX0']), ('MarFluxTwoV', fl4 and fl4['margin_FLUX2'])]: M[k] = fmt(val(v, False), 4, False)
mg = A4['margins']
M['MarSzeroV'] = fmt(val(mg['S0v_vs_0.0153'], False), 4, False); M['MarSoneV'] = fmt(val(mg['S1v_vs_0.6'], False), 4, False)
M['MarMoneV'] = fmt(val(mg['M1core'], False), 4, False); M['MarMzeroV'] = fmt(val(mg['M0core'], False), 4, False)
cm = T['core']['margins']
M['MarGV'] = fmt(val(cm['Ghat22'], False), 4, False)
M['MarRampMoneV'] = fmt(val(cm['ramp M1'], False), 4, False); M['MarRampZoneV'] = fmt(val(cm['ramp z1'], False), 4, False); M['MarRampMzeroV'] = fmt(val(cm['ramp M0'], False), 4, False)
fr = T['core']['flux_margin']['2']
M['MarFluxRefV'] = '%s,\\ %s' % (fmt(val(fr['FLUX0_over_K0'], False), 4, False), fmt(val(fr['FLUX2_over_K2'], False), 4, False))
M['gcV'] = fmt(val(T['pointterm']['gc_28'], True), 5, True) if 'pointterm' in T else fmt(val(A['pointterm']['gc_28'], True), 5, True)
M['PhiV'] = fmt(val(A4['info'][2]['Phi'], True), 5, True) if False else None
M['BhV'] = tup([val(x, True) for x in U['supersolution']['Bhat_num']], 5, True)
um = U['margins']
M['UmarV'] = 'late M1 $\\times%s$, late M0 $\\times%s$, present M1 $\\ge\\times%s$, present M0 $\\ge\\times%s$' % (fmt(val(um['late_M1'], False), 4, False), fmt(val(um['late_M0'], False), 4, False), fmt(Dc(repr(um['present_M1_min']['value'])), 4, False), fmt(Dc(repr(um['present_M0_min']['value'])), 4, False))
M['BhFineV'] = tup([Dc(repr(x)) for x in U['ufine']['Bhat']], 5, True)
M['KsixV'] = fmt(val(U['supersolution']['EA'][0], False) / Dc('1.01e-6'), 4, False)
if __name__ == '__main__':
    for k, v in M.items(): print(k, '=>', v)
