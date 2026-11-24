# v39_values.py -- every record value that the v39 texts print, read from tools/records_v39.json (written by gen_records_v39.py from
# disk) and rounded outward: upper bounds up, lower bounds and margins down.  No typed literal of a record value appears in the builder.
import json
from decimal import Decimal as Dc, ROUND_CEILING, ROUND_FLOOR
D = '/home/claude/s1spec/writeup/lemmaU/bulk/fix17_texts_v39/'
R = json.load(open(D + 'tools/records_v39.json'))
I = R['args']['inputs']; res = R['results']
def _r(x, k, mode):
    v = Dc(repr(float(x))) if not isinstance(x, Dc) else x
    e = v.adjusted(); q = Dc(1).scaleb(e - k + 1)
    return v.quantize(q, rounding=mode)
def tex(v):
    """Decimal -> TeX string, scientific with \\cdot10^{e} below 1e-3."""
    if v != 0 and abs(v) < Dc('0.001'):
        e = v.adjusted(); m = v.scaleb(-e)
        s = format(m.normalize(), 'f') if m != m.to_integral() else str(int(m))
        return '%s\\cdot10^{%d}' % (s, e)
    s = format(v, 'f')
    return s
def up(x, k): return tex(_r(x, k, ROUND_CEILING))
def dn(x, k): return tex(_r(x, k, ROUND_FLOOR))
def sh(name): return I[name]['sha256']
def sh8(name): return I[name]['sha256'][:8] + '...'
ifc = res['interface_proofs']
V = {
 'man': sh('manifest_v38'), 'man8': sh8('manifest_v38'), 'drv': sh('driver_v38'), 'rep': sh('repin_report'),
 'man33': sh('manifest_v33'), 'man32': sh('manifest_v32'), 'excB8': sh('exc_B8'), 'excV4c': sh('exc_v4c'),
 'kit3': sh('kit3'), 'kit3b': I['kit3']['bytes'], 'kit4b': sh('kit4b'), 'kit4bb': I['kit4b']['bytes'],
 'cn34': sh('cn34'), 'cn34pre': sh('cn34_pre'), 's4p': sh('step4_proof'),
 'n1v38': sh('n1_v38'), 'n1chain': sh('n1_chain'),
 'nst': len(R['stages']), 'paths': R['pins']['hashes_verified_at_end']['paths'], 'scripts': R['pins']['scripts'], 'entries': R['pins']['entries'],
 'carried': R['pins']['hashes_verified_at_end']['carried_stages'], 'newst': R['pins']['hashes_verified_at_end']['rerun_or_new_stages'],
 'v32st': R['checks']['v32_manifest_at_end']['stages'], 'v32paths': R['checks']['v32_manifest_at_end']['paths'],
 'v33st': R['checks']['v33_manifest_at_end']['stages'], 'v33paths': R['checks']['v33_manifest_at_end']['paths'],
 # interface certificates (stages 22-27)
 'ke': up(ifc['C1_ke']['sup_ke'], 5), 'keMar': dn(Dc('0.5') / Dc(repr(ifc['C1_ke']['sup_ke'])), 3),
 'gsKu': up(ifc['GS_tau1']['sup_ku'], 5), 'gsUp': up(ifc['GS_tau1']['sup_uprime'], 5),
 'gsKuMar': dn(Dc('0.5') / Dc(repr(ifc['GS_tau1']['sup_ku'])), 3), 'gsUpMar': dn(Dc(1) / Dc(repr(ifc['GS_tau1']['sup_uprime'])), 3),
 'gsXS': repr(ifc['GS_tau1']['XS']),
 'gbKu': up(ifc['GS_boot_v4c']['sup_ku'], 5), 'gbUp': up(ifc['GS_boot_v4c']['sup_uprime'], 5),
 # stage 8 (H7), 9 (notes), 18 (far H4), 20 (f12), 21 (flux)
 'abar': up(res['H7_fixedL']['abar'], 6), 'ratio': up(res['H7_fixedL']['ratio_upper'], 8), 'need': up(res['H7_fixedL']['needed'], 5),
 'tstar': dn(res['H7_fixedL']['t_star_lower'], 6), 'h7mar': dn(res['H7_fixedL']['margin_tstar_over_ratio'], 4),
 'h6mar': dn(res['notes']['H6_task5_min_margin'], 4), 'h4far': dn(res['M9_H4_far']['min_margin'], 4),
 'f12n': res['F12_sliver']['n_checks'], 'f12p': res['F12_sliver']['n_pass'], 'sI2': up(res['F12_sliver']['sI2_prime'], 6),
 'M1p': up(res['F12_sliver']['M1_prime'], 6), 'M1A': res['F12_sliver']['M1_read_by_A'],
 'flux': dn(res['flux_v37']['combined_margin_lower'], 5),
 'ifcn': ifc['n_checks'], 'pinsn': res['v4c_runs_pins']['n_checks'],
 # carried (H2) for f - U_1 and (H2f), arb rows
 'h2u': [up(r['bound'], 4) for r in res['M9_H2_tightened_fU1']['rows']], 'h2uMar': [dn(r['margin'], 4) for r in res['M9_H2_tightened_fU1']['rows']],
 'h2f': [up(r['with_error'], 4) for r in res['M9_H2f']['rows']], 'h2fMar': [dn(r['margin'], 4) for r in res['M9_H2f']['rows']],
 'kit3sum': R['kit3_summary'], 'kit4bsum': R['kit4b_summary'],
}
G1 = json.load(open(I['G1_v38']['path']))
