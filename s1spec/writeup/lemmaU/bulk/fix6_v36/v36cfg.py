"""fix6_v36 (tag v36): single source of the verifier inputs and of the two switches of the (A) -> Lemma U -> zone III chain.

Environment (recorded in every args block through arlib.args_block / t4lib.args_block):
  V36_VER = fin13v4b (default, production) | fin12f4 (regression against v35 / fix4_r3)
  V36_R3  = 1 (default) | 0      R3 correction of fix4_r3 (honest forcing r0', eps-part, c_dU with DL)
  V36_R2  = 1 (default) | 0      level-2 dissipation export r_2 = 1/1000 of (A) (ROUTE_DECISION.md unit G1 / F2')
Production = (fin13v4b, 1, 1).  Regression 1 = (fin12f4, 0, 0) must reproduce v35; regression 2 = (fin12f4, 1, 0) must reproduce fix4_r3.

Verifier constants handed to the base_ibp reference chain (lemmaS_F2 argv) are read from the raw verifier JSON and rounded OUTWARD here:
  J_m, slopes, J4, sJ4, sqrt(Lambda) peak: up to 5 significant digits (Decimal ROUND_CEILING on the exact double);
  tail rho (a decay rate, used as a lower bound): down to 4 significant digits (ROUND_FLOOR);
  abar = sup|(a,b)|, betabar = sup|beta|: the mod_sups strings (6 digits, already rounded up by mod_sups.py) as written.
With V36_VER = fin12f4 this reproduces the strings typed in v35 (run_F10p_ibp.py, capture_F2.py) exactly (asserted)."""
import os, json, decimal

N = '/home/claude/s1spec/writeup/lemmaU/bulk/fix6_v36'
VER = os.environ.get('V36_VER', 'fin13v4b')
R3ON = os.environ.get('V36_R3', '1') == '1'
R2ON = os.environ.get('V36_R2', '1') == '1'
assert VER in ('fin12f4', 'fin13v4b'), VER
assert os.environ.get('V36_R3', '1') in ('0', '1') and os.environ.get('V36_R2', '1') in ('0', '1')
ENVKEYS = ('V36_VER', 'V36_R3', 'V36_R2')

if VER == 'fin12f4':
    VDIR = '/home/claude/s1spec/step3nl_v33/verifier'
    RT = VDIR + '/rate_table_fin12f4_a0.00015.json'
    RT125 = VDIR + '/rate_table_fin12f4_a0.000125.json'
    MODS = VDIR + '/modsups_fin12f4/final_alpha_points_5e-8.json'
    CONSTPY = VDIR + '/constants.py'
    N1CHART = '/home/claude/s1spec/step4/p2/out_v33/n1_chart.json'          # the pinned step-4 v33 stage 19 output (v35 input)
else:
    VDIR = '/home/claude/s1spec/step3nl_v4b/verifier'
    RT = '/home/claude/s1spec/step3nl_v4b/runs/P/rate_table_fin13v4b_a0.00015.json'
    RT125 = '/home/claude/s1spec/step3nl_v4b/runs/P/rate_table_fin13v4b_a0.000125.json'
    MODS = '/home/claude/s1spec/step3nl_v4b/runs/P/modsups_fin13v4b/final_alpha_points_5e-8.json'
    CONSTPY = VDIR + '/constants.py'
    N1CHART = N + '/n1/out/n1_chart_fin13v4b.json'                          # n1_chart.py (copy) rerun on the fin13v4b sups
DPROF = '/home/claude/s1spec/step4/p2/out_v33/b0_dprof.json'                 # verifier-independent (no verifier input; checked in V36_LOG.md)
R2_EXPORT = '1/1000' if R2ON else None
# kappa_J reserve (flux-lemma v2 request, fix7_flux2/FLUX2_DESIGN.md C6): kappa_J used by (A) = certified kappa_J + KJ_PAD (absolute, i.e. J_m gains
# a relative 1e-9 on top of the certified (1 + kappa_J)); production only, so the fin12f4 regressions reproduce v35 / fix4_r3.
KJ_PAD = '1/1000000000' if VER == 'fin13v4b' else '0'

# the FLUX0/FLUX2/S0V/S1V inputs of the verifier run must be the ones (A) is checked against (thresholds typed in assembly_rep.py)
FLUX0, FLUX2, S0V, S1V = '1.64e-20', '8.0e-20', '0.0153', '0.6'


def env_record():
    return {k: os.environ.get(k) for k in ENVKEYS}


def _dec(x):
    return decimal.Decimal(repr(float(x)))


def round_sig(x, sig, up=True):
    """decimal string with sig significant digits, rounded up (ceiling) or down (floor); asserted outward"""
    d = _dec(x); q = decimal.Decimal(1).scaleb(d.adjusted() - (sig - 1))
    r = d.quantize(q, rounding=decimal.ROUND_CEILING if up else decimal.ROUND_FLOOR)
    assert (r >= d) if up else (r <= d)
    e = r.adjusted()
    mant = r.scaleb(-e).normalize()
    s = '%se%d' % (format(mant, 'f'), e)
    assert decimal.Decimal(s) == r
    return s


def verifier_inputs():
    """raw verifier JSON -> outward-rounded strings for the base_ibp argv (alpha_1 = 1.5e-4 row)"""
    rt = json.load(open(RT)); ms = json.load(open(MODS))['1.5e-4']; u = rt['unc']
    assert rt['alpha1'] == 0.00015
    a = rt.get('args')
    if a:   # fin13v4b stores its environment: check the inputs (A) is compared against
        env = a.get('env', {})
        for k, v in (('FLUX0', FLUX0), ('FLUX2', FLUX2), ('S0V', S0V), ('S1V', S1V)):
            if k in env: assert env[k] == v, (k, env[k])
    J = [round_sig(u[k]['const'], 5) for k in ('V0sq', 'V1sq', 'V2sq')]
    SL = [round_sig(u[k]['slope_per_tau'], 5) for k in ('V0sq', 'V1sq', 'V2sq')]
    J4 = round_sig(u['V4sq']['const'], 5); SJ4 = round_sig(u['V4sq']['slope_per_tau'], 5)
    SQ = round_sig(rt['sqrtL_max'], 5); RHO = round_sig(rt['tail']['rho'], 4, up=False)
    assert ms['peak'] == rt['sqrtL_max'] and ms['rho'] == rt['tail']['rho']
    for k in ('V0sq', 'V1sq', 'V2sq', 'V4sq'): assert ms['J'][k] == [u[k]['const'], u[k]['slope_per_tau']], k
    return dict(J=J, slopes=SL, J4=J4, sJ4=SJ4, sqLam=SQ, rho=RHO, abar=sci(ms['sup_ab']), beta=sci(ms['sup_beta']))


def sci(x):
    """exact decimal of the double's shortest repr, written as <mantissa>e<exp> (no rounding)"""
    d = _dec(x); e = d.adjusted(); s = '%se%d' % (format(d.scaleb(-e).normalize(), 'f'), e)
    assert decimal.Decimal(s) == d
    return s


def f2_argv(tag):
    v = verifier_inputs()
    return ['lemmaS_F2.py', '--J', *v['J'], '--handoff', '1e-8', '1e-6', '1e-3', '--sqLam', v['sqLam'], '--alpha1', '1.5e-4',
            '--core', '0.11', '0.52', '--ramp', '0.05', '0.9', '--mhr', '0', '--slopes', *v['slopes'], '--rho', v['rho'], '--M1s', '0.52',
            '--J4', v['J4'], '--sJ4', v['sJ4'], '--rampl4', '--rs33', 'ibp2', '--s33', 'ibp', '--l4', '--abar', v['abar'], '--beta', v['beta'], '--tag', tag]


V35_TYPED = ['--J', '6.0273e-12', '8.4673e-12', '1.8708e-12', '--handoff', '1e-8', '1e-6', '1e-3', '--sqLam', '5.3075e-8', '--alpha1', '1.5e-4',
             '--core', '0.11', '0.52', '--ramp', '0.05', '0.9', '--mhr', '0', '--slopes', '3.4548e-18', '4.8534e-18', '1.0724e-18', '--rho', '2.693e-3', '--M1s', '0.52',
             '--J4', '6.2362e-8', '--sJ4', '3.5746e-14', '--rampl4', '--rs33', 'ibp2', '--s33', 'ibp', '--l4', '--abar', '2.19103e-7', '--beta', '3.33811e-9']
if VER == 'fin12f4':
    assert f2_argv('x')[1:-2] == V35_TYPED, (f2_argv('x'), V35_TYPED)


if __name__ == '__main__':
    print(json.dumps(dict(VER=VER, R3ON=R3ON, R2ON=R2ON, RT=RT, MODS=MODS, N1CHART=N1CHART, inputs=verifier_inputs()), indent=1))
