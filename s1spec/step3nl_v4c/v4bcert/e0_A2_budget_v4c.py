"""fin13v4c certification stage (REVIEW_v4b_v36_proofs.md A-1): the E0 A^2 absorption budget of E0_LEMMA item (iv) with the two pieces
that e0_lemma_consts.py section (10) (fin13v4b) omits:
  (A-1)  lead mismatch of (a,b): the modulation lemma's (a,b) lead (3x3 H form) differs from the 2x2 lead M^-1(<w,L*f1>,<w,L*f2>) used in Qfull
         (bloop_forms.py) by at most (DAB + dlead) A; it couples through alpha ||(gA,gB)|| and through the loop (LOOPC):
             alpha (DAB + CAB0_dlead) (GABN + LOOPC) A^2                               -> per alpha A^2: (DAB + CAB0_dlead)(GABN + LOOPC)
  (A-1b) third-order loop piece alpha^2 |r| * MAB alpha A (n3_123 H33):                 -> per alpha A^2: H33 alpha |r| MAB n3_123 alpha
Criterion (as in fin13v4b):  Pi3_T3 + (need(alpha) - 25) alpha <= PI3 = 0.01607 for |alpha| <= AMAX_E0 = 2.5e-4, where need is nondecreasing in alpha,
so the check at alpha = AMAX_E0 suffices (1.6e-4 for reference).  The verifier constants do not change (PI3, MODA = 25 kept).
Inputs: out/e0_lemma_consts.json (need parts, Pi3 T3, H33, |r| pieces), out/e0v4b_constants.json (N3P, GABN, LOOPC, C_DR_CORE),
out/e2v4b_constants.json (DAB), v4cert/out/constants_v4.json (MAB, CAB0_dlead), agent_mod/bloop/bloop.json (r1, r2).
Every float read from JSON is widened by one ulp in the safe direction (upper quantities up, the slack down). [cert] arb, 200 bit.
Run (chain, after e2_v4b_constants): PYTHONHASHSEED=0 nice -n 19 python3 -B e0_A2_budget_v4c.py -> out/e0_A2_budget_v4c.json"""
import sys, os, json, math, time
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import v4bcommon as VC
from flint import arb, ctx
ctx.prec = 200
t0 = time.time()
P = dict(lc=os.path.join(VC.OUT, 'e0_lemma_consts.json'), e0=os.path.join(VC.OUT, 'e0v4b_constants.json'), e2=os.path.join(VC.OUT, 'e2v4b_constants.json'),
         c4=os.path.join(VC.V4B, 'v4cert', 'out', 'constants_v4.json'), bloop=os.path.join(VC.V4B, 'agent_mod', 'bloop', 'bloop.json'))
J = {k: json.load(open(p)) for k, p in P.items()}
U = lambda x: arb(math.nextafter(float(x), math.inf))          # float upper quantity -> arb, one ulp up
D = lambda x: arb(math.nextafter(float(x), -math.inf))         # float lower quantity -> arb, one ulp down
up = lambda z: float(z.upper()); lo = lambda z: float(z.lower())
lc = J['lc']; E0 = J['e0']['constants']; E2 = J['e2']['constants']; C4 = J['c4']['constants']
N3P = U(E0['N3P']['value']); GABN = U(E0['GABN']['value']); LOOPC = U(E0['LOOPC']['value']); CDR = U(E0['C_DR_CORE']['value'])
DAB = U(E2['DAB']['value']); DLEAD = U(C4['CAB0_dlead']['value']); MAB = U(C4['MAB']['value'])
H33 = U(lc['beta_lead_split']['lemma']['H33'])
assert H33.upper() >= arb(lc['v4b_H3j']['H33'][1]).upper() * (1 - arb('1e-15'))
r1 = U(J['bloop']['alpha2_remainders']['r1']); r2 = U(J['bloop']['alpha2_remainders']['r2']); rn = (r1 * r1 + r2 * r2).sqrt()
PI3 = arb('0.01607'); Pi3T3 = U(lc['Pi3']['treatments']['T3_sym_lemmaG33_eps1'])
slack = PI3 - Pi3T3
assert abs(lo(slack) - lc['v4b_A2_budget']['cases']['2.5e-4']['slack']) < 1e-15        # same slack as fin13v4b (one ulp)
cases = {}
for a in ('1.6e-4', '2.5e-4'):
    al = arb(a)
    base_parts = {k: U(v) for k, v in lc['alpha_form']['need'][a]['parts'].items()}
    base = sum(base_parts.values(), arb(0)) + CDR                                # fin13v4b need (e0_lemma_consts section (10))
    extra = dict(lead_mismatch_ab=(DAB + DLEAD) * (GABN + LOOPC),                # A-1
                 loop_third_order=H33 * al * rn * MAB * N3P * al)                 # A-1b
    need = base + sum(extra.values(), arb(0))
    exc = (need - 25) * al
    assert need.lower() > 25 and slack.lower() > exc.upper(), (a, need, exc, slack)
    exc0 = (base - 25) * al
    cases[a] = dict(need_v4b=up(base), extra={k: up(v) for k, v in extra.items()}, need_total=up(need), excess=up(exc), slack=lo(slack),
                    margin=lo(slack / exc), margin_v4b_form=lo(slack / exc0), ok=True)
out = dict(cases=cases, criterion='Pi3_T3 + (need - 25) alpha <= 0.01607 at alpha = AMAX_E0 = 2.5e-4 (need nondecreasing in alpha); PI3, MODA unchanged',
           review_expected=dict(margin_2p5e4=1.6569484534524035, need_2p5e4=38.723161368864396, extra=dict(lead_mismatch_dab=0.17938241628366314, loop_r_times_MAB=0.0026405498006482025)),
           inputs_widening='floats from JSON widened one ulp (upper up, slack lower)')
out['args'] = VC.args_block(list(P.values()), extra=dict(runtime_s=round(time.time() - t0, 2)))
VC.dump(os.path.join(VC.OUT, 'e0_A2_budget_v4c.json'), out)
print(json.dumps(cases, indent=1))
