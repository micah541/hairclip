"""Step 3 amendment v3.2 (draft): collect every new constant and margin into ../v32_constants.json (inputs: v32_cases.json, v32_h6.json,
ramp_budget.json, scratch_lemmaS Bbox runs).  Run after v32_cases.py, v32_h6.py, the two Bbox runs, ramp_out.py and ramp_budget.py."""
import json
H = '/home/claude/s1spec/writeup/step3_amend_v32'
C = json.load(open(f'{H}/scripts/v32_cases.json')); R = json.load(open(f'{H}/scripts/v32_h6.json')); B = json.load(open(f'{H}/scripts/ramp_budget.json'))
a = json.load(open(f'{H}/scratch_lemmaS/lemmaS_Bbox_1e-05_a0.00015_dv0.15.json')); b = json.load(open(f'{H}/scratch_lemmaS/lemmaS_Bbox_1e-05_a0.00015_dv0.15_frozen.json'))
out = dict(status='adopted as v3.2 (12 Nov 2026); outer-ramp Lemma R open', version='step 3 amendment v3.2',
  families=dict(cone=dict(centres=[200, 380], unit='r/4', X=3, M=0.12, K2=0.18, K3=0.18, changed=False, lemmaO=C['cone X=3 M=0.12 K=0.18']),
                far=dict(centres_from=376, unit=113, X=2, M=0.15, K2=0.3, K3=0.3, changed=True, lemmaO=C['far X=2 M=0.15 K=0.3'],
                         v31=dict(M=1.0, K=1.5, lemmaO=C['far X=2 M=1.0 K=1.5']))),
  H5_data_thresholds=R['thresholds'], far_physical_thresholds_lambda1_units=dict(R['far_physical_thresholds_lambda1_units'], H4_int_h5_sq_cyl_units=0.0225),
  H6_closure=dict(slope=R['slope'], cone={k: v for k, v in R['cone'].items()}, far=R['far'], closure_margins=R['closure_margins'],
                  min_closure_margin=R['min_closure_margin'], data=R['data'], data_vs_thresholds=R['data_vs_thresholds']),
  windows_frozen_profile=[dict(yE=x['yE'], box=y['box'], gamma_prime=y['gp'], EB_v31=x['EB'], EB_frozen=y['EB'], forcing_v31=x['forcing'], forcing_frozen=y['forcing'],
                               derived_v31=x['derived'], derived_frozen=y['derived'], K_v31=x['K'], K_frozen=y['K']) for x, y in zip(a, b)],
  residual_constants=json.load(open(f'{H}/scratch_lemmaS/cres.json')),
  outer_ramp_budget_FLOAT_ESTIMATE=B,
  stale_values=dict(note='step-3 v3.1 numbers computed at the withdrawn ratio 1.2e-7; values at 1.6001e-7 with the stored w = 0.317 (scaled by sqrt(1.6001/1.2))',
                    C0sqrt_t=[0.0021, 0.00237], C1sqrt_t=[0.0063, 0.00726], C2sqrt_t=[0.034, 0.03946], sqrt_ceilX_C3t=[0.0401, 0.04623],
                    H5_cone_threshold=[0.14, 0.1337], slope_margin=[0.1237, 0.1227], recovered_hxx=[0.154, 0.1595], recovered_slope=[0.126, 0.1273]))
json.dump(out, open(f'{H}/v32_constants.json', 'w'), indent=1); print('written', len(json.dumps(out)))
