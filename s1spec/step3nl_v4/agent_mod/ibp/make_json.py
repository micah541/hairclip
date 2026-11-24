import json, os
out = {"description": "IBP modulation lemma constants (F3). Forms: |x_i| <= (lead_i + m_total_i |alpha|) A + q_total_i (A+W1+B)^2 + cPhi_total_i ||chi Phi_F|| + strip; A=||w~||, W1=||w~'||, B=||w~''|| (L2(rho)), jet box |w|,|w'|<=eta on supp chi",
       "symbolic_identities": json.load(open('verify_symbolic.json')),
       "float_crosscheck": [json.load(open(f'float_check_{s}.json')) for s in (1, 2, 3)]}
for tag in ('5e-4', '1e-3'):
    f = f'final_{tag}.json'
    if os.path.exists(f): out[f'alpha_max_{tag}'] = json.load(open(f))
json.dump(out, open('mod_ibp.json', 'w'), indent=1)
