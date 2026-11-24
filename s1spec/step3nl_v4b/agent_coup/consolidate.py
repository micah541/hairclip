import json
F = json.load(open('float_ref.json')); out = {"description": "rigorous UPPER bounds (ball arithmetic) unless noted; 'float' = independent float reference",
 "tip_interval": None, "P": {}}
fm = {'sB': 'sB', 'sA': 'sA', 'a_f': 'af', 'b_f': 'bf', 'IA': 'IA', 'K_rho_16.25': 'K_rho_16.25', 'Phi_34.25': 'Phi_34.25', 'drift_f': 'drift_f'}
for j in range(4): fm[f'd{j}'] = f'd{j}'; fm[f'c{j}'] = f'c{j}'
for P in ('2', '1.5'):
    R = json.load(open(f'coup_P{P}.json')); out['tip_interval'] = R['tip']
    ent = {"rigorous": R['quantities'], "float": {k: F[P][v] for k, v in fm.items()}, "outer": {}, "mu_window": R['mu_window'],
           "runtime_s": R['runtime_s'], "n_pieces": {k: v['n_pieces'] for k, v in R['n_pieces'].items()}}
    for Rr, o in R['outer'].items():
        fo = {'sup|A_zo|': f'Azo_{Rr}', 'sup|B_zo|': f'Bzo_{Rr}', 'W(R-3)': f'W{int(Rr)-3}', 'W(R)': f'W{Rr}'}
        fo.update({f'sup|coef{j}|W^1/2': f'oc{j}_{Rr}' for j in range(4)})
        ent['outer'][Rr] = {"rigorous": o, "float": {k: F[P][v] for k, v in fo.items()}}
    ent['mu_window']['float_design2'] = {"mu_min": F[P]['mu_min'], "mu_max": F[P]['mu_max'], "ratio": F[P]['ratio']}
    out['P'][P] = ent
json.dump(out, open('coupling_constants.json', 'w'), indent=1)
