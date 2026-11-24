"""Assemble mod_forms.json from gram8_g1s4.json (rigorous Gram) and forms_gram8_g1s4.json (mod_forms.py); consistency checks vs lyap_consts2."""
import json
from flint import arb
import mod_forms as MF
G8 = MF.load('gram8_g1s4.json'); R = json.load(open('gram8_g1s4.json')); Fm = json.load(open('forms_gram8_g1s4.json'))
names = R["names"]
M = [[G8[6][6], G8[6][7]], [G8[7][6], G8[7][7]]]
det = M[0][0] * M[1][1] - M[0][1] * M[1][0]; Mi = [[M[1][1] / det, -M[0][1] / det], [-M[1][0] / det, M[0][0] / det]]
tr = M[0][0] + M[1][1]; disc = (tr * tr - 4 * det).sqrt(); lmin = (tr - disc) / 2
Minv_lyap = float(json.load(open('/home/claude/s1spec/step3/lyap_consts2.json'))["Minv"])
Lf = (G8[0][0] + G8[1][1]).sqrt(); Lf_lyap = float(json.load(open('/home/claude/s1spec/step3/lyap_consts2.json'))["Lf"])
iv = lambda x: [float(x.lower()), float(x.upper())]
out = {
 "R1_upper": max(Fm["Q1"]["lam_max_upper"], -Fm["Q1"]["lam_min_lower"]),
 "R1_enclosure_of_sup|Q1|": [max(Fm["Q1"]["lam_max_lower"], -Fm["Q1"]["lam_min_upper"]), max(Fm["Q1"]["lam_max_upper"], -Fm["Q1"]["lam_min_lower"])],
 "Q1_sup_upper": Fm["Q1"]["lam_max_upper"], "Q1_inf_lower": Fm["Q1"]["lam_min_lower"],
 "L2max_upper": Fm["Q2"]["lam_max_upper"], "L2max_lower_witness": Fm["Q2"]["lam_max_lower"],
 "L2min_enclosure": [Fm["Q2"]["lam_min_lower"], Fm["Q2"]["lam_min_upper"]],
 "float_reference_exact_derivatives": {"Q1_eigs": Fm["Q1"]["float_eigs"], "Q2_eigs": Fm["Q2"]["float_eigs"]},
 "vector_order": names, "Gram_8x8": R["G"],
 "M": [[iv(M[i][j]) for j in range(2)] for i in range(2)], "Minv": [[iv(Mi[i][j]) for j in range(2)] for i in range(2)],
 "lambda_min_M": iv(lmin), "Minv_scalar_upper(1/lmin)": float((1 / lmin.lower()).upper()), "Minv_lyap_consts2": Minv_lyap,
 "norms": {names[i]: iv(G8[i][i].sqrt()) for i in range(8)},
 "||(Lf1,Lf2)||": iv(Lf), "||(Lf1,Lf2)||_lyap_consts2": Lf_lyap,
 "||(Df1,Df2)||": iv((G8[4][4] + G8[5][5]).sqrt()), "||(gA,gB)||": iv((G8[2][2] + G8[3][3]).sqrt()),
 "runtime_s": {"gram_group1_nsub4": R["runtime_s"], "forms": "<2"},
 "notes": ["Gram via 2nd-order Taylor models, 51246 tight-profile boxes (x<=tau_tip) + 2400x4 chunk-bound sub-boxes (sigma in [8,32]).",
           "Float references quoted in the task (mask x>0.3, finite differences) UNDERESTIMATE: with the mask removed (or exact ODE derivatives) the float eigenvalues are Q1 [-13.55, 11.33], Q2 [-0.189, 43.32].",
           "Gram enclosure radii ~1e-6 relative; G4 normalised smallest eigenvalue 2e-5 (Q1), 7e-6 (Q2) handled by the Loewner majorant G+ = mid + diag(rad row sums).",
           "psi2 = -sin(theta), cos(theta) = +St*C - Ct*S in the Lyapunov-run frame (checked against float theta; the sign quoted in the task is the opposite)."]}
json.dump(out, open('mod_forms.json', 'w'), indent=1)
for k in ("R1_upper", "R1_enclosure_of_sup|Q1|", "L2max_upper", "L2min_enclosure", "lambda_min_M", "Minv_scalar_upper(1/lmin)", "Minv_lyap_consts2", "||(Lf1,Lf2)||", "||(Lf1,Lf2)||_lyap_consts2", "||(Df1,Df2)||", "||(gA,gB)||", "norms"):
    print(k, out[k])
