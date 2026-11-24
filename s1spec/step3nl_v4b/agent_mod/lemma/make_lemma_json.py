import json
out = {"main": json.load(open('mod_lemma_eta1e-2.json')), "eta1e-3": json.load(open('mod_lemma_eta1e-3.json')), "eta1e-4": json.load(open('mod_lemma_eta1e-4.json')),
       "inputs": {"lem_integrals": json.load(open('lem_integrals.json')), "lyap_extra": json.load(open('lyap_extra.json'))},
       "quantity_names": {"A": "||w~||", "W1": "||w~'||", "B": "||w~''||", "W3": "||w~'''||", "D4": "||w~''''||", "S2B": "sup|w''| * B", "S2W3": "sup|w''| * W3", "S2S2B": "sup|w''|^2 * B",
                          "Phi": "||chi Phi_F||", "S0..S4": "sup over sigma in [29,32] of |w^(i)|, i = 0..4 (raw w, not w~)", "S2S2": "S2^2", "S2S3": "S2 S3", "S2S2S2": "S2^3"}}
json.dump(out, open('mod_lemma.json', 'w'), indent=1)
for tag in ("main", "eta1e-3", "eta1e-4"):
    m = out[tag]; print("==", tag, m["setting"]["jet_box_eta"])
    for n in ("a", "b", "beta"): print(n, {q: float(f"{m['t_coeffs'][n][q]:.3g}") for q in ("A", "W1", "B", "W3", "D4", "S2B", "S2W3", "S2S2B", "Phi")}, {q: f"{m['t_coeffs'][n][q]:.2e}" for q in ("S0", "S1", "S2", "S3", "S4", "S2S2", "S2S3", "S2S2S2")})
    print("vec", {q: float(f"{m['t_coeffs_vec_ab'][q]:.3g}") for q in ("A", "W1", "B", "W3", "D4", "S2B", "S2W3", "S2S2B", "Phi")})
    nn = m["Neumann"]; print("nu", nn["nu_bar"], "n0", nn["n0(alpha)"], "n1", nn["n1(W1)"], "n2", nn["n2(alpha*A)"], "n3", nn["n3(S1)"], "n4", nn["n4(alpha*S0)"])
    print("tau", {q: float(f"{v:.3g}") for q, v in nn["tau_q(max_k t_kq/d_k)"].items() if v > 1e-6}, "d", nn["weights_d"])
m = out["main"]; print("lead", m["lead_coeff_of_A"], "G", m["G"], "G33_lower", m["G33_lower"], "||G^-1||inf", m["norm_G_inverse_inf"])
print("H", m["H=(G^-1)^T"])
# linear corollary: only the leading (Pi2 L* f_j, eps1) terms in t, E-feedback kept
m = out["main"]; nn = m["Neumann"]; d = nn["weights_d"]; lead = [m["lead_coeff_of_A"][k] for k in ("a", "b", "beta")]
tau_lin = max(l / x for l, x in zip(lead, d)); fac = nn["1/(1-nu)"]
cor = {"tau_lin": tau_lin, "A_coeff_at_alpha0": lead, "|(a,b)|_A_coeff": m["lead_coeff_of_A"]["|(a,b)|"],
       "alpha_A_coeff c1 (|x_i| <= lead_i A + c1_i |alpha| A + c2_i W1 A)": [n * tau_lin * fac for n in nn["n0(alpha)"]],
       "W1_A_coeff c2": [n * tau_lin * fac for n in nn["n1(W1)"]]}
out["linear_corollary"] = cor; json.dump(out, open('mod_lemma.json', 'w'), indent=1); print(cor)
