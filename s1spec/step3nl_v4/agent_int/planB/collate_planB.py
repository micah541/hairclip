import json, math
def ub(v):
    if v == 0 or abs(v) < 1e-280: return v
    if isinstance(v, str): return v
    e = math.floor(math.log10(abs(v))) - 3
    return float(f"{math.ceil(v / 10**e) * 10**e:.4g}")
def ubd(d): return {k: (ubd(v) if isinstance(v, dict) else (ub(v) if isinstance(v, (int, float)) else v)) for k, v in d.items()}
L = json.load(open('planB_lint.json')); S = json.load(open('planB_sob.json')); R = json.load(open('planB_res_1.json'))
T3 = json.load(open('planB_T3_0.001_0.01_0.01_1.json'))['sup']; MOD = json.load(open('planB_mod_1.json')); RAMP = json.load(open('planB_ramp.json'))
N = {a: json.load(open(f'planB_n3_{a}_1.json')) for a in ('0.001', '0.0007', '0.0005', '0.0003', '0.0001')}
out = {"design": "Plan B: P=2 (x<=40), 2-S((x-40)/20) on [40,60], 1 (x>=60), no rho blend; zeta: inner S((x-10.25)/6), outer S((260-x)/10), R=260; level-2 weight factor 1+(x/60)^4; |alpha|<=1e-3, jet box 1e-2",
       "MU1": {"global_min": ub(L["global_min_mu1"][0] if False else math.floor(L["global_min_mu1"][0]*1e5)/1e5), "at_x": L["global_min_mu1"][1],
               "per_region": {r: math.floor(v["mu1"]*1e5)/1e5 for r, v in L["regions"].items()}, "asymptotic_(2P-1)/8": 0.125, "profile_cells_2.5": {k: math.floor(v["mu1"]*1e4)/1e4 for k, v in L["cell_profile"].items()}},
       "MU2": {"plain_W_global_min": math.floor(L["global_min_mu2"][0]*1e5)/1e5, "level2_weight_global_min": math.floor(L["global_min_mu2w"][0]*1e5)/1e5, "level2_weight_at_x": L["global_min_mu2w"][1],
               "per_region_plainW": {r: math.floor(v["mu2"]*1e5)/1e5 for r, v in L["regions"].items()}, "per_region_level2weight": {r: math.floor(v["mu2w"]*1e5)/1e5 for r, v in L["regions"].items()},
               "e1_sup": ub(L["global_max_e1"]), "e0_sup": ub(L["global_max_e0"]), "note": "mu2 = mu1 + 1/2 (plain W); with m2 = 1+(x/60)^4 mu2 -> ~1/8 at x=260"},
       "KW": S["W"]["KW"], "PSI1W": S["W"]["Psi1"], "C1W": S["W"]["c1"], "psi1_min_lower_bound": S["W"]["psi1_min_lower_bound"],
       "level2_weight_sobolev": S["W2"], "cutoff_derivative_sups": S["cutoff"], "W_values": S["W_values"],
       "FnW": {"F1W_psi3_res_W (coeff of alpha)": R["psi_level0"]["norm"], "F2W_h2_res_W (coeff of alpha^2)": R["h2_level0"]["norm"],
               "C3W_level0 (alpha^3), alpha_max=1e-3": N["0.001"]["C3W_level0_C4"], "C3W_level0 by alpha_max": {a: ub(v["C3W_level0_C4"]) for a, v in N.items()},
               "C3W_deriv_zeta (alpha^3), alpha_max=1e-3": N["0.001"]["C3W_deriv_zeta"], "C3W_deriv_zeta by alpha_max": {a: ub(v["C3W_deriv_zeta"]) for a, v in N.items()},
               "psi3_res_deriv_zeta": R["psi_level1"]["norm"], "h2_res_deriv_zeta": R["h2_level1"]["norm"],
               "per_region": {"psi3_res": R["psi_level0"]["per_region_norm"], "h2_res": R["h2_level0"]["per_region_norm"], "psi3_res_deriv": R["psi_level1"]["per_region_norm"], "h2_res_deriv": R["h2_level1"]["per_region_norm"],
                              "C3W_deriv_zeta": N["0.001"]["per_region_norm_deriv_zeta"], "sup_level0_c3_bound": N["0.001"]["per_region_sup_level0_bound"], "sup_c3prime": N["0.001"]["per_region_sup_c3prime"]}},
       "T3_per_region": {r: ubd(v) for r, v in T3.items()},
       "MOD_norms_zeta": {k: ub(v["norm_zeta"]) for k, v in MOD["norms"].items()}, "MOD_sup_all": ubd(MOD["sup_all"]), "MOD_sup_per_region": {r: ubd(v) for r, v in MOD["sup_per_region"].items()},
       "RAMP_linear_deviation": {"region": "[10.25,16.25]", "table": ubd(RAMP["table"]), "alpha_split": ubd(RAMP["alpha_split"])}}
out["outer_ramp_check"] = {"W(260)/W(10.25)": S["W_values"]["W(260)/W(10.25)"], "W(260)": S["W_values"]["W(260)"], "sqrt(W(250))": ub(math.sqrt(6.426e-9)), "coeff_sups_[250,260]": ubd(T3["[250,260]"]),
                            "int_W_250_260": S["W_values"]["integral_W_250_260_upper"], "C3W_deriv_zeta_[250,260]": N["0.001"]["per_region_norm_deriv_zeta"]["[250,260]"]}
json.dump(out, open('planB_constants.json', 'w'), indent=1); print(json.dumps(out["FnW"], indent=0)[:1500]); print(out["MU1"]["global_min"], out["MU2"])
