import json, math
def ub(v):
    if v == 0: return 0.0
    e = math.floor(math.log10(abs(v))) - 3
    return float(f"{math.ceil(v / 10**e) * 10**e:.4g}")
def ubd(d): return {k: (ubd(v) if isinstance(v, dict) else ub(v)) for k, v in d.items()}
N = json.load(open('v3_n3d_0.001_0.25.json')); Rr = json.load(open('v3_res_0.25.json')); G = json.load(open('v3_f23_0.5.json'))
out = {"setting": {"P": 2, "alpha_max": 1e-3, "region": "[10.25, R]", "functions": "C^4 versions psi3~, h2~ (c4fun.py); C^3 versions have u'''' jumps 1e-7..4e-6 at all 20 interfaces (delta masses in derivative-level quantities)",
                   "zeta": "box sup of zeta (inner ramp S((x-10.25)/6), outer ramp S((R-x)/3)) in all 'zeta' norms"}}
for R in ("40", "60"):
    o = {}
    n = N[f"R={R}"]
    o["F1a_C3W_deriv_zeta"] = ub(n["C3W_deriv_zeta"]); o["F1a_C3W_deriv_nozeta"] = ub(n["C3W_deriv_nozeta"]); o["F1a_sup_c3prime"] = ub(n["sup_c3prime"])
    o["F1a_C3W_level0_for_C4_functions"] = ub(n["C3W_level0_C4"])
    o["F1a_per_region_norm_zeta"] = ubd(n["per_region_norm"]); o["F1a_per_region_sup_c3prime"] = ubd(n["per_region_sup_c3prime"])
    for w, nm in (("psi", "psi3"), ("h2", "h2")):
        for lev in (0, 1):
            r = Rr[f"{w}_level{lev}_R={R}"]
            o[f"F1b_{nm}_residual_level{lev}" + ("_zeta" if lev else "")] = ub(r["norm"]); o[f"F1b_{nm}_level{lev}_per_region"] = ubd(r["per_region_norm"])
            o[f"F1b_{nm}_level{lev}_max_weighted_pointwise"] = [ub(r["max_weighted_pointwise"]), r["at"]]
    g = G[f"R={R}"]
    o["F3_norms_zeta"] = {k: ub(v["norm_zeta"]) for k, v in g.items()}; o["F3_norms_nozeta"] = {k: ub(v["norm_nozeta"]) for k, v in g.items()}
    o["F2_sup"] = ubd(G[f"F2_sup_R={R}"])
    out[f"R={R}"] = o
out["F3_sup_per_region"] = {r: ubd(v) for r, v in G["sup_per_region"].items()}
json.dump(out, open('int_constants_v3.json', 'w'), indent=1)
for R in ("R=40", "R=60"):
    print(R); [print(' ', k, v) for k, v in out[R].items() if not k.endswith("per_region") and 'pointwise' not in k]
print(out["R=60"]["F1a_per_region_norm_zeta"], out["R=60"]["F1b_psi3_level1_per_region"], out["R=60"]["F1b_h2_level1_per_region"])
print(out["F3_sup_per_region"])
