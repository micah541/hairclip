"""collate the final runs into int_constants.json (rigorous upper bounds, outward-rounded to 4 significant digits upward)"""
import json, math
def r4(v):
    if v == 0: return 0.0
    e = math.floor(math.log10(abs(v))) - 3
    return math.ceil(v / 10**e) * 10**e if v > 0 else -math.ceil(-v / 10**e) * 10**e   # upward for positive; for negative: towards zero is NOT conservative, see below
def ub(v):   # round UP (valid for upper bounds, also for negative numbers)
    if v == 0: return 0.0
    e = math.floor(math.log10(abs(v))) - 3
    return math.ceil(v / 10**e) * 10**e
T1 = json.load(open('int_T1_0.001_0.25.json'))
T2 = json.load(open('int_T2_W_0.25.json'))
T3 = json.load(open('int_T3_0.001_0.01_0.01_0.5.json'))['sup']
T4 = json.load(open('int_T4_0.25.json'))['sup']
out = {"setting": {"alpha_max_T1_T3": 1e-3, "jet_box_T3": {"|w|": 1e-2, "|w'|": 1e-2}, "weight": "W=exp(-psi)=(1+x^2)^-3 for x<=30.25, blend to rho-type; half-line norms",
                   "box_widths": "T1: sigma 0.005/0.00625/0.00625/0.0125; T2: 0.00125-0.005; T3: 0.01/0.0125/0.0125/0.025; T4: 0.0025/0.005"}}
regs = ["[12,20]", "[20,30.25]", "[30.25,40]", "[40,60]"]
for Rr, nreg in (("40", 3), ("60", 4)):
    R = f"R={Rr}"; o = {}
    t1 = T1[R]
    o["T1_C3W"] = ub(t1["C3W"]); o["T1_sup_g3over6"] = ub(t1["sup_c3"])
    o["T1_per_region_norm_sq"] = {k: ub(v) for k, v in t1["per_region_sq"].items()}
    o["T1_per_region_sup_g3over6"] = {k: ub(v) for k, v in t1["per_region_sup_c3"].items()}
    o["T2_psi3_res_W"] = ub(T2[f"psi_R={Rr}"]["norm"]); o["T2_h2_res_W"] = ub(T2[f"h2_R={Rr}"]["norm"])
    o["T2_psi3_per_region_norm"] = {k: ub(v) for k, v in T2[f"psi_R={Rr}"]["per_region_norm"].items()}
    o["T2_h2_per_region_norm"] = {k: ub(v) for k, v in T2[f"h2_R={Rr}"]["per_region_norm"].items()}
    o["T2_psi3_max_weighted_pointwise"] = [ub(T2[f"psi_R={Rr}"]["max_weighted_pointwise"]), T2[f"psi_R={Rr}"]["at"]]
    o["T2_h2_max_weighted_pointwise"] = [ub(T2[f"h2_R={Rr}"]["max_weighted_pointwise"]), T2[f"h2_R={Rr}"]["at"]]
    keys = set().union(*[T3[r].keys() for r in regs[:nreg]])
    o["T3_sup_over_[12,R]"] = {k: ub(max(T3[r][k] for r in regs[:nreg])) for k in keys}
    o["T4_sup_over_[12,R]"] = {k: ub(max(T4[r][k] for r in regs[:nreg])) for k in ("-mu2", "|e1|", "|e0|", "-mu1")}
    out[R] = o
out["T3_per_region"] = {r: {k: ub(v) for k, v in T3[r].items()} for r in regs}
out["T4_per_region"] = {r: {k: ub(v) for k, v in T4[r].items()} for r in regs}
out["T2_parts"] = {w: T2[w]["max_weighted_parts_by_region"] for w in ("psi_R=60", "h2_R=60")}
json.dump(out, open('int_constants.json', 'w'), indent=1)
print(json.dumps({k: out[k] for k in ("R=40", "R=60")}, indent=1))
