import json, math
def ub(v):
    if v == 0: return 0.0
    e = math.floor(math.log10(abs(v))) - 3
    return float(f"{math.ceil(v / 10**e) * 10**e:.4g}")
def ubd(d): return {k: (ubd(v) if isinstance(v, dict) else ub(v)) for k, v in d.items()}
T1 = json.load(open('v2_T1_0.001_0.25.json')); T2 = json.load(open('v2_T2_W_0.25.json'))
T3 = json.load(open('v2_T3_0.001_0.01_0.01_0.5.json'))['sup']; V56 = json.load(open('v2_V56_1.json'))
regs = ["[10.25,16.25]", "[16.25,20]", "[20,30.25]", "[30.25,40]", "[40,60]"]
out = {"design": "v2: zeta=0 for x<=10.25, C3 smoothstep to 1 at 16.25, outer ramp [R-3,R]; psi'=(1-s)2P x/(1+x^2)+s x(1+x^2)^(-1/3)/3, W=(1+x^2)^-P for x<=30.25",
       "settings": {"alpha_max": 1e-3, "jet_box": 1e-2, "region": "[10.25,R], sigma >= 8.0009 everywhere", "zeta_in_norms": "V5 only (box sup of zeta); V1, V2 are plain W-norms"}}
for pk in ("2", "1.5"):
    o = {}
    for Rr, nreg in (("40", 4), ("60", 5)):
        R = f"R={Rr}"; r = {}
        t1 = T1[f"P={pk}_{R}"]; r["V1_C3W"] = ub(t1["C3W"]); r["V1_sup_g3over6"] = ub(t1["sup_c3"])
        r["V1_per_region_norm"] = {k: ub(math.sqrt(v)) for k, v in t1["per_region_norm_sq"].items()}
        r["V1_per_region_sup_g3over6"] = ubd(t1["per_region_sup_c3"])
        for w, nm in (("psi", "V2_psi3_res"), ("h2", "V2_h2_res")):
            t2 = T2[f"P={pk}_{w}_{R}"]; r[nm] = ub(t2["norm"]); r[nm + "_per_region"] = ubd(t2["per_region_norm"])
            r[nm + "_max_weighted_pointwise"] = [ub(t2["max_weighted_pointwise"]), t2["at"]]
        keys = T3[regs[0]].keys()
        r["V3_sup"] = {k: ub(max(T3[g][k] for g in regs[:nreg])) for k in keys}
        r["V4_sup"] = None
        v5 = V56["V5"][f"P={pk}_{R}"]; r["V5_norms_zeta"] = {k: ub(v["norm_zeta"]) for k, v in v5.items()}
        r["V5_norms_nozeta"] = {k: ub(v["norm_nozeta"]) for k, v in v5.items()}
        r["V6"] = ubd(V56["V6"][f"P={pk}_{R}"]) if False else V56["V6"][f"P={pk}_{R}"]
        t4 = json.load(open(f'v2_T4_P{pk}_0.25.json'))['sup']
        r["V4_sup"] = {k: ub(max(t4[g][k] for g in regs[:nreg])) for k in ("-mu2", "|e1|", "|e0|", "-mu1")}
        o[R] = r
    t4 = json.load(open(f'v2_T4_P{pk}_0.25.json'))['sup']; o["V4_per_region"] = {g: ubd(t4[g]) for g in regs}
    out[f"P={pk}"] = o
out["V3_per_region"] = {g: ubd(T3[g]) for g in regs}; out["V5_sup_per_region"] = {g: ubd(v) for g, v in V56["V5"]["sup_per_region"].items()}
out["V2_parts_P2_R60"] = {w: T2[f"P=2_{w}_R=60"]["max_weighted_parts_by_region"] for w in ("psi", "h2")}
json.dump(out, open('int_constants_v2.json', 'w'), indent=1)
for pk in ("2", "1.5"):
    for R in ("R=40", "R=60"):
        r = out[f"P={pk}"][R]; print(pk, R, {k: r[k] for k in ("V1_C3W", "V1_sup_g3over6", "V2_psi3_res", "V2_h2_res", "V4_sup")}); print('  V5', r["V5_norms_zeta"]); print('  V6', {k: ub(v) if not isinstance(v, list) else v for k, v in r["V6"].items()})
print(out["V3_per_region"]["[10.25,16.25]"])
for pk in ("2", "1.5"): print(pk, out[f"P={pk}"]["R=60"]["V1_per_region_norm"], out[f"P={pk}"]["R=60"]["V2_psi3_res_per_region"], out[f"P={pk}"]["R=60"]["V2_h2_res_per_region"])
