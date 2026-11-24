"""Add-on: linear-deviation coefficients on the INNER RAMP [10.25, 16.25] (u0 = alpha psi3~ + alpha^2 h2~, C^4 jets):
   d001 = sup|1 - J^-4|   (w'''' coefficient),  d010 (w'''),  d100 (w''),  A1 (w'),  A0 (w)     [coefficients of D = F(u0+w) - F(u0) - L w, linear part in w]
 tables over alpha_max in {2.5e-4, 5e-4, 1e-3} and jet box eta = |w| = |w'| in {0, 1e-3, 1e-2}.  (sup over the alpha-box and the w-box.)
 alpha-split: c(alpha_max, eta=0) <= a alpha_max + b alpha_max^2 holds at the tabulated alpha_max with a := c(2.5e-4,0)/2.5e-4, b := max_tab (c - a alpha)/alpha^2
 (rigorous at the tabulated alpha_max only; for alpha in between use the next larger tabulated alpha_max since the sup is monotone);  w-part := c(alpha_max, eta) - c(alpha_max, 0)."""
import sys, json
argv0 = sys.argv[:]; sys.argv = ['x', '1e-3', '1e-2', '1e-2', '1']
import pb_nlin as M
from pb_common import *
from c4fun import PSI, H2
KEYS = {"d001_w4": (0, 0, 1), "d010_w3": (0, 1, 0), "d100_w2": (1, 0, 0), "A1_w1": "A1", "A0_w0": "A0"}
bl = make_segments([('x', 10.25), ('x', 16.25)], [0.005])
jets = []
for sa, sb in bl:
    b = make_box(sa, sb); y = PSI.jets_any(b.X, b.xc, 5); h = H2.jets_any(b.X, b.xc, 5); jets.append((b, y, h))
tab = {}
for am in (2.5e-4, 5e-4, 1e-3):
    for eta in (0.0, 1e-3, 1e-2):
        M.A = arb(0, am); M.E0b = arb(0, eta); M.E1b = arb(0, eta)
        sup = {k: 0.0 for k in KEYS}; sup_u1 = 0.0
        for b, y, h in jets:
            d = M.box_vals(y, h, b.p, b.q, b.k, b.m)
            for name, key in KEYS.items(): sup[name] = max(sup[name], up(d[key]))
            sup_u1 = max(sup_u1, up(d["|u0'|"]))
        tab[f"alpha={am:g},eta={eta:g}"] = dict(sup, **{"sup|u0'|": sup_u1})
        print(am, eta, {k: float('%.4g' % v) for k, v in sup.items()}, flush=True)
split = {}
for name in KEYS:
    a = tab["alpha=0.00025,eta=0"][name] / 2.5e-4
    b_ = max((tab[f"alpha={am:g},eta=0"][name] - a * am) / am**2 for am in (2.5e-4, 5e-4, 1e-3))
    wpart = {f"alpha={am:g},eta={eta:g}": tab[f"alpha={am:g},eta={eta:g}"][name] - tab[f"alpha={am:g},eta=0"][name] for am in (2.5e-4, 5e-4, 1e-3) for eta in (1e-3, 1e-2)}
    split[name] = dict(a=a * (1 + 1e-12), b=max(b_, 0.0) * (1 + 1e-12), w_part=wpart)
json.dump({"table": tab, "alpha_split": split, "region": "[10.25,16.25]", "box_width_x": 0.005}, open('planB_ramp.json', 'w'), indent=1)
print(json.dumps(split, indent=0)[:1800])
