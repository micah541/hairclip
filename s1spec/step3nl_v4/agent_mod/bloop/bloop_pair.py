"""Signed pairings <v, f3>, v in (Lf1, Lf2, gA, gB, g0 = chi psi3^, h2c = chi h2^, f1, f2), f3 = chi F - a_c eta1 - b_c eta3.
Tight Taylor-model integrals of rho*v (bg2.py = mod_gram2 machinery, same boxes) times the hull of F (outA5.txt rows):
  tip boxes (prof_tight, x <= tau_tip): F-hull = union of hulls of the MC rows whose x-range meets the box;
  far boxes: the MB row's ws0 hull.  Error terms: Fm * int(rho v) + rad(F) * int|rho v| (abs via hull sup).  Tip bump -a eta1 - b eta3 added by L1 bound."""
import sys, json, time
from flint import arb, arb_series
sys.argv_saved = sys.argv; sys.argv = sys.argv[:1]
import bg2 as B
sys.argv = sys.argv_saved
MG, HR = B.MG, B.HR
import gen_hs3 as G3
BUMP = (0.8324 * 9.0923e-8 + 0.4245 * 2.5707e-7) * 1.01
NAMES = ["Lf1", "Lf2", "gA", "gB", "g0", "h2c", "f1", "f2"]
ab = lambda z: arb(abs(z).upper())
rowsA = []   # (xlo, xhi, Fhull, tag)
CtB, StB, Ct2h, St2h = MG.frame_data(); T0 = arb(G3.T0)
for line in open('/home/claude/s1spec/step3/outA5.txt'):
    tag = line[:2]
    if tag not in ("MB", "MC"): continue
    v = list(map(float, line.split()[1:])); ta, tb = v[0], v[1]
    h = [MG.A(v[2 + 2 * i], v[3 + 2 * i]) for i in range((len(v) - 2) // 2)]
    if tag == "MC":
        st = h[4]; cc = (8 - st) / T0; Xv = MG.hull(cc * (T0 - arb(tb)), cc * (T0 - arb(ta)))
        rowsA.append((float(Xv.lower()), float(Xv.upper()), h[9], tag, ta >= float(T0) - G3.LB - 1e-12))
    else: rowsA.append((None, None, h[5], tag, False))
MBrows = [r for r in rowsA if r[3] == "MB"]; MCrows = [r for r in rowsA if r[3] == "MC"]
state = {"phase": "far", "idx": 0, "Xb": None}
res = {n: arb(0) for n in NAMES}; resmid = {n: arb(0) for n in NAMES}; l1b = {n: arb(0) for n in NAMES}
_psi = B.psi_series
def psi_hook(Xb, n=3):
    state["Xb"] = Xb; return _psi(Xb, n)
B.psi_series = psi_hook
_fun = B.functions
def fun_hook(st, chi, Xs, yser, n=B.N, nout=3):
    out, rho = _fun(st, chi, Xs, yser, n, nout)
    Xball = Xs.coeffs()[0]; J = HR.jets_any(HR.Hf, Xball, float(Xball.mid()))
    h = arb_series([J[0], J[1], J[2] / 2], prec=3)
    out["h2c"] = B.tr(B.tr(chi, nout) * h, nout); return out, rho
B.functions = fun_hook
def acc_hook(G, uc, uh, Rc, Rh, L, I1, I2):
    if state["phase"] == "far":
        row = MBrows[state["idx"]]; state["idx"] += 1; Fh = row[2]; bump = False
        # chi in strip is inside the u's only for gA..; f3 carries chi too: multiply by hull [0,1] on strip rows (rho <= 3e-11, negligible)
    else:
        Xb = state["Xb"]; lo, hi = float(Xb.lower()) - 1e-7, float(Xb.upper()) + 1e-7
        Fh = None; bump = False
        for r in MCrows:
            if r[0] <= hi and r[1] >= lo:
                Fh = r[2] if Fh is None else Fh.union(r[2]); bump = bump or r[4]
        assert Fh is not None
    for n in NAMES:
        Hc = (Rc * uc[n]).coeffs() + [arb(0)] * 3; Hh = (Rh * uh[n]).coeffs() + [arb(0)] * 3
        T = L * Hc[0] + I1 * Hc[1] + I2 * Hh[2]
        sup = ab(Hh[0])                                  # rho*v hull on box (abs sup of hull value)
        Fm = arb(Fh.mid()); Fr = arb(Fh.rad())
        res[n] += Fm * T; resmid[n] += Fm * T
        # rigorous abs bound int|rho v| <= L * sup|rho v| (hull, over box length |L|)
        res[n] += arb(0, float((Fr * ab(L) * sup).upper()))
        if bump: l1b[n] += ab(L) * sup
B.accumulate = acc_hook
if __name__ == "__main__":
    t0 = time.time(); NF = 7; G = [[arb(0)] * NF for _ in range(NF)]
    state["phase"] = "far"; B.far_region(G, 1)
    state["phase"] = "tip"; B.tip_region(G, 3)
    assert state["idx"] == len(MBrows)
    for n in NAMES: res[n] += arb(0, float((l1b[n] * BUMP).upper()))
    out = {n: [float(res[n].lower()), float(res[n].upper())] for n in NAMES}
    for n in NAMES: print(n, out[n], "rel rad %.2e" % (float(res[n].rad() / abs(res[n].mid()))))
    out['_mid_F_nonrigorous'] = {n: [float(resmid[n].lower()), float(resmid[n].upper())] for n in NAMES}
    json.dump(out, open('pair_f3.json', 'w'), indent=1); print(f"[{time.time()-t0:.0f}s]")
