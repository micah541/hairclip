"""Hull-box signed inner products <v, f3> (f3 = chi F - a_c eta1 - b_c eta3, F = ws0 hull of outA5.txt) for
v in (L*f1, L*f2, gA, gB, chi psi3, chi h2, f1, f2).  First pass: plain hulls (as psi3ip.py)."""
import sys, json, time
from flint import arb
sys.argv_saved = sys.argv; sys.argv = sys.argv[:1]
sys.path.insert(0, '/home/claude/s1spec/step3nl/agent_mod'); sys.path.insert(0, '/home/claude/s1spec/step3nl')
import mod_gram as MG
HR = MG.HR
sys.argv = sys.argv_saved
import gen_hs3 as G3
BUMP = (0.8324 * 9.0923e-8 + 0.4245 * 2.5707e-7) * 1.01
ab = lambda z: arb(abs(z).upper())
def run():
    t0 = time.time(); names = ["Lf1", "Lf2", "gA", "gB", "psi3", "h2", "f1", "f2", "F"]
    acc = {n: arb(0) for n in names}; l1 = {n: arb(0) for n in names}
    T0 = G3.T0
    lines = [l for l in open('/home/claude/s1spec/step3/outA5.txt') if l[:2] in ("MB", "MC")]
    for b, line in zip(MG.boxes(), lines):
        v = list(map(float, line.split()[1:])); tag = line[:2]
        h = [MG.A(v[2 + 2 * i], v[3 + 2 * i]) for i in range((len(v) - 2) // 2)]
        ws0 = h[5] if tag == "MB" else h[9]
        chi = b['chi']; Xv = b['X']; xc = float(Xv.mid())
        y = HR.jets_any(HR.F, Xv, xc)[:1]; hh = HR.jets_any(HR.Hf, Xv, xc)[:1]
        vals = dict(b['u']); vals["psi3"] = chi * y[0]; vals["h2"] = chi * hh[0]; vals["F"] = chi * ws0
        f3 = chi * ws0
        bump = (tag == "MC" and b['ta'] >= T0 - G3.LB - 1e-12)
        for n in names:
            acc[n] += b['w'] * vals[n] * f3
            if bump: l1[n] += b['w'] * ab(vals[n])
    for n in names:
        acc[n] = acc[n] + arb(0, float((l1[n] * BUMP).upper()))
    out = {n: [float(acc[n].lower()), float(acc[n].upper())] for n in names}
    for n in names: print(n, out[n], " mid %.6g rad %.2e" % (float(acc[n].mid()), float(acc[n].rad())))
    json.dump(out, open('hull_f3_pairings.json', 'w'), indent=1); print(f"[{time.time()-t0:.0f}s]")
if __name__ == "__main__": run()
