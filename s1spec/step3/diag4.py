import numpy as np
from gram_post4 import tipraw, naux, GEN
import gram_float as GF
sol, up = GF.run(); C, S, ws, za, zb, Ia, Ib, yy, II = up(sol.y[:, -1])
gi_idx = {g: i for i, g in enumerate(GEN)}
print("ZM runs (W13, m1, m3)")
for gi in GEN:
    t = tipraw(f"out4Z{gi}.txt"); o = 4 + naux([gi]); W = t[o:o+6]; mm = t[o+6:o+10]
    i = gi_idx[gi]
    m = [za[r]*Ib[i] - zb[r]*Ia[i] for r in range(4)]
    W13f = za[1]*zb[3] - zb[1]*za[3]
    for nm, enc, fv in [("W13", W[4], W13f), ("m1", mm[1], m[1]), ("m3", mm[3], m[3])]:
        mid, rad = (enc[0]+enc[1])/2, (enc[1]-enc[0])/2
        print(f"  {gi} {nm}: mid={mid:+.5e} rad={rad:.2e} (rel {rad/abs(mid):.1e})   float*scale?={fv:+.5e} ratio={mid/fv:.6f}")
print("Q runs (y1, y3, I)")
for gj in GEN:
    for gi in GEN:
        t = tipraw(f"out4Q{gj}{gi}.txt"); o = 4 + naux([gj, gi]); yv = t[o:o+4]; Iv = t[o+4]
        j, i = gi_idx[gj], gi_idx[gi]
        s = []
        for nm, enc, fv in [("y1", yv[1], yy[j, 1]), ("y3", yv[3], yy[j, 3]), ("I", Iv, II[j, i])]:
            mid, rad = (enc[0]+enc[1])/2, (enc[1]-enc[0])/2
            s.append(f"{nm}={mid:+.4e}±{rad:.1e}[f {fv:+.4e}]")
        print(f"  Q{gj}{gi}: " + "  ".join(s))
