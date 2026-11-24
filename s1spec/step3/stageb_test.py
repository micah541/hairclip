import os, subprocess, gen_hs as G
from flint import arb
res, T = G.tail_data()
eta = res['eta']
prof, X0, E0, R0 = G.common_init()
var = ["X", "E"] + [f"W{i}{j}" for (i, j) in G.PAIRS]
full = ["p","q","k","m"] + var
fun_full = G.PROF + ["-1", "2*X*E^4/3"] + G.pluecker(var[2:])
funB = G.stageB(fun_full)[4:]
def pl(i, j):
    return [T[a, i]*T[b, j] - T[b, i]*T[a, j] for (a, b) in G.PAIRS]
import math
w23 = pl(2, 3)
scale = 2.0 ** -math.floor(math.log2(max(abs(float(w.mid())) for w in w23)))
x, box, cols, r0 = [], [], [], []
for v in (X0, E0): m, b = G.mid_box(v); x.append((m, m)); box.append(G.interval_pair(b))
extra = [arb(0)]*6
for colv, e in [(pl(0, 3), eta), (pl(1, 3), eta), (pl(0, 2), eta), (pl(1, 2), eta), (pl(0, 1), eta*eta)]:
    cvec = [0.0]*8
    for ij in range(6):
        m, b = G.mid_box(colv[ij]*scale); cvec[2+ij] = m; extra[ij] += abs(b)*e
    cols.append(cvec); r0.append((-float(e.upper())*(1+1e-12), float(e.upper())*(1+1e-12)))
for ij in range(6):
    m, b = G.mid_box(w23[ij]*scale); x.append((m, m)); box.append(G.interval_pair(b + arb(0, 1)*extra[ij]))
h = float(os.environ.get('HSTEP', '0.05'))
with open('sb.txt', 'w') as f:
    pars = G.params()
    f.write(f"8 {len(pars)}\n")
    for nm, lo, hi in pars: f.write(f"{nm} {lo!r} {hi!r}\n")
    vf = G.vf_string(var, funB)
    f.write(vf + "\n" + vf + "\n24.0\n")
    ch = G.chunks()
    f.write(f"{len(ch)} 4 kk mm pp dq\n")
    for ln, hh, pv in ch: f.write(f"{ln!r} {h!r} " + " ".join(f"{a!r} {b!r}" for a, b in pv) + "\n")
    for lo, hi in x: f.write(f"{lo!r} {hi!r}\n")
    f.write(f"{len(cols)}\n")
    for i in range(8): f.write(" ".join(repr(float(c[i])) for c in cols) + "\n")
    for lo, hi in r0: f.write(f"{lo!r} {hi!r}\n")
    for lo, hi in box: f.write(f"{lo!r} {hi!r}\n")
