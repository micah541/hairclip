import json, numpy as np
from flint import arb, ctx
ctx.prec = 128
GEN = ["p", "C", "S", "w"]
NA = {"p": 0, "C": 2, "S": 2, "w": 4}
def tipraw(path):
    for l in open(path):
        if l.startswith("TIP "):
            v = list(map(float, l.split()[1:])); return [(v[2*i], v[2*i+1]) for i in range(len(v)//2)]
    raise RuntimeError(path)
def A(lo, hi): return arb((lo+hi)/2) + arb(0, 1)*((hi-lo)/2*(1+1e-15) + 1e-300)
def isect(a, b): return (max(a[0], b[0]), min(a[1], b[1]))
def naux(gs):
    s = set()
    if "C" in gs or "S" in gs: s |= {"C", "S"}
    n = (2 if ("C" in gs or "S" in gs) else 0) + (4 if "w" in gs else 0)
    return n
W13, m1, m3, Wb = {}, {}, {}, {}
for gi in GEN:
    t = tipraw(f"out4Z{gi}.txt"); na = naux([gi]); o = 4 + na
    W = t[o:o+6]; mm = t[o+6:o+10]
    W13[gi] = A(*W[4]); m1[gi] = A(*mm[1]); m3[gi] = A(*mm[3])
y = {}; I = {}
for gj in GEN:
    enc = None
    for gi in GEN:
        t = tipraw(f"out4Q{gj}{gi}.txt"); na = naux([gj, gi]); o = 4 + na
        yy = t[o:o+4]; I[(gj, gi)] = A(*t[o+4])
        enc = yy if enc is None else [isect(a, b) for a, b in zip(enc, yy)]
    assert all(a <= b for a, b in enc)
    y[gj] = [A(*e) for e in enc]
G4 = [[(y[gj][1]*m3[gi] - y[gj][3]*m1[gi])/W13[gi] + I[(gj, gi)] for gj in GEN] for gi in GEN]
for i in range(4):
    print("  " + "  ".join(f"{float(G4[i][j].mid()):+.6e}(+-{float(G4[i][j].rad()):.1e})" for j in range(4)))
print("symmetry consistent:", all((G4[i][j]-G4[j][i]).contains(0) for i in range(4) for j in range(4)))
Zc = tipraw("out4ZC.txt"); Ct, St = A(*Zc[4]), A(*Zc[5])
c1, c2 = arb(float(Ct.mid())), arb(float(St.mid()))
R = [[arb(1), arb(0), arb(0), arb(0)], [arb(0), c1, c2, arb(0)], [arb(0), arb(0), arb(0), arb(1)]]
G3 = [[sum(R[a][i]*G4[i][j]*R[b][j] for i in range(4) for j in range(4)) for b in range(3)] for a in range(3)]
Gs = [[(G3[a][b] + G3[b][a])/2 for b in range(3)] for a in range(3)]
print("C_tip, S_tip:", Ct, St)
for a in range(3): print("  G3:", [f"{float(Gs[a][b].mid()):+.6e}(+-{float(Gs[a][b].rad()):.1e})" for b in range(3)])
N = [[-Gs[a][b] for b in range(3)] for a in range(3)]
mi1 = N[0][0]; mi2 = N[0][0]*N[1][1] - N[0][1]*N[1][0]
mi3 = (N[0][0]*(N[1][1]*N[2][2]-N[1][2]*N[2][1]) - N[0][1]*(N[1][0]*N[2][2]-N[1][2]*N[2][0]) + N[0][2]*(N[1][0]*N[2][1]-N[1][1]*N[2][0]))
print("minors of -G3:", mi1, mi2, mi3, " NEGATIVE DEFINITE:", bool(mi1 > 0 and mi2 > 0 and mi3 > 0))
Gn = np.array([[float(Gs[a][b].mid()) for b in range(3)] for a in range(3)]); D = np.sqrt(abs(np.diag(Gn)))
print("eig mid:", np.linalg.eigvalsh(Gn), " normalized:", np.linalg.eigvalsh(Gn/np.outer(D, D)))
