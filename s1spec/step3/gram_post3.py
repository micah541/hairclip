"""Rigorous Gram matrix via basis-invariant minors:
   G_ij = <A^{-1} f_j, f_i> = [y1^(j) m_{3,i} - y3^(j) m_{1,i}] / W13 + I_i^(j)
with W = Pluecker vector of the admissible plane of A, m_{r,i} = za_r I_b,i - zb_r I_a,i (both from run ZM, same frame),
y^(j), I^(j) the forced solution for f_j with zero data at sigma = 32 (run QS_j).  Generators: p, cos K, sin K, w*.
Constraints: f1 = p, f2 = c1 cosK + c2 sinK (c = tip midpoints: f2 = -sin(theta + const), const tiny), f3 = psi~3* = w*/rho."""
import json
import numpy as np
from flint import arb, ctx
ctx.prec = 128

def tip(path):
    for l in open(path):
        if l.startswith("TIP "):
            v = list(map(float, l.split()[1:]))
            return [arb((v[2 * i] + v[2 * i + 1]) / 2) + arb(0, 1) * ((v[2 * i + 1] - v[2 * i]) / 2 * (1 + 1e-15) + 1e-300)
                    for i in range(len(v) // 2)]
    raise RuntimeError(path)

ZM = tip("out3ZM.txt")
W = ZM[10:16]; W01, W13 = W[0], W[4]
m = [[ZM[16 + 4 * r + i] for i in range(4)] for r in range(4)]
Ctip, Stip, ws1 = ZM[4], ZM[5], ZM[7]
G4 = [[None] * 4 for _ in range(4)]
for j in range(4):
    Q = tip(f"out3QS{j}.txt")
    y, I = Q[10:14], Q[14:18]
    for i in range(4):
        G4[i][j] = (y[1] * m[3][i] - y[3] * m[1][i]) / W13 + I[i]
print("W13(tip) =", W13, "  W01(tip) =", W01)
print("G4 (generators p, cosK, sinK, w*):")
for i in range(4):
    print("  " + "  ".join(f"{float(G4[i][j].mid()):+.6e}(+-{float(G4[i][j].rad()):.1e})" for j in range(4)))
symok = all((G4[i][j] - G4[j][i]).contains(0) for i in range(4) for j in range(4))
print("symmetry consistent:", symok)
print("psi*'(tip) ~ w*_1(tip) =", ws1, "  (vs |w*| at tip ~", ZM[6], ")")
c1, c2 = arb(float(Ctip.mid())), arb(float(Stip.mid()))
R = [[arb(1), arb(0), arb(0), arb(0)], [arb(0), c1, c2, arb(0)], [arb(0), arb(0), arb(0), arb(1)]]
G3 = [[sum(R[a][i] * G4[i][j] * R[b][j] for i in range(4) for j in range(4)) for b in range(3)] for a in range(3)]
Gs = [[(G3[a][b] + G3[b][a]) / 2 for b in range(3)] for a in range(3)]
print("C_tip, S_tip =", Ctip, Stip)
print("G3 (symmetrized):")
for a in range(3):
    print("  ", [f"{float(Gs[a][b].mid()):+.6e}(+-{float(Gs[a][b].rad()):.1e})" for b in range(3)])
N = [[-Gs[a][b] for b in range(3)] for a in range(3)]
m1 = N[0][0]
m2 = N[0][0] * N[1][1] - N[0][1] * N[1][0]
m3 = (N[0][0] * (N[1][1] * N[2][2] - N[1][2] * N[2][1]) - N[0][1] * (N[1][0] * N[2][2] - N[1][2] * N[2][0])
      + N[0][2] * (N[1][0] * N[2][1] - N[1][1] * N[2][0]))
neg = bool(m1 > 0 and m2 > 0 and m3 > 0)
print("leading principal minors of -G3:", m1, m2, m3)
print("G3 NEGATIVE DEFINITE (rigorous):", neg)
Gn = np.array([[float(Gs[a][b].mid()) for b in range(3)] for a in range(3)])
D = np.sqrt(abs(np.diag(Gn)))
print("eig(G3 mid) =", np.linalg.eigvalsh(Gn), " normalized:", np.linalg.eigvalsh(Gn / np.outer(D, D)))
json.dump(dict(neg=neg, G3=[[float(Gs[a][b].mid()) for b in range(3)] for a in range(3)],
               G3rad=[[float(Gs[a][b].rad()) for b in range(3)] for a in range(3)],
               minors=[str(m1), str(m2), str(m3)], c=[float(c1), float(c2)]), open("gram_result.json", "w"), indent=1)
