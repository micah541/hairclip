"""Rigorous assembly of the Gram matrix G = [<A^{-1} f_j, f_i>] from the CAPD tip enclosures.
A = -H_sym - omega (even class, half line).  chi_j = b1 Z1 + b2 Z2 + Q_j with the even BCs u'(tip) = y4(tip) = 0."""
import sys
from flint import arb, arb_mat, ctx
ctx.prec = 128

def tip(path):
    for l in open(path):
        if l.startswith("TIP "):
            v = list(map(float, l.split()[1:]))
            return [arb((v[2 * i] + v[2 * i + 1]) / 2) + arb(0, 1) * ((v[2 * i + 1] - v[2 * i]) / 2 * (1 + 1e-15) + 1e-300)
                    for i in range(len(v) // 2)]
    raise RuntimeError("no TIP in " + path)

Z = tip("out2Z.txt")
za, zb = Z[14:18], Z[18:22]
Ia, Ib = Z[22:27], Z[27:32]
Ctip, Stip = Z[4], Z[5]
wa1, wb1 = Z[7], Z[11]
MZ = arb_mat([[za[1], zb[1]], [za[3], zb[3]]])
MZi = MZ.inv()
G5 = [[None] * 5 for _ in range(5)]
for j in range(5):
    Q = tip(f"out2Q{j + 1}.txt")
    y, IQ = Q[14:18], Q[18:23]
    rhs = arb_mat([[-y[1]], [-y[3]]])
    b = MZi * rhs
    for i in range(5):
        G5[i][j] = b[0, 0] * Ia[i] + b[1, 0] * Ib[i] + IQ[i]
print("G5 (generators p, cosK, sinK, wa/rho, wb/rho):")
for i in range(5):
    print("  " + "  ".join(f"{float(G5[i][j].mid()):+.5e}(+-{float(G5[i][j].rad()):.1e})" for j in range(5)))
asym = max(float(abs(G5[i][j] - G5[j][i]).upper()) / max(1e-300, float(abs(G5[i][j]).upper())) for i in range(5) for j in range(5))
symok = all((G5[i][j] - G5[j][i]).contains(0) for i in range(5) for j in range(5))
print(f"symmetry consistent (G5_ij - G5_ji contains 0 for all i,j): {symok};  max relative |G_ij - G_ji| bound {asym:.2e}")
# constraint functions: f1 = p, f2 = c1 cosK + c2 sinK (c = midpoints of (C, S) at the tip, i.e. f2 = -sin(theta) up to
# a constant phase of size |K_tip enclosure|), f3 = alpha wa/rho + beta wb/rho with psi*'(tip) ~ 0
c1, c2 = arb(float(Ctip.mid())), arb(float(Stip.mid()))
al, be = arb(float(wb1.mid())), arb(-float(wa1.mid()))
R = [[arb(1), 0, 0, 0, 0], [0, c1, c2, 0, 0], [0, 0, 0, al, be]]
G3 = [[sum(R[a][i] * G5[i][j] * R[b][j] for i in range(5) for j in range(5) if R[a][i] != 0 and R[b][j] != 0)
       for b in range(3)] for a in range(3)]
Gs = [[(G3[a][b] + G3[b][a]) / 2 for b in range(3)] for a in range(3)]
print("C_tip, S_tip =", Ctip, Stip, " (cos/sin of total turning K_tip)")
print("residual BC of f3 at tip: alpha*wa1 + beta*wb1 in", al * wa1 + be * wb1, " vs |alpha wa1| ~", abs(al * wa1))
print("G3 (sym):")
for a in range(3): print("  ", [f"{float(Gs[a][b].mid()):+.6e}(+-{float(Gs[a][b].rad()):.1e})" for b in range(3)])
# negative definiteness of G3 <=> leading principal minors of -G3 positive
N = [[-Gs[a][b] for b in range(3)] for a in range(3)]
m1 = N[0][0]
m2 = N[0][0] * N[1][1] - N[0][1] * N[1][0]
m3 = (N[0][0] * (N[1][1] * N[2][2] - N[1][2] * N[2][1]) - N[0][1] * (N[1][0] * N[2][2] - N[1][2] * N[2][0])
      + N[0][2] * (N[1][0] * N[2][1] - N[1][1] * N[2][0]))
print("leading minors of -G3:", m1, m2, m3)
print("G3 NEGATIVE DEFINITE (rigorous):", bool(m1 > 0 and m2 > 0 and m3 > 0))
# normalized smallest eigenvalue (diagnostic)
import numpy as np
Gn = np.array([[float(Gs[a][b].mid()) for b in range(3)] for a in range(3)])
D = np.sqrt(abs(np.diag(Gn)))
print("eig(G3) =", np.linalg.eigvalsh(Gn), " normalized:", np.linalg.eigvalsh(Gn / np.outer(D, D)))
