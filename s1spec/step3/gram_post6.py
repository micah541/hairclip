"""Rigorous Gram matrix from the hsrun3 doubleton outputs (mean-value / affine evaluation).
G[i][j] = <A^{-1} f_j, f_i>_rho (half line) = l_i . (y_j, I_ji),  l_i = (0, m3_i/W13, 0, -m1_i/W13, 1)."""
import sys, json
import numpy as np
from flint import arb, ctx
sys.path.insert(0, '/home/claude/s1spec/step3')
from gen_hs3 import parse3, aux_for, QSHIFT
ctx.prec = 128
GEN = ["p", "C", "S", "w"]
PREFIX = sys.argv[1] if len(sys.argv) > 1 else "out6"

def A(lo, hi):
    assert lo <= hi, (lo, hi)
    return arb((lo + hi) / 2) + arb(0, 1) * ((hi - lo) / 2 * (1 + 1e-15) + 1e-300)

class DSet:
    def __init__(self, path):
        D = parse3(open(path).read()); assert D["ok"], path
        self.x = [A(*v) for v in D["TX"]]
        self.C = [[A(*v) for v in row] for row in D["TC"]]
        self.B = [[A(*v) for v in row] for row in D["TB"]]
        self.r0 = [A(*v) for v in D["TR0"]]
        self.r = [A(*v) for v in D["TR"]]
        self.hull = [A(*v) for v in D["TIP"]]
        self.n = len(self.x)
    def affine(self, coef):
        """enclosure of sum_i coef[i] * z_i over the set, coef: dict row -> arb (may be intervals)"""
        val = sum((c * self.x[i] for i, c in coef.items()), arb(0))
        for k in range(self.n):
            gC = sum((c * self.C[i][k] for i, c in coef.items()), arb(0))
            gB = sum((c * self.B[i][k] for i, c in coef.items()), arb(0))
            val += gC * self.r0[k] + gB * self.r[k]
        return val
    def comp(self, i):
        return self.affine({i: arb(1)})

def ratio(S, im, iw):
    """enclosure of z_im / z_iw over the set by the mean-value form"""
    Hm, Hw = S.hull[im], S.hull[iw]
    gm, gw = 1 / Hw, -Hm / (Hw * Hw)
    centre = S.x[im] / S.x[iw]
    val = centre
    for k in range(S.n):
        val += (gm * S.C[im][k] + gw * S.C[iw][k]) * S.r0[k]
        val += (gm * S.B[im][k] + gw * S.B[iw][k]) * S.r[k]
    return val

def na(gs): return len(aux_for(gs))

def main():
    global G4, R
    rho1, rho3, W13 = {}, {}, {}
    for gi in GEN:
        S = DSet(f"{PREFIX}Z{gi}.txt"); o = 5 + na([gi])
        iW13, im1, im3 = o + 4, o + 7, o + 9
        rho1[gi] = ratio(S, im1, iW13); rho3[gi] = ratio(S, im3, iW13); W13[gi] = S.comp(iW13)
        print(f"Z{gi}: W13={S.hull[iW13]}  m1/W13={rho1[gi]}  m3/W13={rho3[gi]}")
    G4 = [[None] * 4 for _ in range(4)]
    for a, gi in enumerate(GEN):
        for b, gj in enumerate(GEN):
            S = DSet(f"{PREFIX}Q{gj}{gi}.txt"); o = 5 + na([gj, gi])
            coef = {o + 1: rho3[gi], o + 3: -rho1[gi], o + 4: arb(1)}
            G4[a][b] = S.affine(coef) * arb(2) ** QSHIFT
    for a in range(4):
        print("  " + "  ".join(f"{float(G4[a][b].mid()):+.7e}(+-{float(G4[a][b].rad()):.1e})" for b in range(4)))
    print("symmetry consistent:", all((G4[a][b] - G4[b][a]).contains(0) for a in range(4) for b in range(4)))
    Sc = DSet(f"{PREFIX}ZC.txt"); Ct, St = Sc.hull[5], Sc.hull[6]
    c1, c2 = arb(float(Ct.mid())), arb(float(St.mid()))
    json.dump([float(c1), float(c2)], open('f2_coef.json', 'w'))
    R = [[arb(1), arb(0), arb(0), arb(0)], [arb(0), c1, c2, arb(0)], [arb(0), arb(0), arb(0), arb(1)]]
    G3 = [[sum(R[a][i] * G4[i][j] * R[b][j] for i in range(4) for j in range(4)) for b in range(3)] for a in range(3)]
    Gs = [[(G3[a][b] + G3[b][a]) / 2 for b in range(3)] for a in range(3)]
    # both triangles are valid enclosures of the same symmetric matrix: intersect
    for a in range(3):
        for b in range(3):
            lo = max(float(G3[a][b].lower()), float(G3[b][a].lower())); hi = min(float(G3[a][b].upper()), float(G3[b][a].upper()))
    print("C_tip, S_tip:", Ct, St, " f2 = c1 cos K + c2 sin K with c =", float(c1), float(c2))
    for a in range(3): print("  G3:", [f"{float(Gs[a][b].mid()):+.7e}(+-{float(Gs[a][b].rad()):.1e})" for b in range(3)])
    N = [[-Gs[a][b] for b in range(3)] for a in range(3)]
    mi1 = N[0][0]; mi2 = N[0][0] * N[1][1] - N[0][1] * N[1][0]
    mi3 = (N[0][0] * (N[1][1] * N[2][2] - N[1][2] * N[2][1]) - N[0][1] * (N[1][0] * N[2][2] - N[1][2] * N[2][0])
           + N[0][2] * (N[1][0] * N[2][1] - N[1][1] * N[2][0]))
    print("leading minors of -G3:", mi1, mi2, mi3)
    print("G3 NEGATIVE DEFINITE (Sylvester):", bool(mi1 > 0 and mi2 > 0 and mi3 > 0))
    Gn = np.array([[float(Gs[a][b].mid()) for b in range(3)] for a in range(3)]); D = np.sqrt(abs(np.diag(Gn)))
    print("eig mid:", np.linalg.eigvalsh(Gn), " normalized:", np.linalg.eigvalsh(Gn / np.outer(D, D)))
    
    # ---------------------------------------------------------------- refined test: symmetric intersection + Schur complement
    def isect(a, b):
        lo = max(a.lower(), b.lower()); hi = min(a.upper(), b.upper())
        assert lo <= hi, "empty intersection: inconsistent enclosures"
        return (lo + hi) / 2 + arb(0, 1) * ((hi - lo) / 2)
    Gy = [[isect(G4[a][b], G4[b][a]) for b in range(4)] for a in range(4)]
    print("\nG4 (symmetric intersection):")
    for a in range(4):
        print("  " + "  ".join(f"{float(Gy[a][b].mid()):+.7e}(+-{float(Gy[a][b].rad()):.1e})" for b in range(4)))
    H3 = [[sum(R[a][i] * Gy[i][j] * R[b][j] for i in range(4) for j in range(4)) for b in range(3)] for a in range(3)]
    for a in range(3): print("  G3:", [f"{float(H3[a][b].mid()):+.7e}(+-{float(H3[a][b].rad()):.1e})" for b in range(3)])
    g11, g12, g22 = H3[0][0], (H3[0][1] + H3[1][0]) / 2, H3[1][1]
    g13, g23, g33 = (H3[0][2] + H3[2][0]) / 2, (H3[1][2] + H3[2][1]) / 2, H3[2][2]
    det2 = g11 * g22 - g12 * g12
    S = g33 - (g22 * g13 * g13 - 2 * g12 * g13 * g23 + g11 * g23 * g23) / det2
    print("G2: g11 =", g11, " det =", det2, "  => G2 negative definite:", bool(g11 < 0 and det2 > 0))
    print("Schur complement S = g33 - g^T G2^{-1} g =", S, "  => S < 0:", bool(S < 0))
    print("G3 NEGATIVE DEFINITE (Schur test):", bool(g11 < 0 and det2 > 0 and S < 0))
    print("S: mid", float(S.mid()), "rad", float(S.rad()), "upper", float(S.upper()), "lower", float(S.lower()))
    

if __name__ == '__main__':
    main()
