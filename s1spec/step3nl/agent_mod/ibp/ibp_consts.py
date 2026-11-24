"""Global constants of the modulation lemma that are NOT changed by the integration by parts:
Gram G of (psi_i, f_j), H = (G^-1)^T, leading linear coefficients (Pi-projected Gram of (L*f1, L*f2), eps1), forcing coefficients,
E-matrix (feedback of x = (a,b,beta) into the constraint derivative).  Copied from ../lemma/lemma_consts.py (absolute paths, no output file)."""
import sys, json
import numpy as np
from flint import arb, arb_mat, ctx
ctx.prec = 150
sys.path.insert(0, '/home/claude/s1spec/step3nl/agent_mod')
import mod_forms as MF
LEM = '/home/claude/s1spec/step3nl/agent_mod/lemma/'
mag = lambda z: arb(abs(z).upper())
up = lambda z: float(z.upper())
I = json.load(open(LEM + 'lem_integrals.json')); LY = json.load(open('/home/claude/s1spec/step3/lyap_consts2.json')); LX = json.load(open(LEM + 'lyap_extra.json'))
P3 = json.load(open('/home/claude/s1spec/step3nl/psi3ip.json')); H2 = json.load(open('/home/claude/s1spec/step3nl/h2res.json'))
GR = MF.load('/home/claude/s1spec/step3nl/agent_mod/gram8_g1s4.json')
iv = lambda lo, hi: MF.ball(lo, hi)
SP = I['SP']
M = [[GR[6][6], GR[6][7]], [GR[7][6], GR[7][7]]]
nf3 = arb(LX['nf3_hi'])
g13 = LY['psif3'][0] + SP['1_F']; g23 = LY['psif3'][1] + SP['2_F']
g31 = iv(P3['psi3_f1'][0] - SP['y_1'], P3['psi3_f1'][1] + SP['y_1']); g32 = iv(P3['psi3_f2'][0] - SP['y_2'], P3['psi3_f2'][1] + SP['y_2'])
g33 = iv(P3['psi3_f3'][0] - SP['y_F'], P3['psi3_f3'][1] + SP['y_F'])
Gm = arb_mat([[M[0][0], M[0][1], arb(0, g13)], [M[1][0], M[1][1], arb(0, g23)], [g31, g32, g33]])
det_ = M[0][0] * M[1][1] - M[0][1] * M[1][0]; Mi_ = [[M[1][1] / det_, -M[0][1] / det_], [-M[1][0] / det_, M[0][0] / det_]]
cvec = [arb(0, g13), arb(0, g23)]; rvec = [g31, g32]
Mc = [sum((Mi_[a][b] * cvec[b] for b in range(2)), arb(0)) for a in range(2)]; rM = [sum((rvec[b] * Mi_[b][a] for b in range(2)), arb(0)) for a in range(2)]
s_ = g33 - sum((rvec[a] * Mc[a] for a in range(2)), arb(0)); assert s_.lower() > 0
Ginv = [[Mi_[a][b] + Mc[a] * rM[b] / s_ for b in range(2)] + [-Mc[a] / s_] for a in range(2)] + [[-rM[0] / s_, -rM[1] / s_, 1 / s_]]
Gi = arb_mat(Ginv)
chk = Gm * Gi
for i in range(3):
    for j in range(3): assert chk[i, j].contains(1 if i == j else 0), (i, j, chk[i, j])
H = [[Ginv[j][i] for j in range(3)] for i in range(3)]          # H = (G^-1)^T
det = M[0][0] * M[1][1] - M[0][1] * M[1][0]; Mi = [[M[1][1] / det, -M[0][1] / det], [-M[1][0] / det, M[0][0] / det]]
Gp = [[GR[a][b] - sum((GR[a][6 + i] * Mi[i][j] * GR[6 + j][b] for i in range(2) for j in range(2)), arb(0)) for b in range(2)] for a in range(2)]
eps1 = arb(LY['eps1'])
def qform(h, Gm_):
    s = sum((h[a] * h[b] * Gm_[a][b] for a in range(len(h)) for b in range(len(h))), arb(0)); return arb(max(s.upper(), 0)).sqrt()
def lam_max_sym2(S):
    t = S[0][0] + S[1][1]; d = S[0][0] * S[1][1] - S[0][1] * S[1][0]; disc = arb(max((t * t - 4 * d).upper(), 0))
    return arb((t.upper() + disc.sqrt().upper()) / 2)
lead = [qform([H[i][0], H[i][1]], Gp) + mag(H[i][2]) * eps1 for i in range(3)]
C2 = [[H[0][0], H[0][1]], [H[1][0], H[1][1]]]
S2 = [[sum((C2[a][k] * Gp[k][l] * C2[b][l] for k in range(2) for l in range(2)), arb(0)) for b in range(2)] for a in range(2)]
lead_vec = lam_max_sym2(S2).sqrt() + (H[0][2] * H[0][2] + H[1][2] * H[1][2]).sqrt() * eps1
Gf = [[M[0][0], M[0][1], arb(0, g13)], [M[1][0], M[1][1], arb(0, g23)], [arb(0, g13), arb(0, g23), nf3 * nf3]]
cF = [qform([H[i][j] for j in range(3)], Gf) for i in range(3)]
S3 = [[sum((H[a][k] * Gf[k][l] * H[b][l] for k in range(3) for l in range(3)), arb(0)) for b in range(2)] for a in range(2)]
cF_vec = lam_max_sym2(S3).sqrt()
# c_j: per-pairing leading constants  |<w~, L* f_j>| <= c_j A  (w~ orthogonal to f1,f2: projected norms for j=1,2; eps1 for j=3)
c_j = [arb(max(Gp[0][0].upper(), 0)).sqrt(), arb(max(Gp[1][1].upper(), 0)).sqrt(), eps1]
# ---- E matrix pieces (feedback), unchanged
h2n = arb(H2['h2_norm2_half'][1]).sqrt()
nfc = {1: arb(up(M[0][0].sqrt()) * (1 + 1e-12)), 2: arb(up(M[1][1].sqrt()) * (1 + 1e-12)), 3: nf3}
Ea = [[arb(I['e1_%d' % j]) for j in (1, 2, 3)], [arb(I['e2_%d' % j]) for j in (1, 2, 3)], [2 * h2n * nfc[j] for j in (1, 2, 3)]]
Ew1 = [[arb(I['Pq1_%d' % j]) for j in (1, 2, 3)], [arb(I['Pc1_%d' % j]) for j in (1, 2, 3)], [arb(0)] * 3]
Ew2 = [[arb(I['Pq0_%d' % j]) for j in (1, 2, 3)], [arb(I['Pc0_%d' % j]) for j in (1, 2, 3)], [arb(0)] * 3]
Es1 = [[arb(I['ISq_%d' % j]) for j in (1, 2, 3)], [arb(I['ISc_%d' % j]) for j in (1, 2, 3)], [arb(0)] * 3]
Es0 = [[arb(I['ISqk_%d' % j]) for j in (1, 2, 3)], [arb(I['ISck_%d' % j]) for j in (1, 2, 3)], [arb(0)] * 3]
Hab = lambda i, j: mag(H[i][j])
def Nmat(E): return [[sum((Hab(i, j) * E[k][j] for j in range(3)), arb(0)) for k in range(3)] for i in range(3)]
N0, N1, N2, N3, N4 = (Nmat(E) for E in (Ea, Ew1, Ew2, Es1, Es0))
# strip commutator coefficients kappa_{j,i} (coefficients of S_i = sup_strip |w^(i)|, i = 0..3): unchanged ([L,chi] w pairing)
KAP = [I[f'kappa_{j}'] for j in (1, 2, 3)]
IS = [I[f'IS_{j}'] for j in (1, 2, 3)]
if __name__ == "__main__":
    f = lambda z: float(z.mid())
    print("H mid:", [[f(H[i][j]) for j in range(3)] for i in range(3)])
    print("lead:", [up(x) for x in lead], up(lead_vec), " cF:", [up(x) for x in cF], up(cF_vec))
    print("c_j:", [up(x) for x in c_j]); print("kappa:", KAP)
    print("N0 (alpha)", [[up(x) for x in r] for r in N0]); print("N1 (W1)", [[up(x) for x in r] for r in N1])
    print("N2", [[up(x) for x in r] for r in N2]); print("N3", [[up(x) for x in r] for r in N3]); print("N4", [[up(x) for x in r] for r in N4])
