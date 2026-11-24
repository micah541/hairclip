"""Modulation lemma constants (see MOD_LEMMA.md).  Unknowns x = (a, b, beta);  G^T-hat x = R,  G-hat = G + E,
 x = H R - H E^T x,  H = (G^-1)^T.   |x| <= t + N |x|  (componentwise),  N_ik = sum_j |H_ij| |E_kj|.
Usage: python3 lemma_consts.py [nlin_json]   (default nlin_E2_0.001_0.01_0.01.json = jet box eta = 1e-2)"""
import sys, json
import numpy as np
from flint import arb, arb_mat, ctx
ctx.prec = 150
sys.path.insert(0, '/home/claude/s1spec/step3nl/agent_mod')
import mod_forms as MF
NL = sys.argv[1] if len(sys.argv) > 1 else 'nlin_E2_0.001_0.01_0.01.json'
TAG = sys.argv[2] if len(sys.argv) > 2 else 'eta1e-2'
AMAX = 1e-3; W1REF = 1e-3; AREF = 1e-3; SREF = 1e3
mag = lambda z: arb(abs(z).upper())
up = lambda z: float(z.upper())
I = json.load(open('lem_integrals.json')); LY = json.load(open('/home/claude/s1spec/step3/lyap_consts2.json')); LX = json.load(open('lyap_extra.json'))
P3 = json.load(open('/home/claude/s1spec/step3nl/psi3ip.json')); H2 = json.load(open('/home/claude/s1spec/step3nl/h2res.json'))
GR = MF.load('/home/claude/s1spec/step3nl/agent_mod/gram8_g1s4.json')
iv = lambda lo, hi: MF.ball(lo, hi)
# ---------------------------------------------------------------- G (rows psi_i, cols f_j)
SP = I['SP']
M = [[GR[6][6], GR[6][7]], [GR[7][6], GR[7][7]]]
nf3 = arb(LX['nf3_hi'])
g13 = LY['psif3'][0] + SP['1_F']; g23 = LY['psif3'][1] + SP['2_F']      # |<chi psi_i, f3>| <= psif3 + strip difference
g31 = iv(P3['psi3_f1'][0] - SP['y_1'], P3['psi3_f1'][1] + SP['y_1']); g32 = iv(P3['psi3_f2'][0] - SP['y_2'], P3['psi3_f2'][1] + SP['y_2'])
g33 = iv(P3['psi3_f3'][0] - SP['y_F'], P3['psi3_f3'][1] + SP['y_F'])
Gm = arb_mat([[M[0][0], M[0][1], arb(0, g13)], [M[1][0], M[1][1], arb(0, g23)], [g31, g32, g33]])
# block inverse (avoids the dependency blow-up of ball Gaussian elimination):  G = [[M, c], [r, g33]]
det_ = M[0][0] * M[1][1] - M[0][1] * M[1][0]; Mi_ = [[M[1][1] / det_, -M[0][1] / det_], [-M[1][0] / det_, M[0][0] / det_]]
cvec = [arb(0, g13), arb(0, g23)]; rvec = [g31, g32]
Mc = [sum((Mi_[a][b] * cvec[b] for b in range(2)), arb(0)) for a in range(2)]; rM = [sum((rvec[b] * Mi_[b][a] for b in range(2)), arb(0)) for a in range(2)]
s_ = g33 - sum((rvec[a] * Mc[a] for a in range(2)), arb(0)); assert s_.lower() > 0
Ginv = [[Mi_[a][b] + Mc[a] * rM[b] / s_ for b in range(2)] + [-Mc[a] / s_] for a in range(2)] + [[-rM[0] / s_, -rM[1] / s_, 1 / s_]]
Gi = arb_mat(Ginv)
chk = Gm * Gi                                   # must contain the identity
for i in range(3):
    for j in range(3): assert chk[i, j].contains(1 if i == j else 0), (i, j, chk[i, j])
H = [[Ginv[j][i] for j in range(3)] for i in range(3)]
for i in range(3):
    for j in range(3): assert H[i][j].is_finite()
# ---------------------------------------------------------------- Pi2-projected Gram of (L*f1, L*f2)
det = M[0][0] * M[1][1] - M[0][1] * M[1][0]; Mi = [[M[1][1] / det, -M[0][1] / det], [-M[1][0] / det, M[0][0] / det]]
Gp = [[GR[a][b] - sum((GR[a][6 + i] * Mi[i][j] * GR[6 + j][b] for i in range(2) for j in range(2)), arb(0)) for b in range(2)] for a in range(2)]
eps1 = arb(LY['eps1'])
def qform(h, Gm_):
    s = sum((h[a] * h[b] * Gm_[a][b] for a in range(len(h)) for b in range(len(h))), arb(0)); return arb(max(s.upper(), 0)).sqrt()
def lam_max_sym2(S):
    t = S[0][0] + S[1][1]; d = S[0][0] * S[1][1] - S[0][1] * S[1][0]; disc = arb(max((t * t - 4 * d).upper(), 0))
    return arb((t.upper() + disc.sqrt().upper()) / 2)
lead = [qform([H[i][0], H[i][1]], Gp) + mag(H[i][2]) * eps1 for i in range(3)]
C = [[H[0][0], H[0][1]], [H[1][0], H[1][1]]]
S2 = [[sum((C[a][k] * Gp[k][l] * C[b][l] for k in range(2) for l in range(2)), arb(0)) for b in range(2)] for a in range(2)]
lead_vec = lam_max_sym2(S2).sqrt() + (H[0][2] * H[0][2] + H[1][2] * H[1][2]).sqrt() * eps1
Gf = [[M[0][0], M[0][1], arb(0, g13)], [M[1][0], M[1][1], arb(0, g23)], [arb(0, g13), arb(0, g23), nf3 * nf3]]
cF = [qform([H[i][j] for j in range(3)], Gf) for i in range(3)]
S3 = [[sum((H[a][k] * Gf[k][l] * H[b][l] for k in range(3) for l in range(3)), arb(0)) for b in range(2)] for a in range(2)]
cF_vec = lam_max_sym2(S3).sqrt()
# ---------------------------------------------------------------- D-term coefficients (nlin_E2)
NLd = json.load(open(NL))['sup']
def dco(reg):
    s = NLd[reg]
    return {"D4": s['(0, 0, 1)'], "W3": s['(0, 1, 0)'], "B": s['(1, 0, 0)'], "A": s['A0'], "W1": s['A1'], "S2B": s['(2, 0, 0)'], "S2W3": s['(1, 1, 0)'], "S2S2B": s['(3, 0, 0)']}
DC = {"tip": dco('x<=12'), "far": dco('12<x<=34.25')}
QN = ["A", "W1", "B", "W3", "D4", "S2B", "S2W3", "S2S2B", "Phi", "S0", "S1", "S2", "S3", "S4", "S2S2", "S2S3", "S2S2S2"]
nfj = {j: {"tip": I[f'nf_{j}']['x<=12'], "far": I[f'nf_{j}']['x>12']} for j in (1, 2, 3)}
nfj[1] = {r: v * 1.0 for r, v in nfj[1].items()}
tc = [{q: arb(0) for q in QN} for _ in range(3)]; tcv = {q: [arb(0), arb(0)] for q in QN}
Hab = lambda i, j: mag(H[i][j])
for i in range(3):
    tc[i]["A"] += lead[i]
    for j in (1, 2, 3):
        for reg in ("tip", "far"):
            for q, cf in DC[reg].items(): tc[i][q] += Hab(i, j - 1) * arb(nfj[j][reg] * (1 + 1e-12)) * arb(cf)
        # strip pieces: commutator kappa_{j,i'} S_i'  and  IS_j * far coefficients
        kap = I[f'kappa_{j}']
        for ip in range(4): tc[i][f"S{ip}"] += Hab(i, j - 1) * arb(kap[ip])
        far = DC["far"]; ISj = arb(I[f'IS_{j}'])
        for q, key in (("S4", "D4"), ("S3", "W3"), ("S2", "B"), ("S0", "A"), ("S1", "W1"), ("S2S2", "S2B"), ("S2S3", "S2W3"), ("S2S2S2", "S2S2B")):
            tc[i][q] += Hab(i, j - 1) * ISj * arb(far[key])
    tc[i]["Phi"] += cF[i]
# vector (a, b): leading via lead_vec, the rest by components
vec = {q: arb((tc[0][q] * tc[0][q] + tc[1][q] * tc[1][q]).upper()).sqrt() for q in QN}
vec["A"] = arb(((tc[0]["A"] - lead[0]) ** 2 + (tc[1]["A"] - lead[1]) ** 2).upper()).sqrt() + lead_vec
vec["Phi"] = cF_vec
# ---------------------------------------------------------------- E matrix pieces  N = |H| |E|^T
h2n = arb(H2['h2_norm2_half'][1]).sqrt()
nfc = {1: arb(up(M[0][0].sqrt()) * (1 + 1e-12)), 2: arb(up(M[1][1].sqrt()) * (1 + 1e-12)), 3: nf3}
Ea = [[arb(I['e1_%d' % j]) for j in (1, 2, 3)], [arb(I['e2_%d' % j]) for j in (1, 2, 3)], [2 * h2n * nfc[j] for j in (1, 2, 3)]]          # coefficient of |alpha|
Ew1 = [[arb(I['Pq1_%d' % j]) for j in (1, 2, 3)], [arb(I['Pc1_%d' % j]) for j in (1, 2, 3)], [arb(0)] * 3]                              # coefficient of W1
Ew2 = [[arb(I['Pq0_%d' % j]) for j in (1, 2, 3)], [arb(I['Pc0_%d' % j]) for j in (1, 2, 3)], [arb(0)] * 3]                              # |alpha| A
Es1 = [[arb(I['ISq_%d' % j]) for j in (1, 2, 3)], [arb(I['ISc_%d' % j]) for j in (1, 2, 3)], [arb(0)] * 3]                              # S1
Es0 = [[arb(I['ISqk_%d' % j]) for j in (1, 2, 3)], [arb(I['ISck_%d' % j]) for j in (1, 2, 3)], [arb(0)] * 3]                            # |alpha| S0
def Nmat(E): return [[sum((Hab(i, j) * E[k][j] for j in range(3)), arb(0)) for k in range(3)] for i in range(3)]
N0, N1, N2, N3, N4 = (Nmat(E) for E in (Ea, Ew1, Ew2, Es1, Es0))
Nref = [[N0[i][k] * AMAX + N1[i][k] * W1REF + N2[i][k] * AMAX * AREF + N3[i][k] * SREF + N4[i][k] * AMAX * SREF for k in range(3)] for i in range(3)]
Nf = np.array([[float(Nref[i][k].upper()) for k in range(3)] for i in range(3)])
ev, V = np.linalg.eig(Nf); v = np.abs(V[:, np.argmax(ev.real)].real); d = v + 1e-3 * v.max()
dA = [arb(float(x)) for x in d]
nu = max(float(sum((Nref[i][k] * dA[k] for k in range(3)), arb(0)).upper()) / float(d[i]) for i in range(3))
nu = nu * (1 + 1e-12); assert nu < 1, nu
tau = {q: max(float((tc[k][q] if q != "A" else tc[k][q]).upper()) / float(d[k]) for k in range(3)) for q in QN}
nd = lambda Nm: [sum((Nm[i][k] * dA[k] for k in range(3)), arb(0)) for i in range(3)]
n0, n1, n2, n3, n4 = (nd(Nm) for Nm in (N0, N1, N2, N3, N4))
fac = 1 / (1 - nu)
out = {"setting": {"alpha_max": AMAX, "jet_box_eta": TAG, "nlin_json": NL, "W1_ref": W1REF, "A_ref": AREF, "S_ref(strip sups)": SREF, "alpha_g_min": I['ag_min']},
       "G": [[[float(Gm[i, j].lower()), float(Gm[i, j].upper())] for j in range(3)] for i in range(3)],
       "H=(G^-1)^T": [[[float(H[i][j].lower()), float(H[i][j].upper())] for j in range(3)] for i in range(3)],
       "norm_G_inverse_inf": max(float(sum((mag(Gi[i, j]) for j in range(3)), arb(0)).upper()) for i in range(3)),
       "G33_lower": float(g33.lower()), "eps1": float(eps1), "Gpi_2x2_(Pi2 L*f)": [[[float(Gp[a][b].lower()), float(Gp[a][b].upper())] for b in range(2)] for a in range(2)],
       "lead_coeff_of_A": {"a": up(lead[0]), "b": up(lead[1]), "beta": up(lead[2]), "|(a,b)|": up(lead_vec)},
       "t_coeffs": {n: {q: up(tc[i][q]) for q in QN} for i, n in enumerate(("a", "b", "beta"))}, "t_coeffs_vec_ab": {q: up(vec[q]) for q in QN},
       "Neumann": {"nu_bar": nu, "weights_d": [float(x) for x in d], "tau_q(max_k t_kq/d_k)": tau, "1/(1-nu)": float(fac),
                   "n0(alpha)": [up(x) for x in n0], "n1(W1)": [up(x) for x in n1], "n2(alpha*A)": [up(x) for x in n2], "n3(S1)": [up(x) for x in n3], "n4(alpha*S0)": [up(x) for x in n4]},
       "DC": DC, "nf": nfj, "h2_norm": up(h2n)}
# feedback bounds g_i <= (|alpha| n0 + W1 n1 + |alpha| A n2 + S1 n3 + |alpha| S0 n4)_i * (sum_q tau_q q) / (1 - nu)
# At alpha = AMAX, W1 = W1REF, A = AREF, strip sups = 0: numeric size of the feedback coefficient
fb = [float((n0[i] * AMAX + n1[i] * W1REF + n2[i] * AMAX * AREF).upper()) * float(fac) for i in range(3)]
out["feedback_factor_at_ref(per unit tau_q q)"] = fb
json.dump(out, open(f'mod_lemma_{TAG}.json', 'w'), indent=1)
print(json.dumps({k: out[k] for k in ("G33_lower", "norm_G_inverse_inf", "lead_coeff_of_A")}, indent=1)); print("H =", out["H=(G^-1)^T"])
for n in ("a", "b", "beta"): print(n, {q: float(f"{v:.4g}") for q, v in out["t_coeffs"][n].items() if v > 1e-30})
print("vec", {q: float(f"{v:.4g}") for q, v in out["t_coeffs_vec_ab"].items() if v > 1e-30})
print("nu_bar", nu, "d", d, "n0", [up(x) for x in n0], "n1", [up(x) for x in n1], "fb", fb)
print("tau", {q: float(f"{v:.4g}") for q, v in tau.items() if v > 1e-30})
