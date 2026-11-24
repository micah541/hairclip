"""FLOAT feasibility of the COUPLED bootstrap  Lyap = E + mu E_I  (alpha handled by the comparison ODE), replacing the CFAR
placeholder of inner_closure_float.py by the intermediate energy E_I = ||v||_W^2 + epsI ||v''||_W^2, v = zeta w' (design of
BOOTSTRAP.md 'Far-field coupling: DECIDED DESIGN', P = 2).  NOT a proof.  Constants marked [P3] are proxies from the old P = 3
intermediate run (agent_int), [EST] are estimates, [DES] from design/coupling_design2.py (float); all others certified.
Worst case at fixed Lyap over the split lam = E/Lyap and the shapes (A,B,D4) inner, (V0,V2,V4) intermediate."""
import numpy as np, json, sys, os
from scipy.integrate import solve_ivp
from scipy.interpolate import RegularGridInterpolator
EPS = float(sys.argv[1]) if len(sys.argv) > 1 else 1e-5
ALPHA1 = float(sys.argv[2]) if len(sys.argv) > 2 else 1e-3
W1_0 = float(sys.argv[3]) if len(sys.argv) > 3 else 1e-6        # ||w~|| at hand-off
VI_0 = float(sys.argv[4]) if len(sys.argv) > 4 else 1e-6        # ||v||_W at hand-off
MU = float(os.environ.get("MU", 0.02)); EPSI = float(os.environ.get("EPSI", 3e-5))
c1 = 0.5804; PHX = 1.083; K2 = 2.078; K12 = 992.4; KR16 = 1 / np.exp(-(1 + 16.25**2) ** (2 / 3) / 4)   # 1/rho(16.25)
LAMT = -0.2561532; PI3 = 0.01607; G33 = 0.031105; EPS1 = 6.645e-5
a = dict(a00=55.6, a01=58.1, a02=30.8, a10=29.8, a11=38.2, a12=0.5, a20=5.9, a21=7.6, a22=6.2)
MOD_ALPHA = 6.19; MOD_E2 = 43.33;   # MOD_ALPHA: Q1 + beta-loop combined on {f1,f2,f3}^perp (agent_mod/bloop)
C3 = 12.104
CAB = float(os.environ.get('CAB', 7.57))   # |(a,b)| <= CAB A (agent_mod/lemma: leading 3.708, with nonlinear feedback at jet box 1e-4: 7.57)
cubj = json.load(open('cubic_E0_0.001_0.01_0.01.json'))["sup"]
cub = {"tip": cubj["x<=2"], "mid": cubj["2<x<=12"], "far": cubj["12<x<=31.25"]}
cubf = json.load(open('cubic_E0fine_0.001_0.01_0.01.json'))["sup"]       # finer regions (cubic_E0fine.py)
KREG = {"0<x<=1": 1.488, "1<x<=2": 2.078, "2<x<=4": 5.222, "4<x<=8": 56.92, "8<x<=12": 992.4, "12<x<=16.25": 3.0208e4, "16.25<x<=31.25": None}
# intermediate design (P = 2)
MU1, MU2 = 0.1552, 0.6555                 # L-int level 0 (certified, XA = 10.25); level-2 Gaarding -mu2 <= -0.6555 (agent_int/v2)
KW = 2.323e6; PSI1W = 1.305; C1W = 0.8477  # K_W [30,34.25], Psi1, c1W for P = 2, R = 60 (agent_sob/sob_constants_P2.json)
PHI34 = 9.8547e6                          # int_16.25^34.25 W^-1 (agent_coup)
SB, SA, IA = 3.5955e-5, 1.3372e-4, 9.5997e-12   # strip constants (agent_coup, rigorous)
AF, BF = 6.9555e-2, 3.9359e-2             # inner ramp flux (agent_coup, rigorous)
CF = (0.4421, 0.5259, 0.4235, 0.4250)     # level-2 ramp flux c3, c2, c1, c0 (agent_coup, rigorous)
DS = (6.841e-3, 1.5418e-2, 2.3648e-2, 4.5627e-2)  # level-2 strip d3..d0 (agent_coup)
W3125 = (1 + 31.25**2) ** -2.0
C3W = 0.1823                              # ||N0(u0) - a^2 Q2||_W <= C3W |a|^3 (agent_int/v2, P = 2)
TQ = dict(d001=1.85e-2, d010=1.5e-3, d100=1.1e-4, d200=1.5e-2, d110=0.98, d300=3.26, A0=1.7e-5, A1=7.1e-5)  # agent_int/v2 V3 (R = 60)
QPSI = 0.273 + 0.0159                     # ||zeta A1(u0)'||_W/|a| + ||zeta B1(u0)'||_W/|a| (agent_int/v2)
def Fn(al):
    al = abs(al); return 7.93e-7 * al + 2.769e-6 * al ** 2 + C3 * al ** 3      # C^4 profiles (agent_int/c4switch)
def FnW(al):
    al = abs(al); return 1.09e-7 * al + 6.5e-9 * al ** 2 + C3W * al ** 3       # W-residuals of psi3^, h2^ (agent_int/v2)
def rhs(al, A, B, D4, V0, V2, V4):
    al = abs(al)
    W1 = np.sqrt(A * B + c1 * A * A); W3 = np.sqrt(B * D4 + c1 * B * B)
    V1 = np.sqrt(V0 * V2 + C1W * V0 * V0); V3 = np.sqrt(V2 * V4 + C1W * V2 * V2)
    S = {}
    for reg, K in (("tip", K2), ("mid", K12)):
        S[reg] = (np.sqrt(K * (PHX * A * A + 2 * A * W1)), np.sqrt(K * (PHX * W1 * W1 + 2 * W1 * B)), np.sqrt(K * (PHX * B * B + 2 * B * W3)))
    sr = (np.sqrt(KR16 * (PHX * A * A + 2 * A * W1)), np.sqrt(KR16 * (PHX * W1 * W1 + 2 * W1 * B)), np.sqrt(KR16 * (PHX * B * B + 2 * B * W3)))
    sI1 = np.sqrt(KW * (PSI1W * V0 * V0 + 2 * V0 * V1)); sI2 = np.sqrt(KW * (PSI1W * V1 * V1 + 2 * V1 * V2))
    sI0 = sr[0] + np.sqrt(PHI34) * V0
    S["far"] = (np.maximum(sr[0], sI0), np.maximum(sr[1], sI1), np.maximum(sr[2], sI2))
    L2 = (A, W1, B)
    # ---- E0 (inner) ----
    lin = -(1 / 32 - PI3) * A * A - (1 / 9) * B * B
    if os.environ.get("ALIN") == "spec":      # alin_E0_spec.json: max over regions of the (w,w') spectral radius and cross terms
        alin = al * (59.3 * (A * A + W1 * W1) + 35.9 * A * B + 7.9 * W1 * B + 6.2 * B * B)
    else:
        alin = al * (a['a00'] * A * A + (a['a01'] + a['a10']) * A * W1 + (a['a02'] + a['a20']) * A * B + a['a11'] * W1 * W1
                     + (a['a12'] + a['a21']) * W1 * B + a['a22'] * B * B)
    mod = (MOD_ALPHA + float(os.environ.get('BLOOP', 0.0))) * al * A * A      # BLOOP: beta-loop alpha-form |beta_alpha| n3 / alpha
    forc = 5.72 * Fn(al) * A
    cubic = 0.0
    if os.environ.get("CUBFINE"):
        for reg, K in KREG.items():
            if K is None: s0, s1, s2 = sI0, sI1, sI2
            else: s0, s1, s2 = (np.sqrt(K * (PHX * A * A + 2 * A * W1)), np.sqrt(K * (PHX * W1 * W1 + 2 * W1 * B)), np.sqrt(K * (PHX * B * B + 2 * B * W3)))
            for key, C in cubf[reg].items():
                ab = eval(key); mono = [0] * ab[0] + [1] * ab[1] + [2] * ab[2]; mono.sort(reverse=True)
                keep, rest = mono[:2], mono[2:]
                val = C * L2[keep[0]] * L2[keep[1]]
                for d in rest: val = val * (s0, s1, s2)[d]
                cubic = cubic + 0.5 * val
    for reg in (() if os.environ.get("CUBFINE") else ("tip", "mid", "far")):
        s0, s1, s2 = S[reg]
        for key, C in cub[reg].items():
            ab = eval(key); mono = [0] * ab[0] + [1] * ab[1] + [2] * ab[2]; mono.sort(reverse=True)
            keep, rest = mono[:2], mono[2:]
            val = C * L2[keep[0]] * L2[keep[1]]
            for d in rest: val = val * (s0, s1, s2)[d]
            cubic = cubic + 0.5 * val
    strip = SB * V0 * V0 + SA * V0 * V0 + 2 * IA * KR16 * (A * A + 2 * A * W1)
    R0 = lin + alin + mod + forc + cubic + strip
    # ---- E2 (inner) ----
    supw2 = np.maximum(S["tip"][2], S["mid"][2]); supw2f = S["far"][2]
    ND = (6.9e-3 * D4 + 1.2e-2 * W3 + 4.2e-2 * B + 0.54 * A + 0.62 * W1 + 4.75 * supw2 * B + 0.136 * supw2 * W3 + 3.05 * supw2 ** 2 * B
          + 4.3e-3 * D4 + 6.9e-4 * W3 + 2.9e-5 * B + 0.014 * supw2f * B + 0.46 * supw2f * W3 + 3.06 * supw2f ** 2 * B + 8.6e-6 * A + 8.6e-6 * W1)
    wstrip = np.sqrt(6 * W3125) * (np.sqrt(KR16 * (PHX * A * A + 2 * A * W1)) + np.sqrt(PHI34) * V0)
    R2 = (DS[0] * V2 + DS[1] * V1 + DS[2] * V0 + DS[3] * wstrip) * (D4 + 2.2 * W3 + 1.2 * B) + (-0.5 * D4 * D4 + 14.34 * B * B + 16.3 * W1 * B + 10.7 * A * B + ND * (D4 + 2.2 * W3 + 1.2 * B)
          + MOD_E2 * A * A + Fn(al) * (D4 + 2.2 * W3 + 1.2 * B))
    # ---- intermediate level 0 ----
    supv1 = sI2                                     # sup |w''| on the far region (= sup |v'|)
    I0 = (-MU1 * V0 * V0 - 0.5 * V2 * V2 + AF * W1 * W1 + BF * B * B
          + TQ['d001'] * V2 * V2 + (TQ['d010'] * V2 + TQ['d100'] * V1 + TQ['A1'] * V0) * (V1 + 1.67 * V0)
          + (TQ['d200'] * supv1 * V1 + TQ['d110'] * supv1 * V2 + TQ['d300'] * supv1 ** 2 * V1) * (V1 + 1.67 * V0)
          + CAB * A * al * QPSI * V0 + CAB * A * (38.5 + 0.66) * V0 * V0 + CAB * A * 2.5 * W1 * W1 + FnW(al) * (V1 + 1.67 * V0))
    # ---- intermediate level 2 ----
    D2v = V4 + 2.61 * V3 + 1.7 * V2
    supv2 = np.sqrt(KW * (PSI1W * V2 * V2 + 2 * V2 * V3))          # sup |w'''| far
    I2 = (-0.5 * V4 * V4 - MU2 * V2 * V2 + (CF[0] * D4 + CF[1] * W3 + CF[2] * B + CF[3] * W1) * D2v
          + TQ['d001'] * V4 * V4 + (TQ['d010'] * V3 + TQ['d110'] * supv1 * V3 + 2 * TQ['d200'] * supv1 * V3 + TQ['d110'] * supv2 * V2
          + TQ['d300'] * supv1 ** 2 * V2) * D2v + 1.0 * abs(al) ** 3 * D2v
          + CAB * A * (38.5 + 0.66) * V2 * V2 + 62 * A * V2 * V3 + CAB * A * al * 0.3 * D2v)
    tot = R0 + EPS * R2 + MU * (I0 + EPSI * I2)
    return tot, dict(lin=lin, alin=alin, mod=mod, cubic=cubic, strip=strip, E2=EPS * R2, I0=MU * I0, I2=MU * EPSI * I2)
TH = np.linspace(1e-3, np.pi / 2 - 1e-3, 13)
def worst(al, Ly):
    """max over lam, shapes of (1/2) dLyap/dtau  (vectorised grid)."""
    lam = np.linspace(0.02, 0.98, 9)[:, None, None, None, None, None]
    th = TH[None, :, None, None, None, None]; thI = TH[None, None, None, :, None, None]
    E = lam * Ly; EI = (1 - lam) * Ly / MU
    A = np.sqrt(E) * np.cos(th); B = np.sqrt(E / EPS) * np.sin(th)
    V0 = np.sqrt(EI) * np.cos(thI); V2 = np.sqrt(EI / EPSI) * np.sin(thI)
    d4 = np.concatenate([[0.0], np.logspace(-1, 3, 9)])[None, None, :, None, None, None]
    v4 = np.concatenate([[0.0], np.logspace(-1, 3, 9)])[None, None, None, None, :, None]
    D4 = B * d4; V4 = V2 * v4
    tot, parts = rhs(al, A, B, D4, V0, V2, V4)
    i = np.unravel_index(np.argmax(tot), tot.shape)
    shp = tot.shape
    return tot[i], {k: float(np.broadcast_to(v, shp)[i]) for k, v in parts.items()}, i
if __name__ == "__main__":
    print(f"eps = {EPS}, epsI = {EPSI}, mu = {MU}, alpha1 = {ALPHA1}, ||w~||_1 = {W1_0}, ||v||_1 = {VI_0}")
    for Ly in (1e-14, 1e-12, 1e-10, 1e-9):
        for al in (ALPHA1, 0.0):
            b, pb, i = worst(al, Ly)
            print(f"  Lyap={Ly:.0e} alpha={al:.1e}: rate {b/Ly:+.4f}  [" + " ".join(f"{k} {v/Ly:+.3g}" for k, v in pb.items()) + f"] idx {i}")
    LE = np.linspace(-18, -8, 21); AL = np.linspace(0, ALPHA1 * 1.05, 6)
    TAB = np.array([[worst(al, 10 ** le)[0] / 10 ** le for le in LE] for al in AL])
    RT = RegularGridInterpolator((AL, LE), TAB, bounds_error=False, fill_value=None)
    def f(t, y):
        al, Ly = y
        b = RT([[abs(al), min(max(np.log10(max(Ly, 1e-30)), -18), -8)]])[0] * Ly
        beta = (0.00244 + 63.5 * abs(al)) * np.sqrt(max(Ly, 0)) + 0.626 * Fn(al)      # agent_mod/lemma (linear + forcing parts)
        return [LAMT * al + np.sign(al) * beta, 2 * b]
    Ly1 = W1_0 ** 2 * (1 + EPS * 9) + MU * VI_0 ** 2 * (1 + EPSI * 9)
    sol = solve_ivp(f, (0, 60), [ALPHA1, Ly1], max_step=1.0, rtol=1e-3)
    print("  max sqrt(Lyap):", np.sqrt(sol.y[1].max()), " final:", np.sqrt(sol.y[1][-1]), " (start", np.sqrt(Ly1), ")")
