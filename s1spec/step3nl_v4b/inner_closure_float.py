"""FLOAT feasibility check of the inner bootstrap (E = ||w~||^2 + eps ||w~''||^2, half-line rho norms), using the certified
constants (alin_E0 1e-3, nlin_E2, garding_cert, dmod_norms, n3_forcing, cubic_E0, lyap_consts2).  NOT a proof.
Worst case over the 'shape' (A = ||w~||, B = ||w~''||, D4 = ||w~''''||) at fixed E; comparison ODE in tau with
alpha' = lamt alpha + beta, |beta| <= (eps1 A + 0.63 Fn)/G33-type bound.  The far-region pointwise bounds (x >= 12) are
modelled as sup|w^(j)| <= CFAR sqrt(E) (to be supplied by the intermediate energy)."""
import numpy as np, json, sys
from scipy.integrate import solve_ivp
EPS = float(sys.argv[1]) if len(sys.argv) > 1 else 1e-5
ALPHA1 = float(sys.argv[2]) if len(sys.argv) > 2 else 1e-3
W1_0 = float(sys.argv[3]) if len(sys.argv) > 3 else 1e-6     # ||w~|| at hand-off
CFAR = float(sys.argv[4]) if len(sys.argv) > 4 else 300.0
c1 = 0.5804; PHX = 1.083; K2 = 2.078; K12 = 992.4     # agent_sob/sob_constants.json
LAMT = -0.2561532; PI3 = 0.01607; G33 = 0.031105; EPS1 = 6.645e-5
a = dict(a00=55.6, a01=58.1, a02=30.8, a10=29.8, a11=38.2, a12=0.5, a20=5.9, a21=7.6, a22=6.2)
import os
MOD_ALPHA = float(os.environ.get("MODA", 268.0))          # |(a,b)| ||chi(A,B)(u0)|| / |alpha|
MOD_E2 = float(os.environ.get("MODE2", 458.4))             # E2-level modulation coupling
C3 = 12.34
cub = {"tip": json.load(open('cubic_E0_0.001_0.01_0.01.json'))["sup"]["x<=2"],
       "mid": json.load(open('cubic_E0_0.001_0.01_0.01.json'))["sup"]["2<x<=12"],
       "far": json.load(open('cubic_E0_0.001_0.01_0.01.json'))["sup"]["12<x<=31.25"]}
def Fn(al):
    al = abs(al); return 1.3e-6 * al + 1.6e-6 * al ** 2 + C3 * al ** 3
def rhs(al, A, B, D4, E):
    W1 = np.sqrt(A * B + c1 * A * A); W3 = np.sqrt(B * D4 + c1 * B * B)
    S = {}
    for reg, K in (("tip", K2), ("mid", K12)):
        S[reg] = (np.sqrt(K * (PHX * A * A + 2 * A * W1)), np.sqrt(K * (PHX * W1 * W1 + 2 * W1 * B)), np.sqrt(K * (PHX * B * B + 2 * B * W3)))
    S["far"] = (CFAR * np.sqrt(E),) * 3
    L2 = (A, W1, B)
    # E0
    lin = -(1 / 32 - PI3) * A * A - (1 / 9) * B * B
    alin = abs(al) * (a['a00'] * A * A + (a['a01'] + a['a10']) * A * W1 + (a['a02'] + a['a20']) * A * B + a['a11'] * W1 * W1
                      + (a['a12'] + a['a21']) * W1 * B + a['a22'] * B * B)
    mod = MOD_ALPHA * abs(al) * A * A
    forc = 5.72 * Fn(al) * A
    cubic = 0.0
    for reg in ("tip", "mid", "far"):
        s0, s1, s2 = S[reg]
        for key, C in cub[reg].items():
            ab = eval(key); mono = [0] * ab[0] + [1] * ab[1] + [2] * ab[2]
            mono.sort(reverse=True)                 # keep the two highest derivatives in L^2
            keep, rest = mono[:2], mono[2:]
            val = C * L2[keep[0]] * L2[keep[1]]
            for d in rest: val *= (s0, s1, s2)[d]
            cubic += 0.5 * val
    R0 = lin + alin + mod + forc + cubic
    # E2
    supw2 = S["tip"][2]
    ND = 6.9e-3 * D4 + 1.2e-2 * W3 + 4.2e-2 * B + 0.54 * A + 0.62 * W1 + 4.75 * supw2 * B + 0.136 * supw2 * W3 + 3.05 * supw2 ** 2 * B
    R2 = (-0.5 * D4 * D4 + 14.34 * B * B + 16.3 * W1 * B + 10.7 * A * B + ND * (D4 + 2.2 * W3 + 1.2 * B)
          + MOD_E2 * A * A + Fn(al) * (D4 + 2.2 * W3 + 1.2 * B))
    return R0 + EPS * R2, dict(lin=lin, alin=alin, mod=mod, forc=forc, cubic=cubic, E2=EPS * R2)
def worst(al, E):
    best = -np.inf; arg = None
    for th in np.linspace(1e-3, np.pi / 2 - 1e-3, 31):          # A = sqrt(E) cos th, B = sqrt(E/eps) sin th
        A = np.sqrt(E) * np.cos(th); B = np.sqrt(E / EPS) * np.sin(th)
        for D4 in np.concatenate([[0.0], B * np.logspace(-1, 3, 17)]):
            v, parts = rhs(al, A, B, D4, E)
            if v > best: best, arg = v, (th, D4 / max(B, 1e-300), parts)
    return best, arg
if __name__ == "__main__":
    print(f"eps = {EPS}, alpha1 = {ALPHA1}, ||w~||_1 = {W1_0}, CFAR = {CFAR}")
    for E in (1e-14, 1e-12, 1e-10, 1e-9, 1e-8):
        for al in (ALPHA1, ALPHA1 / 4, 0.0):
            b, (th, d4r, parts) = worst(al, E)
            print(f"  E={E:.0e} alpha={al:.1e}: max 1/2 E' = {b:+.3e}  (rate {b/E:+.4f}; worst B/A~{np.tan(th)/np.sqrt(EPS):.1f}, D4/B~{d4r:.2g}) "
                  f"[alin {parts['alin']/E:+.3f} mod {parts['mod']/E:+.3f} cub {parts['cubic']/E:+.2e} E2 {parts['E2']/E:+.3f} lin {parts['lin']/E:+.3f}]")
    # comparison ODE
    LE = np.linspace(-18, -8, 41); AL = np.linspace(0, ALPHA1 * 1.05, 11)
    TAB = np.array([[worst(al, 10 ** le)[0] / 10 ** le for le in LE] for al in AL])     # rate table
    from scipy.interpolate import RegularGridInterpolator
    RT = RegularGridInterpolator((AL, LE), TAB, bounds_error=False, fill_value=None)
    def f(t, y):
        al, E = y
        b = RT([[abs(al), min(max(np.log10(max(E, 1e-30)), -18), -8)]])[0] * E
        A = np.sqrt(max(E, 0))
        beta = (EPS1 * A + 0.63 * Fn(al)) / G33
        return [LAMT * al + np.sign(al) * beta, 2 * b]
    E1 = W1_0 ** 2 * (1 + EPS * 9)        # B1 ~ 3 A1
    sol = solve_ivp(f, (0, 60), [ALPHA1, E1], max_step=1.0, rtol=1e-3)
    for t, al, E in zip(sol.t[::max(1, len(sol.t)//12)], sol.y[0][::max(1, len(sol.t)//12)], sol.y[1][::max(1, len(sol.t)//12)]):
        print(f"  tau={t:5.1f}  alpha={al:+.3e}  sqrt(E)={np.sqrt(max(E,0)):.3e}")
    print("  max sqrt(E) over the run:", np.sqrt(sol.y[1].max()), " final:", np.sqrt(sol.y[1][-1]))
