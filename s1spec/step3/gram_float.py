"""Non-rigorous float64 replica of the rigorous Gram pipeline (gen_hs2 split runs + gram_post4).
Same equations, same initial data (tail plane T at sigma=32, w* from wstar_coef.json, C=1, S=0, y=0),
generators truncated at sigma = 32.  Reversed time tau = 32 - sigma.
Returns G4 (4x4 in generators p, C, S, w) and G3 (constraints p, -sin(theta), w*)."""
import sys, json, pickle, numpy as np
from scipy.integrate import solve_ivp
sys.path.insert(0, '/home/claude/s1spec')
from evans import PROF, STIP

OMEGA = 3 / 64
LAMT = -0.25626
SIG0 = 32.0

def prof(s, freeze_tail=True):
    if s <= 22.0:
        return PROF(s)
    return np.array([0.0, s, 0.0, 0.0])

def coef(s, lam):
    X = s - STIP
    E = (1 + X * X) ** (-1 / 3)
    P1 = X * E / 3
    P2 = E / 3 - 2 * X * X * E ** 4 / 9
    P3 = -2 * X * E ** 4 / 3 + 16 * X ** 3 * E ** 7 / 27
    P4 = -2 * E ** 4 / 3 + 32 * X * X * E ** 7 / 9 - 224 * X ** 4 * E ** 10 / 81
    R = np.exp(-(1 + X * X) ** (2 / 3) / 4)
    p, q, k, m = prof(s)
    A = 2 * (P1 ** 2 - P2) + k * k
    V = (-P1 ** 4 / 2 + 3 * P1 ** 2 * P2 - P1 ** 2 * k * k / 2 - 2 * P1 * P3 - P1 * k * m / 2 - P1 * q / 8
         - 3 * P2 ** 2 / 2 + P2 * k * k / 2 + P4 / 2 - k * p / 2 - 3 * m * m / 2 + 3 / 8)
    return dict(P1=P1, A=A, VML=V - lam, R=R, p=p, q=q, k=k, m=m)

def Mfwd(c):
    return np.array([[0, 1, 0, 0], [0, 0, 1, 0], [0, -c['A'], c['P1'], -1], [-c['VML'], 0, 0, c['P1']]], float)

def MLT(c):
    """forward-x matrix of L^T w = lamt w on (w, w', w'', w'''):  w'''' = -k^2 w'' + (k m + q/4) w' + (1/2 - lamt) w"""
    k, m, q = c['k'], c['m'], c['q']
    return np.array([[0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1], [(0.5 - LAMT), (k * m + q / 4), -k * k, 0]], float)

def lt_frame():
    q = 32.0
    roots = np.roots([1, 0, 0, -q / 4, -(0.5 - LAMT)])
    z = [r for r in roots if r.real < -0.5][0]
    v = np.array([1, z, z * z, z ** 3])
    return v.real, v.imag

def wstar_init():
    wa, wb = lt_frame()
    al, be = json.load(open('/home/claude/s1spec/step3/wstar_coef.json'))
    return al * wa + be * wb

GEN = ["p", "C", "S", "w"]
CUT = True

def run(lam=-OMEGA, T=None, w0=None, rtol=1e-11, sig_end=None, lamt=None, gens_fun=None):
    global LAMT
    if lamt is not None: LAMT = lamt
    if T is None: T = pickle.load(open('/home/claude/s1spec/step3/tail_mid.pkl', 'rb'))['Tm']
    za0, zb0 = T[:, 2].copy(), T[:, 3].copy()
    if w0 is None: w0 = wstar_init()
    sig_end = STIP if sig_end is None else sig_end
    # state: C,S | ws(4) | za(4) zb(4) | Ia(4) Ib(4) [int za0*g_i, zb0*g_i] | y_j (16) | I_{j,i} (16)
    def unpack(Y):
        C, S = Y[0], Y[1]; ws = Y[2:6]; za = Y[6:10]; zb = Y[10:14]
        Ia = Y[14:18]; Ib = Y[18:22]; yy = Y[22:38].reshape(4, 4); II = Y[38:54].reshape(4, 4)
        return C, S, ws, za, zb, Ia, Ib, yy, II
    def f(tau, Y):
        s = SIG0 - tau
        c = coef(s, lam)
        C, S, ws, za, zb, Ia, Ib, yy, II = unpack(Y)
        M = Mfwd(c)
        dC, dS = -S * c['k'], C * c['k']
        dws = -(MLT(c) @ ws) - c['P1'] * ws
        dza = -(M @ za); dzb = -(M @ zb)
        u = min(max((SIG0 - s) / 3.0, 0.0), 1.0); chi = u**4 * (35 - 84*u + 70*u*u - 20*u**3) if CUT else 1.0
        fv = chi * np.array([c['p'], C, S, ws[0]])
        gv = fv * c['R']
        dIa = za[0] * gv; dIb = zb[0] * gv
        dyy = np.zeros((4, 4)); dII = np.zeros((4, 4))
        for j in range(4):
            dyy[j] = -(M @ yy[j]); dyy[j, 3] += fv[j]     # forward y4' = ... - f  ->  reversed: +f
            dII[j] = yy[j, 0] * gv
        return np.concatenate([[dC, dS], dws, dza, dzb, dIa, dIb, dyy.ravel(), dII.ravel()])
    Y0 = np.zeros(54); Y0[0] = 1.0; Y0[2:6] = w0; Y0[6:10] = za0; Y0[10:14] = zb0
    sol = solve_ivp(f, [0, SIG0 - sig_end], Y0, method='DOP853', rtol=rtol, atol=1e-14, dense_output=True)
    return sol, unpack

def gram(sol, unpack):
    C, S, ws, za, zb, Ia, Ib, yy, II = unpack(sol.y[:, -1])
    G4 = np.zeros((4, 4))
    for j in range(4):
        # u_j = y_j + al za + be zb,  BC at tip: index 1 and 3 vanish
        Mbc = np.array([[za[1], zb[1]], [za[3], zb[3]]])
        al, be = np.linalg.solve(Mbc, -np.array([yy[j, 1], yy[j, 3]]))
        for i in range(4):
            G4[i, j] = II[j, i] + al * Ia[i] + be * Ib[i]
    Rm = np.array([[1, 0, 0, 0], [0, C, S, 0], [0, 0, 0, 1]])
    G3 = Rm @ G4 @ Rm.T
    return G4, G3, dict(C=C, S=S, ws=ws, za=za, zb=zb)

if __name__ == "__main__":
    sol, up = run()
    G4, G3, tip = gram(sol, up)
    np.set_printoptions(linewidth=150, precision=6)
    print("G4 =\n", G4)
    print("asym:", np.abs(G4 - G4.T).max())
    Gs = (G3 + G3.T) / 2
    D = np.sqrt(abs(np.diag(Gs)))
    print("G3 =\n", Gs)
    print("eig:", np.linalg.eigvalsh(Gs), " normalized:", np.linalg.eigvalsh(Gs / np.outer(D, D)))
    print("tip:", tip)
