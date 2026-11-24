"""Inputs for hsrun3: tip parameter st as a state, rescaled stage C ending exactly at the tip, doubleton output.
Stage B (sigma = 32 - t, t in [0,24]):   X = sigma - st.
Stage C (t in [0,T0], sigma = 8 - c t, c = (8-st)/T0):   X = c (T0 - t);  every component multiplied by c.
Runs:  Z<g>  : admissible plane (Pluecker W) + minors m_r for generator g        (rows: st | aux | W(6) | m(4))
       Q<j><i>: forced solution y for generator j + integral against generator i  (rows: st | aux | y(4) | I)"""
import sys, math, os, json, subprocess
sys.path.insert(0, '/home/claude/s1spec/step3')
import numpy as np
from flint import arb, ctx
import gen_hs2 as G2
from gen_hs2 import dn, up, PAIRS, mid_box, ipair, tail_data, chunks
import hstail as H
from gen_evans import profile_part

ctx.prec = 200
OMEGA = G2.OMEGA
LAM = -OMEGA
FORM = G2.FORM
LAMT = -0.2561532224577463          # refined (refine_wstar.py): adjoint solution has F'(0) = F'''(0) = 0 in float
ST_MID = float(H.STIP.mid())
ST_RAD = float(H.STIP.rad()) * (1 + 1e-9) + 1e-18
T0 = 8.0 - ST_MID                        # exact double; c = (8 - st)/T0 ~ 1
LB = 4.0                                 # support length of the tip bumps (in rescaled time = X/c)
TS = T0 - LB                             # stage C2 = last LB units of rescaled time (tip bumps live there)
PARS = ["lam", "lamt", "T0", "kk", "mm", "pp", "dq", "ca", "cb"]
CC = "((8-st)/T0)"
NCUT = 3          # generators are multiplied by chi(sigma): chi = 1 for sigma <= 29, C^3 smoothstep to 0 at sigma = 32
U = f"(t/{NCUT})"
CHI_B = f"(ca+cb*{U}^4*(35-84*{U}+70*{U}^2-20*{U}^3))"
ZSHIFT = int(os.environ.get("ZSHIFT", "10"))
QSHIFT = int(os.environ.get("QSHIFT", "4"))

def weight_terms(stage):
    X = "((32-t)-st)" if stage == "B" else f"({CC}*(T0-t))"      # stages C and C2
    E = f"(1+{X}^2)^(-1/3)"
    P1 = f"({X}*{E}/3)"
    P2 = f"({E}/3-2*{X}^2*{E}^4/9)"
    P3 = f"(-2*{X}*{E}^4/3+16*{X}^3*{E}^7/27)"
    P4 = f"(-2*{E}^4/3+32*{X}^2*{E}^7/9-224*{X}^4*{E}^10/81)"
    R = f"exp(-(1+{X}^2)^(2/3)/4)"
    return P1, P2, P3, P4, R

def coeffs(stage):
    P1, P2, P3, P4, R = weight_terms(stage)
    if stage == "B":
        p, q, k, m = "pp", "((32-t)+dq)", "kk", "mm"   # frozen tail profile (interval parameters)
    else:
        p, q, k, m = "p", "q", "k", "m"
    A = f"(2*({P1}^2-{P2})+{k}^2)"
    V = (f"(-{P1}^4/2+3*{P1}^2*{P2}-{P1}^2*{k}^2/2-2*{P1}*{P3}-{P1}*{k}*{m}/2-{P1}*{q}/8-3*{P2}^2/2+{P2}*{k}^2/2"
         f"+{P4}/2-{k}*{p}/2-3*{m}^2/2+3/8)")
    if stage == "C2":
        sv = f"((T0-t)/{LB!r})"; zeta = f"(1-{sv}^2)^4"
        e1, e3 = f"({LB!r}*{sv}*{zeta})", f"({LB**3!r}*{sv}^3*{zeta}/6)"
    else:
        e1, e3 = "0", "0"
    VML = f"({V}-lam)"
    if FORM: A, VML = f"(({FORM})*{A})", f"(({FORM})*{VML})"
    return dict(P1=P1, A=A, VML=VML, R=R, p=p, q=q, k=k, m=m, chi=(CHI_B if stage == "B" else "1"), e1=e1, e3=e3)

PROF = ["k*q", "-1-k*p", "-m", "-p/4"]

def vf(var, fun):
    return "time:t;par:" + ",".join(PARS) + ";var:" + ",".join(var) + ";fun:" + ",".join(fun) + ";"

def scaleC(fun):
    return [f"{CC}*({f})" if f != "0" else "0" for f in fun]

def params():
    return [("lam", dn(LAM), up(LAM)), ("lamt", LAMT, LAMT), ("T0", T0, T0),
            ("kk", 0.0, 0.0), ("mm", 0.0, 0.0), ("pp", 0.0, 0.0), ("dq", 0.0, 0.0), ("ca", 1.0, 1.0), ("cb", 0.0, 0.0)]

def write(path, varB, funB, funC, funD, x, cols, r0, box, mon=(), dsub=0.1):
    """varB/funB/funC exclude st (prepended here, row 0); cols are nB-1 vectors (st column added)."""
    tau1, xlo, xhi, Ct, r0lo, r0hi = profile_part()
    varB = ["st"] + varB
    funB = ["0"] + funB
    funC = PROF + ["0"] + funC
    funD = PROF + ["0"] + funD
    varC = ["p", "q", "k", "m"] + varB
    x = [(ST_MID, ST_MID)] + list(x)
    box = [(0.0, 0.0)] + list(box)
    nB = len(varB)
    cols = [[1.0] + [0.0] * (nB - 1)] + [[0.0] + list(c) for c in cols]
    r0 = [(-ST_RAD, ST_RAD)] + list(r0)
    with open(path, "w") as f:
        ps = params()
        f.write(f"{nB} {len(ps)}\n")
        for nm, lo, hi in ps: f.write(f"{nm} {lo!r} {hi!r}\n")
        f.write(vf(varB, funB) + "\n" + vf(varC, scaleC(funC)) + "\n" + vf(varC, scaleC(funD)) + "\n")
        ch = chunks()
        f.write(f"{len(ch)} 6 kk mm pp dq ca cb\n")
        for c, (ln, h, pv) in enumerate(ch):
            assert ln == 1.0
            cut = [(0.0, 0.0), (1.0, 1.0)] if c < NCUT else [(1.0, 1.0), (0.0, 0.0)]
            f.write(f"{ln!r} {h!r} " + " ".join(f"{a!r} {b!r}" for a, b in list(pv) + cut) + "\n")
        for lo, hi in x: f.write(f"{lo!r} {hi!r}\n")
        f.write(f"{len(cols)}\n")
        for i in range(nB): f.write(" ".join(repr(float(cc[i])) for cc in cols) + "\n")
        for lo, hi in r0: f.write(f"{lo!r} {hi!r}\n")
        for lo, hi in box: f.write(f"{lo!r} {hi!r}\n")
        for i in range(4): f.write(f"{float(xlo[i])!r} {float(xhi[i])!r}\n")
        for i in range(4): f.write(" ".join(repr(float(Ct[i][j])) for j in range(4)) + "\n")
        for i in range(4): f.write(f"{float(r0lo[i])!r} {float(r0hi[i])!r}\n")
        f.write(f"{T0!r} {TS!r}\n")
        f.write(f"{len(mon)} " + " ".join(str(i) for i in mon) + f"\n{dsub!r}\n")
    return varC

# ---------------------------------------------------------------- generators
GEN = ["p", "C", "S", "w"]
GEN_E = ["e1", "e3"]
WSTAR = None
def lt_frame():
    """frozen decaying pair of L^T w = lamt w at sigma = 32 (scaled frame initial data; exact doubles once computed)"""
    roots = np.roots([1, 0, 0, -32.0 / 4, -(0.5 - LAMT)])
    z = [r for r in roots if r.real < -0.5][0]
    v = np.array([1, z, z * z, z ** 3])
    return v.real, v.imag
def wstar_init():
    wa, wb = lt_frame()
    al = json.load(open('/home/claude/s1spec/step3/wstar_refined.json'))["alpha"]
    return [float(al * a - b) for a, b in zip(wa, wb)]

def ltr(c, w):
    w0, w1, w2, w3 = w
    k, m, q, P1 = c['k'], c['m'], c['q'], c['P1']
    return [f"-{w1}-{P1}*{w0}", f"-{w2}-{P1}*{w1}", f"-{w3}-{P1}*{w2}",
            f"-(-{k}^2*{w2}+({k}*{m}+{q}/4)*{w1}+(1/2-lamt)*{w0})-{P1}*{w3}"]

def jac(c, y, forcing="0"):
    y0, y1, y2, y3 = y
    return [f"-{y1}", f"-{y2}", f"{c['A']}*{y1}-{c['P1']}*{y2}+{y3}", f"{c['VML']}*{y0}-{c['P1']}*{y3}+({forcing})"]

def aux_for(gs):
    rows = []
    if "C" in gs or "S" in gs: rows += ["C", "S"]
    if "w" in gs: rows += ["ws0", "ws1", "ws2", "ws3"]
    return rows
def aux_fun(c, rows):
    out = []
    if "C" in rows: out += [f"-S*{c['k']}", f"C*{c['k']}"]
    if "ws0" in rows: out += ltr(c, ["ws0", "ws1", "ws2", "ws3"])
    return out
def aux_init(rows):
    x = []
    if "C" in rows: x += [(1.0, 1.0), (0.0, 0.0)]
    if "ws0" in rows: x += [(cc, cc) for cc in wstar_init()]
    return x
def gexpr(c, g):
    """generator times weight; e1, e3 are the tip bumps (no far-field cutoff needed)"""
    if g in ("e1", "e3"):
        return "0" if c[g] == "0" else f"{c[g]}*{c['R']}"
    chi = f"{c['chi']}*" if c['chi'] != "1" else ""
    return chi + {"p": f"{c['p']}*{c['R']}", "C": f"C*{c['R']}", "S": f"S*{c['R']}", "w": f"ws0*{c['R']}"}[g]
def fexpr(c, g):
    if g in ("e1", "e3"): return c[g]
    chi = f"{c['chi']}*" if c['chi'] != "1" else ""
    return chi + {"p": c['p'], "C": "C", "S": "S", "w": "ws0"}[g]

def build_Z(path, gi):
    res, T = tail_data(); eta = res['eta']
    Wn = [f"W{i}{j}" for (i, j) in PAIRS]
    Mn = [f"m{r}" for r in range(4)]
    aux = aux_for([gi])
    var = aux + Wn + Mn
    def fun(c):
        g = gexpr(c, gi)
        Mrev = [["0", "-1", "0", "0"], ["0", "0", "-1", "0"], ["0", c['A'], f"(-{c['P1']})", "1"], [c['VML'], "0", "0", f"(-{c['P1']})"]]
        idx = {pq: n for n, pq in enumerate(PAIRS)}
        out = aux_fun(c, aux) + G2.pluecker(c, Wn)
        for r in range(4):
            terms = [f"({Mrev[r][k]})*{Mn[k]}" for k in range(4) if Mrev[r][k] != "0"]
            wr0 = "0" if r == 0 else f"(-{Wn[idx[(0, r)]]})"
            if not (g == "0" or wr0 == "0"): terms.append(f"{g}*{wr0}")
            out.append("(" + "+".join(terms) + ")")
        return out
    def pl(i, j):
        return [T[a, i] * T[b, j] - T[b, i] * T[a, j] for (a, b) in PAIRS]
    w23 = pl(2, 3)
    scale = 2.0 ** -math.floor(math.log2(max(abs(float(w.mid())) for w in w23)))
    scale *= 2.0 ** -ZSHIFT
    nB = len(var); na = len(aux)
    x = aux_init(aux); box = [(0.0, 0.0)] * na
    cols, r0, extra = [], [], [arb(0)] * 6
    for colv, e in [(pl(0, 3), eta), (pl(1, 3), eta), (pl(0, 2), eta), (pl(1, 2), eta), (pl(0, 1), eta * eta)]:
        cvec = [0.0] * nB
        for ij in range(6):
            m, b = mid_box(colv[ij] * scale); cvec[na + ij] = m; extra[ij] += abs(b) * e
        cols.append(cvec); r0.append((-float(e.upper()) * (1 + 1e-12), float(e.upper()) * (1 + 1e-12)))
    for ij in range(6):
        m, b = mid_box(w23[ij] * scale); x.append((m, m)); box.append(ipair(b + arb(0, 1) * extra[ij]))
    x += [(0.0, 0.0)] * 4; box += [(0.0, 0.0)] * 4
    varC = write(path, var, fun(coeffs("B")), fun(coeffs("C")), fun(coeffs("C2")), x, cols, r0, box)
    return varC

def build_Q(path, gj, gi):
    aux = aux_for([gj, gi])
    var = aux + [f"y{i}" for i in range(4)] + ["I0"]
    def fun(c):
        fx = fexpr(c, gj); gx = gexpr(c, gi)
        forcing = "0" if fx == "0" else f"{2.0 ** -QSHIFT!r}*" + fx
        return aux_fun(c, aux) + jac(c, ["y0", "y1", "y2", "y3"], forcing=forcing) + ["0" if gx == "0" else f"y0*{gx}"]
    x = aux_init(aux) + [(0.0, 0.0)] * 5
    box = [(0.0, 0.0)] * len(x)
    varC = write(path, var, fun(coeffs("B")), fun(coeffs("C")), fun(coeffs("C2")), x, [], [], box)
    return varC

def build_W(path):
    """adjoint solution alone (diagnostic)"""
    aux = aux_for(["w"])
    def fun(c): return aux_fun(c, aux)
    return write(path, aux, fun(coeffs("B")), fun(coeffs("C")), fun(coeffs("C2")), aux_init(aux), [], [], [(0.0, 0.0)] * 4)

def run(path, order=20, env=None, exe="./hsrun3"):
    e = dict(os.environ); e.update(env or {})
    return subprocess.run([exe, path, str(order)], capture_output=True, text=True, env=e, cwd='/home/claude/s1spec/step3').stdout

# ---------------------------------------------------------------- parsing of the doubleton output
def parse3(out):
    D = {}
    for l in out.splitlines():
        s = l.split()
        if not s: continue
        if s[0] in ("TIP", "TX", "TR0", "TR", "STAGEB"):
            v = list(map(float, s[1:])); D[s[0]] = [(v[2 * i], v[2 * i + 1]) for i in range(len(v) // 2)]
        elif s[0] in ("TC", "TB"):
            r, cN = int(s[1]), int(s[2]); v = list(map(float, s[3:]))
            D[s[0]] = [[(v[2 * (i * cN + j)], v[2 * (i * cN + j) + 1]) for j in range(cN)] for i in range(r)]
    D["ok"] = "RESULT OK" in out
    return D

if __name__ == "__main__":
    import time
    tag = sys.argv[1]; t0 = time.time()
    path = f"in6{tag}.txt"
    if tag == "W": varC = build_W(path)
    elif tag.startswith("Z"): varC = build_Z(path, tag[1:])
    else: varC = build_Q(path, tag[1], tag[2])
    out = run(path)
    open(f"out6{tag}.txt", "w").write(out)
    D = parse3(out)
    print(tag, "OK" if D["ok"] else out[-600:], f"({time.time() - t0:.1f}s)")

# ---------------------------------------------------------------- run for the Lyapunov constants
A_INT = ["IFF", "IpF", "ICF", "ISF", "Ipp", "IpC", "IpS", "ICC", "ICS", "ISS"]
def build_A(path):
    """aux rows C, S, ws0..ws3 plus integrals of rho * (products) over [sigma_tip, 32]:
       IFF = int rho (chi F)^2, IpF = int rho p chi F, ICF, ISF likewise with C, S;
       Ixy = int rho chi x y for x, y in {p, C, S} (the matrix <psi_i, f_j>)."""
    aux = aux_for(["C", "w"])
    def fun(c):
        chi = c['chi']; R = c['R']; P = c['p']
        F = f"({chi})*ws0" if chi != "1" else "ws0"
        ch = f"({chi})*" if chi != "1" else ""
        ints = [f"{F}^2*{R}", f"{P}*{F}*{R}", f"C*{F}*{R}", f"S*{F}*{R}",
                f"{ch}{P}^2*{R}", f"{ch}{P}*C*{R}", f"{ch}{P}*S*{R}", f"{ch}C^2*{R}", f"{ch}C*S*{R}", f"{ch}S^2*{R}"]
        return aux_fun(c, aux) + ints
    x = aux_init(aux) + [(0.0, 0.0)] * len(A_INT)
    return write(path, aux + A_INT, fun(coeffs("B")), fun(coeffs("C")), fun(coeffs("C2")), x, [], [], [(0.0, 0.0)] * len(x))

# ---------------------------------------------------------------- Lyapunov-constant run, version 2 (angle restarted at sigma = 8)
A2_INT = ["BFF", "BpF", "BCF", "BSF", "Bpp", "BpC", "BpS", "BCC", "BCS", "BSS",
          "CFF", "CpF", "CCF", "CSF", "Cpp", "CpC", "CpS", "CCC", "CCS", "CSS"]
def build_A2(path):
    """rows: C, S (angle K from sigma=32; frozen in stage C), Ct2, St2 (angle from sigma=8; frozen in stage B),
       ws0..ws3, then integrals over [8,32] (B*) and over [sigma_tip, 8] (C*) of rho * products, as in build_A."""
    aux = ["C", "S", "Ct2", "St2", "ws0", "ws1", "ws2", "ws3"]
    def fun(c, stage):
        k = c['k']; chi = c['chi']; R = c['R']; P = c['p']
        if stage == "B":
            out = [f"-S*{k}", f"C*{k}", "0", "0"]
            cs = ("C", "S")
        else:
            out = ["0", "0", f"-St2*{k}", f"Ct2*{k}"]
            cs = ("Ct2", "St2")
        out += ltr(c, ["ws0", "ws1", "ws2", "ws3"])
        F = f"({chi})*ws0" if chi != "1" else "ws0"
        ch = f"({chi})*" if chi != "1" else ""
        Cv, Sv = cs
        ints = [f"{F}^2*{R}", f"{P}*{F}*{R}", f"{Cv}*{F}*{R}", f"{Sv}*{F}*{R}",
                f"{ch}{P}^2*{R}", f"{ch}{P}*{Cv}*{R}", f"{ch}{P}*{Sv}*{R}", f"{ch}{Cv}^2*{R}", f"{ch}{Cv}*{Sv}*{R}", f"{ch}{Sv}^2*{R}"]
        z = ["0"] * 10
        return out + (ints + z if stage == "B" else z + ints)
    x = [(1.0, 1.0), (0.0, 0.0), (1.0, 1.0), (0.0, 0.0)] + [(cc, cc) for cc in wstar_init()] + [(0.0, 0.0)] * len(A2_INT)
    return write(path, aux + A2_INT, fun(coeffs("B"), "B"), fun(coeffs("C"), "C"), fun(coeffs("C2"), "C"), x, [], [], [(0.0, 0.0)] * len(x))
