"""Inputs for hsrun2 (time-dependent vector fields).  sigma = 32 - t (stage B), 8 - t (stage C);
x = sigma - stip;  weight terms explicit in x;  stage B profile terms -> interval parameters (step-1 tail bounds).
Modes:  P (Pluecker of the H_sym admissible plane at lam = -omega: Morse count),
        Z (admissible frames + integrals against the constraint generators),  Q1..Q5 (forced solutions)."""
import sys, math, os, re, subprocess
sys.path.insert(0, '/home/claude/s1spec/lin'); sys.path.insert(0, '/home/claude/s1_cap'); sys.path.insert(0, '/home/claude/s1spec/cap')
import numpy as np
from flint import arb, ctx
import lintail2 as LT
from gen_input import dn, up
from gen_evans import profile_part, EA, EB
import hstail as H

ctx.prec = 200
_OM = os.environ.get('OMEGA_Q', '3/64').split('/'); OMEGA = arb(int(_OM[0])) / int(_OM[1])
LAM = -OMEGA
FORM = os.environ.get('FORM_MULT', '')        # e.g. '9/8' = 1/(1-theta) for the theta-form (1-theta) u''^2 - A u'^2 - (V+omega) u^2
LAMT = -0.25626
TLO, THI = 10.249019, 10.249127
PAIRS = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
PARS = ["lam", "lamt", "stip", "kk", "mm", "pp", "dq"]

def weight_terms(SIG):
    X = f"({SIG}-stip)"
    E = f"(1+{X}^2)^(-1/3)"
    P1 = f"({X}*{E}/3)"
    P2 = f"({E}/3-2*{X}^2*{E}^4/9)"
    P3 = f"(-2*{X}*{E}^4/3+16*{X}^3*{E}^7/27)"
    P4 = f"(-2*{E}^4/3+32*{X}^2*{E}^7/9-224*{X}^4*{E}^10/81)"
    R = f"exp(-(1+{X}^2)^(2/3)/4)"
    Ri = f"exp((1+{X}^2)^(2/3)/4)"
    return P1, P2, P3, P4, R, Ri

def coeffs(stage):
    SIG = "(32-t)" if stage == "B" else "(8-t)"
    P1, P2, P3, P4, R, Ri = weight_terms(SIG)
    if stage == "B":
        p, q, k, m = "pp", f"({SIG}+dq)", "kk", "mm"
    else:
        p, q, k, m = "p", "q", "k", "m"
    A = f"(2*({P1}^2-{P2})+{k}^2)"
    V = (f"(-{P1}^4/2+3*{P1}^2*{P2}-{P1}^2*{k}^2/2-2*{P1}*{P3}-{P1}*{k}*{m}/2-{P1}*{q}/8-3*{P2}^2/2+{P2}*{k}^2/2"
         f"+{P4}/2-{k}*{p}/2-3*{m}^2/2+3/8)")
    VML = f"({V}-lam)"
    if FORM: A, VML = f"(({FORM})*{A})", f"(({FORM})*{VML})"
    return dict(P1=P1, A=A, VML=VML, R=R, Ri=Ri, p=p, q=q, k=k, m=m)

def jac(c, y, forcing="0"):
    y0, y1, y2, y3 = y
    return [f"-{y1}", f"-{y2}", f"{c['A']}*{y1}-{c['P1']}*{y2}+{y3}", f"{c['VML']}*{y0}-{c['P1']}*{y3}+({forcing})"]

def ltr(c, w):
    # scaled frame wt = (w, dw, d2w, d3w)/rho of L^T w = lamt w :  d wt/dtau = M_LT wt - P1 wt
    w0, w1, w2, w3 = w
    k, m, q, P1 = c['k'], c['m'], c['q'], c['P1']
    return [f"-{w1}-{P1}*{w0}", f"-{w2}-{P1}*{w1}", f"-{w3}-{P1}*{w2}",
            f"-(-{k}^2*{w2}+({k}*{m}+{q}/4)*{w1}+(1/2-lamt)*{w0})-{P1}*{w3}"]

def pluecker(c, Wn):
    M = [["0", "1", "0", "0"], ["0", "0", "1", "0"], ["0", f"(-{c['A']})", c['P1'], "(-1)"], [f"(-{c['VML']})", "0", "0", c['P1']]]
    idx = {pq: n for n, pq in enumerate(PAIRS)}
    out = []
    for (i, j) in PAIRS:
        terms = []
        for kk in range(4):
            for coef, (a, b) in ((M[i][kk], (kk, j)), (M[j][kk], (i, kk))):
                if coef == "0" or a == b: continue
                s = 1
                if a > b: a, b, s = b, a, -1
                terms.append(f"{'+' if s > 0 else '-'}{coef}*{Wn[idx[(a, b)]]}")
        out.append("-(" + "".join(terms) + ")")
    return out

PROF = ["k*q", "-1-k*p", "-m", "-p/4"]

def vf(var, fun):
    return "time:t;par:" + ",".join(PARS) + ";var:" + ",".join(var) + ";fun:" + ",".join(fun) + ";"

def chunks():
    rows = []
    hfix = os.environ.get('HSTEP')
    for c in range(24):
        sb, sa = 32 - c, 31 - c
        kb, mb, pb, qs = H.profile_bounds(arb(sa))
        h = float(hfix) if hfix else 0.05
        rows.append((1.0, h, [(-up(kb), up(kb)), (-up(mb), up(mb)), (-up(pb), up(pb)), (-up(qs), up(qs))]))
    return rows

def params():
    return [("lam", dn(LAM), up(LAM)), ("lamt", LAMT, LAMT), ("stip", dn(H.STIP), up(H.STIP)),
            ("kk", 0.0, 0.0), ("mm", 0.0, 0.0), ("pp", 0.0, 0.0), ("dq", 0.0, 0.0)]

def mid_box(v):
    m = float(v.mid())
    return m, arb(v - arb(m))

def ipair(v):
    return (dn(v), up(v))

def write(path, varB, funB, varC, funC, x, cols, r0, box, mon=(), dsub=0.1):
    tau1, xlo, xhi, Ct, r0lo, r0hi = profile_part()
    nB = len(varB)
    with open(path, "w") as f:
        ps = params()
        f.write(f"{nB} {len(ps)}\n")
        for nm, lo, hi in ps: f.write(f"{nm} {lo!r} {hi!r}\n")
        f.write(vf(varB, funB) + "\n" + vf(varC, funC) + "\n")
        ch = chunks()
        f.write(f"{len(ch)} 4 kk mm pp dq\n")
        for ln, h, pv in ch: f.write(f"{ln!r} {h!r} " + " ".join(f"{a!r} {b!r}" for a, b in pv) + "\n")
        for lo, hi in x: f.write(f"{lo!r} {hi!r}\n")
        f.write(f"{len(cols)}\n")
        for i in range(nB): f.write(" ".join(repr(float(cc[i])) for cc in cols) + "\n")
        for lo, hi in r0: f.write(f"{lo!r} {hi!r}\n")
        for lo, hi in box: f.write(f"{lo!r} {hi!r}\n")
        for i in range(4): f.write(f"{float(xlo[i])!r} {float(xhi[i])!r}\n")
        for i in range(4): f.write(" ".join(repr(float(Ct[i][j])) for j in range(4)) + "\n")
        for i in range(4): f.write(f"{float(r0lo[i])!r} {float(r0hi[i])!r}\n")
        f.write(f"{TLO!r} {THI!r}\n")
        f.write(f"{len(mon)} " + " ".join(str(i) for i in mon) + f"\n{dsub!r}\n")

_TAIL = None
def tail_data():
    global _TAIL
    if _TAIL is None:
        res = H.tail_theorem(LAM)
        T, zu, zs, P1v = H.plane_at_x0(LAM, res['eta'])
        _TAIL = (res, T)
    return _TAIL

def build_P(path, dsub=0.004):
    res, T = tail_data(); eta = res['eta']
    Wn = [f"W{i}{j}" for (i, j) in PAIRS]
    cB, cC = coeffs("B"), coeffs("C")
    funB = pluecker(cB, Wn)
    funC = PROF + pluecker(cC, Wn)
    def pl(i, j):
        return [T[a, i] * T[b, j] - T[b, i] * T[a, j] for (a, b) in PAIRS]
    w23 = pl(2, 3)
    scale = 2.0 ** -math.floor(math.log2(max(abs(float(w.mid())) for w in w23)))
    cols, r0, extra = [], [], [arb(0)] * 6
    for colv, e in [(pl(0, 3), eta), (pl(1, 3), eta), (pl(0, 2), eta), (pl(1, 2), eta), (pl(0, 1), eta * eta)]:
        cvec = [0.0] * 6
        for ij in range(6):
            m, b = mid_box(colv[ij] * scale); cvec[ij] = m; extra[ij] += abs(b) * e
        cols.append(cvec); r0.append((-float(e.upper()) * (1 + 1e-12), float(e.upper()) * (1 + 1e-12)))
    x, box = [], []
    for ij in range(6):
        m, b = mid_box(w23[ij] * scale); x.append((m, m)); box.append(ipair(b + arb(0, 1) * extra[ij]))
    write(path, Wn, funB, ["p", "q", "k", "m"] + Wn, funC, x, cols, r0, box, mon=(4, 5), dsub=dsub)
    return res

def lt_frame():
    q = 32.0
    roots = np.roots([1, 0, 0, -q / 4, -(0.5 - LAMT)])
    z = [r for r in roots if r.real < -0.5][0]
    v = np.array([1, z, z * z, z ** 3]); sc = 1.0
    return [float(c) for c in (sc * v.real)], [float(c) for c in (sc * v.imag)]

def build_ZQ(path, mode):
    res, T = tail_data(); eta = res['eta']
    wa, wb = lt_frame()
    base = ["C", "S"] + [f"wa{i}" for i in range(4)] + [f"wb{i}" for i in range(4)]
    def basefun(c):
        return [f"-S*{c['k']}", f"C*{c['k']}"] + ltr(c, base[2:6]) + ltr(c, base[6:10])
    def gens(c):
        return [f"{c['p']}*{c['R']}", f"C*{c['R']}", f"S*{c['R']}", f"wa0*{c['R']}", f"wb0*{c['R']}"]
    def forc(c):
        return [c['p'], "C", "S", "wa0", "wb0"]
    x, box = [(1.0, 1.0), (0.0, 0.0)], [(0.0, 0.0), (0.0, 0.0)]
    for cc in wa + wb: x.append((cc, cc)); box.append((0.0, 0.0))
    if mode == 'Z':
        Zn = [f"za{i}" for i in range(4)] + [f"zb{i}" for i in range(4)]
        In = [f"Ia{i}" for i in range(5)] + [f"Ib{i}" for i in range(5)]
        var = base + Zn + In
        def fun(c):
            return basefun(c) + jac(c, Zn[:4]) + jac(c, Zn[4:]) + [f"za0*{g}" for g in gens(c)] + [f"zb0*{g}" for g in gens(c)]
        nB = len(var)
        cols, r0, extra = [], [], [arb(0)] * 8
        for (target, src) in [(0, 0), (0, 1), (1, 0), (1, 1)]:
            cvec = [0.0] * nB
            for i in range(4):
                m, b = mid_box(T[i, src]); cvec[10 + 4 * target + i] = m; extra[4 * target + i] += abs(b) * eta
            cols.append(cvec); r0.append((-float(eta.upper()) * (1 + 1e-12), float(eta.upper()) * (1 + 1e-12)))
        for target, src in [(0, 2), (1, 3)]:
            for i in range(4):
                m, b = mid_box(T[i, src]); x.append((m, m)); box.append(ipair(b + arb(0, 1) * extra[4 * target + i]))
        x += [(0.0, 0.0)] * 10; box += [(0.0, 0.0)] * 10
    else:
        j = mode[1]
        Qn = [f"y{i}" for i in range(4)]; In = [f"I{i}" for i in range(5)]
        var = base + Qn + In
        def fun(c):
            return basefun(c) + jac(c, Qn, forcing=forc(c)[j - 1]) + [f"y0*{g}" for g in gens(c)]
        cols, r0 = [], []
        x += [(0.0, 0.0)] * 9; box += [(0.0, 0.0)] * 9
    write(path, var, fun(coeffs("B")), ["p", "q", "k", "m"] + var, PROF + fun(coeffs("C")), x, cols, r0, box)
    return var

def run(path, order=20, env=None):
    e = dict(os.environ); e.update(env or {})
    return subprocess.run(["./hsrun2", path, str(order)], capture_output=True, text=True, env=e).stdout

def parse(out):
    tip, mon = None, []
    for l in out.splitlines():
        if l.startswith("TIP "):
            v = list(map(float, l.split()[1:])); tip = [(v[2 * i], v[2 * i + 1]) for i in range(len(v) // 2)]
        elif l.startswith("MON "):
            s = l.split(); mon.append((s[1], float(s[2]), float(s[3]), list(map(float, s[4:]))))
    return tip, mon

if __name__ == "__main__":
    import time
    mode = sys.argv[1]; t0 = time.time()
    if mode == "P":
        build_P("in2P.txt"); out = run("in2P.txt")
    else:
        md = 'Z' if mode == 'Z' else ('Q', int(mode[1:]))
        build_ZQ(f"in2{mode}.txt", md); out = run(f"in2{mode}.txt")
    open(f"out2{mode}.txt", "w").write(out)
    print(("RESULT OK" if "RESULT OK" in out else out[-800:]), f"({time.time() - t0:.1f}s)")

# ---------------------------------------------------------------- minor formulation (basis-invariant)
def wstar_init():
    """initial data of the single adjoint solution w* = al*wa + be*wb (al, be from the previous tip midpoints, so that
    psi*'(tip) ~ 0); stored as exact doubles."""
    import json
    wa, wb = lt_frame()
    al, be = json.load(open('wstar_coef.json'))
    return [al * a + be * b for a, b in zip(wa, wb)]

def build_ZM(path):
    res, T = tail_data(); eta = res['eta']
    Wn = [f"W{i}{j}" for (i, j) in PAIRS]
    Mn = [[f"m{r}{i}" for i in range(4)] for r in range(4)]          # m_{r,i}, generators i = p, C, S, w*
    base = ["C", "S"] + [f"ws{i}" for i in range(4)]
    var = base + Wn + [Mn[r][i] for r in range(4) for i in range(4)]
    def fun(c):
        g = [f"{c['p']}*{c['R']}", f"C*{c['R']}", f"S*{c['R']}", f"ws0*{c['R']}"]
        Mrev = [["0", "-1", "0", "0"], ["0", "0", "-1", "0"], ["0", c['A'], f"(-{c['P1']})", "1"], [c['VML'], "0", "0", f"(-{c['P1']})"]]
        idx = {pq: n for n, pq in enumerate(PAIRS)}
        def Wr0(r):          # W_{r0} = za_r zb_0 - zb_r za_0 = -W_{0r}
            return "0" if r == 0 else f"(-{Wn[idx[(0, r)]]})"
        out = [f"-S*{c['k']}", f"C*{c['k']}"] + ltr(c, base[2:6]) + pluecker(c, Wn)
        for r in range(4):
            for i in range(4):
                terms = [f"+{Mrev[r][k]}*{Mn[k][i]}" for k in range(4) if Mrev[r][k] != "0"]
                out.append("(" + "".join(terms) + f"+{g[i]}*{Wr0(r)})")
        return out
    def pl(i, j):
        return [T[a, i] * T[b, j] - T[b, i] * T[a, j] for (a, b) in PAIRS]
    w23 = pl(2, 3)
    scale = 2.0 ** -math.floor(math.log2(max(abs(float(w.mid())) for w in w23)))
    nB = len(var)
    x, box = [(1.0, 1.0), (0.0, 0.0)], [(0.0, 0.0), (0.0, 0.0)]
    for cc in wstar_init(): x.append((cc, cc)); box.append((0.0, 0.0))
    cols, r0, extra = [], [], [arb(0)] * 6
    for colv, e in [(pl(0, 3), eta), (pl(1, 3), eta), (pl(0, 2), eta), (pl(1, 2), eta), (pl(0, 1), eta * eta)]:
        cvec = [0.0] * nB
        for ij in range(6):
            m, b = mid_box(colv[ij] * scale); cvec[6 + ij] = m; extra[ij] += abs(b) * e
        cols.append(cvec); r0.append((-float(e.upper()) * (1 + 1e-12), float(e.upper()) * (1 + 1e-12)))
    for ij in range(6):
        m, b = mid_box(w23[ij] * scale); x.append((m, m)); box.append(ipair(b + arb(0, 1) * extra[ij]))
    x += [(0.0, 0.0)] * 16; box += [(0.0, 0.0)] * 16
    write(path, var, fun(coeffs("B")), ["p", "q", "k", "m"] + var, PROF + fun(coeffs("C")), x, cols, r0, box)
    return var

def build_QS(path, j):
    """forced solution for generator j in (p, C, S, w*) and its integrals against all four generators"""
    base = ["C", "S"] + [f"ws{i}" for i in range(4)]
    Qn = [f"y{i}" for i in range(4)]; In = [f"I{i}" for i in range(4)]
    var = base + Qn + In
    def fun(c):
        g = [f"{c['p']}*{c['R']}", f"C*{c['R']}", f"S*{c['R']}", f"ws0*{c['R']}"]
        f = [c['p'], "C", "S", "ws0"]
        return [f"-S*{c['k']}", f"C*{c['k']}"] + ltr(c, base[2:6]) + jac(c, Qn, forcing=f[j]) + [f"y0*{gg}" for gg in g]
    x, box = [(1.0, 1.0), (0.0, 0.0)], [(0.0, 0.0), (0.0, 0.0)]
    for cc in wstar_init(): x.append((cc, cc)); box.append((0.0, 0.0))
    x += [(0.0, 0.0)] * 8; box += [(0.0, 0.0)] * 8
    write(path, var, fun(coeffs("B")), ["p", "q", "k", "m"] + var, PROF + fun(coeffs("C")), x, [], [], box)
    return var

# ---------------------------------------------------------------- split runs (small Lohner sets)
GEN = ["p", "C", "S", "w"]          # generators f: p, cos K, sin K, w* (= scaled adjoint solution component 0)
def aux_for(gens_needed):
    rows = []
    if "C" in gens_needed or "S" in gens_needed: rows += ["C", "S"]
    if "w" in gens_needed: rows += [f"ws{i}" for i in range(4)]
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
    return {"p": f"{c['p']}*{c['R']}", "C": f"C*{c['R']}", "S": f"S*{c['R']}", "w": f"ws0*{c['R']}"}[g]
def fexpr(c, g):
    return {"p": c['p'], "C": "C", "S": "S", "w": "ws0"}[g]

def build_ZMi(path, gi):
    res, T = tail_data(); eta = res['eta']
    Wn = [f"W{i}{j}" for (i, j) in PAIRS]
    Mn = [f"m{r}" for r in range(4)]
    aux = aux_for([gi])
    var = aux + Wn + Mn
    def fun(c):
        g = gexpr(c, gi)
        Mrev = [["0", "-1", "0", "0"], ["0", "0", "-1", "0"], ["0", c['A'], f"(-{c['P1']})", "1"], [c['VML'], "0", "0", f"(-{c['P1']})"]]
        idx = {pq: n for n, pq in enumerate(PAIRS)}
        out = aux_fun(c, aux) + pluecker(c, Wn)
        for r in range(4):
            terms = [f"+{Mrev[r][k]}*{Mn[k]}" for k in range(4) if Mrev[r][k] != "0"]
            wr0 = "0" if r == 0 else f"(-{Wn[idx[(0, r)]]})"
            out.append("(" + "".join(terms) + f"+{g}*{wr0})")
        return out
    def pl(i, j):
        return [T[a, i] * T[b, j] - T[b, i] * T[a, j] for (a, b) in PAIRS]
    w23 = pl(2, 3)
    scale = 2.0 ** -math.floor(math.log2(max(abs(float(w.mid())) for w in w23)))
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
    write(path, var, fun(coeffs("B")), ["p", "q", "k", "m"] + var, PROF + fun(coeffs("C")), x, cols, r0, box)
    return var

def build_Qji(path, gj, gi):
    aux = aux_for([gj, gi])
    var = aux + [f"y{i}" for i in range(4)] + ["I0"]
    def fun(c):
        return aux_fun(c, aux) + jac(c, ["y0", "y1", "y2", "y3"], forcing=fexpr(c, gj)) + [f"y0*{gexpr(c, gi)}"]
    x = aux_init(aux) + [(0.0, 0.0)] * 5
    box = [(0.0, 0.0)] * len(x)
    write(path, var, fun(coeffs("B")), ["p", "q", "k", "m"] + var, PROF + fun(coeffs("C")), x, [], [], box)
    return var
