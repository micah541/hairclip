"""Build inputs for hsrun (validated inward integration) and post-process.
Modes:  P  - Pluecker coordinates of the H_sym admissible plane (Morse count), lam = -omega
        Z  - admissible frames Z1, Z2 + integrals against the constraint generators
        Q<j> - particular solution for forcing generator j (1..5) + integrals
Constraint generators: f1 = p, f2a = cos K, f2b = sin K (K = int_sigma^32 k), f3a = wa/rho, f3b = wb/rho,
where (wa, wb) solve L^T w = lamt w (inward from sigma = 32, approximately admissible data)."""
import sys, math, subprocess
sys.path.insert(0, '/home/claude/s1spec/lin'); sys.path.insert(0, '/home/claude/s1_cap'); sys.path.insert(0, '/home/claude/s1spec/cap')
import numpy as np
from flint import arb, arb_mat, ctx
import lintail2 as LT
from gen_input import dn, up
from gen_evans import profile_part, EA, EB
import hstail as H

ctx.prec = 200
OMEGA = arb(3) / 64
LAM = -OMEGA
LAMT = -0.25626                      # approximate third eigenvalue (any double; enters only the residual later)
TLO, THI = 10.249019, 10.249127
PAIRS = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]

P1 = "(X*E/3)"; P2 = "(E/3-2*X^2*E^4/9)"; P3 = "(-2*X*E^4/3+16*X^3*E^7/27)"
P4 = "(-2*E^4/3+32*X^2*E^7/9-224*X^4*E^10/81)"
ACO = f"(2*({P1}^2-{P2})+k^2)"
VV = (f"(-{P1}^4/2+3*{P1}^2*{P2}-{P1}^2*k^2/2-2*{P1}*{P3}-{P1}*k*m/2-{P1}*q/8-3*{P2}^2/2+{P2}*k^2/2"
      f"+{P4}/2-k*p/2-3*m^2/2+3/8)")
VML = f"({VV}-lam)"
PROF = ["k*q", "-1-k*p", "-m", "-p/4"]

def jac(y, forcing="0"):
    y0, y1, y2, y3 = y
    return [f"-{y1}", f"-{y2}", f"{ACO}*{y1}-{P1}*{y2}+{y3}", f"{VML}*{y0}-{P1}*{y3}+({forcing})"]

def ltr(w):
    w0, w1, w2, w3 = w
    return [f"-{w1}", f"-{w2}", f"-{w3}", f"-(-k^2*{w2}+(k*m+q/4)*{w1}+(1/2-lamt)*{w0})"]

def pluecker(Wn):
    # forward  W_ij' = sum_k M_ik W_kj + M_jk W_ik ; reversed time: minus
    M = [["0", "1", "0", "0"], ["0", "0", "1", "0"], ["0", f"(-{ACO})", P1, "(-1)"], [f"(-{VML})", "0", "0", P1]]
    idx = {pq: n for n, pq in enumerate(PAIRS)}
    def Wv(a, b):
        if a == b: return None
        s = 1
        if a > b: a, b, s = b, a, -1
        return (s, Wn[idx[(a, b)]])
    out = []
    for (i, j) in PAIRS:
        terms = []
        for kk in range(4):
            for (coef, pair) in ((M[i][kk], (kk, j)), (M[j][kk], (i, kk))):
                if coef == "0": continue
                wv = Wv(*pair)
                if wv is None: continue
                s, name = wv
                terms.append(f"{'+' if s > 0 else '-'}{coef}*{name}")
        out.append("-(" + ("".join(terms) if terms else "0") + ")")
    return out

import re as _re
PARS = ["lam", "lamt", "stip", "kk", "mm", "pp", "dq"]
def vf_string(var, fun, par=PARS):
    return "par:" + ",".join(par) + ";var:" + ",".join(var) + ";fun:" + ",".join(fun) + ";"
def stageB(fun):
    """profile rows frozen (dummy); profile terms -> interval parameters; q -> sigma + dq"""
    out = ["0", "0", "0", "0"]
    for e in fun[4:]:
        e = _re.sub(r"\bp\b", "pp", e); e = _re.sub(r"\bk\b", "kk", e); e = _re.sub(r"\bm\b", "mm", e)
        e = _re.sub(r"\bq\b", "(X+stip+dq)", e)
        out.append(e)
    return out
def chunks():
    rows = []
    for c in range(24):
        sb, sa = 32 - c, 31 - c
        kb, mb, pb, qs = H.profile_bounds(arb(sa))
        h = min(0.1, 3.2 / sb)
        import os
        if os.environ.get('HSTEP'): h = float(os.environ['HSTEP'])
        rows.append((1.0, h, [(-up(kb), up(kb)), (-up(mb), up(mb)), (-up(pb), up(pb)), (-up(qs), up(qs))]))
    return rows

def interval_pair(v):
    return (dn(v), up(v))

def common_init():
    prof = LT.profile_box_at(32, EA, EB)
    X0 = arb(32) - H.STIP
    E0 = (1 + X0 * X0) ** (-arb(1) / 3)
    R0 = (-(1 + X0 * X0) ** (arb(2) / 3) / 4).exp()
    return prof, X0, E0, R0

def lt_frame():
    """approximately admissible pair for L^T at sigma = 32 (k ~ 0, q ~ 32): roots of z^4 = (q/4) z + 1/2 - lamt"""
    q = 32.0
    roots = np.roots([1, 0, 0, -q / 4, -(0.5 - LAMT)])
    z = [r for r in roots if r.real < -0.5][0]
    v = np.array([1, z, z * z, z ** 3])
    sc = 2.0 ** -36
    return [float(c) for c in (sc * v.real)], [float(c) for c in (sc * v.imag)]

def write_input(path, var, fun, x, cols, r0, box, params):
    tau1, xlo, xhi, Ct, r0lo, r0hi = profile_part()
    n = len(var)
    vfC = vf_string(var, fun)
    vfB = vf_string(var, stageB(fun))
    with open(path, "w") as f:
        f.write(f"{n} {len(params)}\n")
        for name, lo, hi in params: f.write(f"{name} {lo!r} {hi!r}\n")
        f.write(vfB + "\n" + vfC + "\n")
        f.write("24.0\n")
        ch = chunks()
        f.write(f"{len(ch)} 4 kk mm pp dq\n")
        for ln, h, pv in ch:
            f.write(f"{ln!r} {h!r} " + " ".join(f"{a!r} {b!r}" for a, b in pv) + "\n")
        for lo, hi in x: f.write(f"{lo!r} {hi!r}\n")
        f.write(f"{len(cols)}\n")
        for i in range(n): f.write(" ".join(repr(float(c[i])) for c in cols) + ("\n" if cols else "\n"))
        for lo, hi in r0: f.write(f"{lo!r} {hi!r}\n")
        for lo, hi in box: f.write(f"{lo!r} {hi!r}\n")
        for i in range(4): f.write(f"{float(xlo[i])!r} {float(xhi[i])!r}\n")
        for i in range(4): f.write(" ".join(repr(float(Ct[i][j])) for j in range(4)) + "\n")
        for i in range(4): f.write(f"{float(r0lo[i])!r} {float(r0hi[i])!r}\n")
        f.write(f"{TLO!r} {THI!r}\n")

def tail_data():
    res = H.tail_theorem(LAM)
    T, zu, zs, P1v = H.plane_at_x0(LAM, res['eta'])
    return res, T

def params():
    lam = LAM
    import os
    st = (dn(H.STIP), up(H.STIP)) if not os.environ.get('STIPPOINT') else (float(H.STIP.mid()), float(H.STIP.mid()))
    return [("lam", dn(lam), up(lam)), ("lamt", LAMT, LAMT), ("stip", st[0], st[1])]

def mid_box(v):
    m = float(v.mid())
    return m, arb(v - arb(m))

def build_P(path, mon_dsub=0.004):
    res, T = tail_data()
    eta = res['eta']
    prof, X0, E0, R0 = common_init()
    var = ["p", "q", "k", "m", "X", "E"] + [f"W{i}{j}" for (i, j) in PAIRS]
    Wn = var[6:]
    fun = PROF + ["-1", "2*X*E^4/3"] + pluecker(Wn)
    def pl(i, j):
        return [T[a, i] * T[b, j] - T[b, i] * T[a, j] for (a, b) in PAIRS]
    w23 = pl(2, 3)
    colsW = [(pl(0, 3), eta), (pl(1, 3), eta), (pl(0, 2), eta), (pl(1, 2), eta), (pl(0, 1), eta * eta)]
    scale = 2.0 ** -math.floor(math.log2(max(abs(float(w.mid())) for w in w23)))
    x, box = [], []
    for v in prof: x.append((0.0, 0.0)); box.append((0.0, 0.0))
    for v in (X0, E0): m, b = mid_box(v); x.append((m, m)); box.append(interval_pair(b))
    cols, r0 = [], []
    extra = [arb(0)] * 6
    for colv, e in colsW:
        cvec = [0.0] * 12
        for ij in range(6):
            m, b = mid_box(colv[ij] * scale)
            cvec[6 + ij] = m
            extra[ij] += abs(b) * e
        cols.append(cvec); r0.append((-float(e.upper()) * (1 + 1e-12), float(e.upper()) * (1 + 1e-12)))
    for ij in range(6):
        m, b = mid_box(w23[ij] * scale)
        x.append((m, m)); box.append(interval_pair(b + arb(0, 1) * extra[ij]))
    write_input(path, var, fun, x, cols, r0, box, params())
    with open(path, "a") as f:
        f.write(f"2 {6 + 0} {6 + 1}\n{mon_dsub!r}\n")
    return res

def build_ZQ(path, mode):
    """mode 'Z' or ('Q', j) with j in 1..5"""
    res, T = tail_data()
    eta = res['eta']
    prof, X0, E0, R0 = common_init()
    wa, wb = lt_frame()
    base = ["p", "q", "k", "m", "X", "E", "R", "Ri", "C", "S"] + [f"wa{i}" for i in range(4)] + [f"wb{i}" for i in range(4)]
    bfun = PROF + ["-1", "2*X*E^4/3", "X*E*R/3", "-X*E*Ri/3", "-S*k", "C*k"] + ltr(base[10:14]) + ltr(base[14:18])
    gens = ["p*R", "C*R", "S*R", "wa0", "wb0"]          # f_i * rho  (integrand weights)
    forc = ["p", "C", "S", "wa0*Ri", "wb0*Ri"]          # forcings f_i
    x, box = [], []
    for v in prof: x.append((0.0, 0.0)); box.append((0.0, 0.0))
    for v in (X0, E0, R0, 1 / R0): m, b = mid_box(v); x.append((m, m)); box.append(interval_pair(b))
    x += [(1.0, 1.0), (0.0, 0.0)]; box += [(0.0, 0.0), (0.0, 0.0)]
    for c in wa + wb: x.append((c, c)); box.append((0.0, 0.0))
    if mode == 'Z':
        Zn = [f"za{i}" for i in range(4)] + [f"zb{i}" for i in range(4)]
        In = [f"Ia{i}" for i in range(5)] + [f"Ib{i}" for i in range(5)]
        var = base + Zn + In
        fun = bfun + jac(Zn[:4]) + jac(Zn[4:]) + [f"za0*{g}" for g in gens] + [f"zb0*{g}" for g in gens]
        n = len(var)
        cols, r0 = [], []
        extra = [arb(0)] * 8
        # Z1 = T e2 + G11 T e0 + G21 T e1 ;  Z2 = T e3 + G12 T e0 + G22 T e1
        for (target, src) in [(0, 0), (0, 1), (1, 0), (1, 1)]:
            cvec = [0.0] * n
            for i in range(4):
                m, b = mid_box(T[i, src])
                cvec[18 + 4 * target + i] = m
                extra[4 * target + i] += abs(b) * eta
            cols.append(cvec); r0.append((-float(eta.upper()) * (1 + 1e-12), float(eta.upper()) * (1 + 1e-12)))
        for target, src in [(0, 2), (1, 3)]:
            for i in range(4):
                m, b = mid_box(T[i, src])
                x.append((m, m)); box.append(interval_pair(b + arb(0, 1) * extra[4 * target + i]))
        x += [(0.0, 0.0)] * 10; box += [(0.0, 0.0)] * 10
    else:
        j = mode[1]
        Qn = [f"y{i}" for i in range(4)]
        In = [f"I{i}" for i in range(5)]
        var = base + Qn + In
        fun = bfun + jac(Qn, forcing=forc[j - 1]) + [f"y0*{g}" for g in gens]
        cols, r0 = [], []
        x += [(0.0, 0.0)] * 9; box += [(0.0, 0.0)] * 9
    write_input(path, var, fun, x, cols, r0, box, params())
    with open(path, "a") as f:
        f.write("0\n0.1\n")
    return var

def run(path, order=20, env=None):
    import os
    e = dict(os.environ); e.update(env or {})
    out = subprocess.run(["./hsrun", path, str(order)], capture_output=True, text=True, env=e).stdout
    return out

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
    mode = sys.argv[1]
    t0 = time.time()
    if mode == "P":
        res = build_P("inP.txt")
        print(f"tail: eta={float(res['eta'].upper()):.4e}", flush=True)
        out = run("inP.txt")
    else:
        md = 'Z' if mode == 'Z' else ('Q', int(mode[1:]))
        build_ZQ(f"in{mode}.txt", md)
        out = run(f"in{mode}.txt")
    open(f"out{mode}.txt", "w").write(out)
    print(out[-600:] if "RESULT OK" not in out else "RESULT OK", f"({time.time() - t0:.1f}s)")
