"""Input generation and run of evans_seg on a lambda-segment [lam_a, lam_b]."""
import sys, math, subprocess, re
sys.path.insert(0, '/home/claude/s1spec/lin'); sys.path.insert(0, '/home/claude/s1_cap')
from flint import arb, acb
import lintail2 as LT
from gen_input import dn, up
from gen_evans import profile_part, EA, EB
S0 = 32
WSCALE = 1e-8
_PROF0 = None

def make(lam_a, lam_b, path, box=None, eta=None):
    """segment run (box=None) or lambda-box run (box = (rl, rh, il, ih); lam_a/lam_b ignored)"""
    global _PROF0
    if box is None:
        la, lb = complex(lam_a), complex(lam_b)
        R = LT.segment_W(la, lb, eta=eta)
        first = f"{la.real!r} {la.real!r} {la.imag!r} {la.imag!r} {(lb - la).real!r} {(lb - la).imag!r}\n"
    else:
        R = LT.box_W(box)
        rl, rh, il, ih = box
        first = f"{float(rl)!r} {float(rh)!r} {float(il)!r} {float(ih)!r} 0.0 0.0\n"
    if _PROF0 is None:
        _PROF0 = LT.profile_box_at(S0, EA, EB)
    eta = R["eta"]
    W0c = R["W0c"]
    mx = max(max(abs(float(w.real.mid())), abs(float(w.imag.mid()))) for w in W0c)
    kap = 2.0 ** math.floor(math.log2(WSCALE / mx))
    rrW = []
    for ij in range(6):
        t = arb(0)
        for dc in R["dcols"]:
            t += abs(dc[ij])
        # det(Gamma) W2dir is carried as two parameter columns (|det Gamma| <= eta^2); only the variation of W2dir
        # over the segment box goes into the box part
        w2c = R["W2dir_c"][ij]
        t = eta * t + eta * eta * abs(R["W2dir"][ij] - acb(w2c.real, w2c.imag)) + abs(R["rem"][ij])
        rrW.append(arb(kap) * t)
    n = 17
    # profile box at sigma0 -> centre + 4 parameter columns (keeps the profile uncertainty out of the Lohner box part)
    X = []
    for i, v in enumerate(_PROF0):
        m = float(v.mid())
        X.append((m, m))
    for ij in range(6):
        for part in (W0c[ij].real, W0c[ij].imag):
            v = arb(kap) * part + arb(0, 1) * rrW[ij]
            X.append((dn(v), up(v)))
    X.append((0.5, 0.5))
    C = [[0.0] * n for _ in range(n)]
    R0 = [(0.0, 0.0)] * n
    for i, v in enumerate(_PROF0):
        m = float(v.mid())
        rad = max(up(v) - m, m - dn(v)) * (1 + 1e-12) + 1e-300
        C[i][i] = 1.0; R0[i] = (-rad, rad)
    ecol = 4
    eu = float(eta.upper()) * (1 + 1e-12)
    for colv in R["cols"]:
        for mult in (1, 1j):
            for ij in range(6):
                z = colv[ij] * mult
                C[4 + 2 * ij][ecol] = kap * z.real
                C[4 + 2 * ij + 1][ecol] = kap * z.imag
            R0[ecol] = (-eu, eu); ecol += 1
    # t column: exact doubles kap*mid(dWt); the rounding error of the column times |t-1/2| <= 1/2 goes into the centre box
    for ij in range(6):
        d = R["dWt"][ij]
        cr, ci = kap * float(d.real.mid()), kap * float(d.imag.mid())
        C[4 + 2 * ij][ecol] = cr; C[4 + 2 * ij + 1][ecol] = ci
        er = (arb(kap) * d.real - arb(cr)) * arb(0, 0.5); ei = (arb(kap) * d.imag - arb(ci)) * arb(0, 0.5)
        a, b = X[4 + 2 * ij]; X[4 + 2 * ij] = (dn(arb(a) + er), up(arb(b) + er))
        a, b = X[4 + 2 * ij + 1]; X[4 + 2 * ij + 1] = (dn(arb(a) + ei), up(arb(b) + ei))
    C[16][ecol] = 1.0; R0[ecol] = (-0.5, 0.5); ecol += 1
    e2 = float((eta * eta).upper()) * (1 + 1e-12)
    for mult in (1, 1j):
        for ij in range(6):
            z = R["W2dir_c"][ij] * mult
            C[4 + 2 * ij][ecol] = kap * z.real
            C[4 + 2 * ij + 1][ecol] = kap * z.imag
        R0[ecol] = (-e2, e2); ecol += 1
    tau1, xlo, xhi, Ct, r0lo, r0hi = profile_part()
    with open(path, "w") as f:
        f.write(first)
        f.write(f"{float(S0 - 8)!r}\n")
        for a, b in X: f.write(f"{float(a)!r} {float(b)!r}\n")
        for i in range(n): f.write(" ".join(repr(float(C[i][j])) for j in range(n)) + "\n")
        for a, b in R0: f.write(f"{float(a)!r} {float(b)!r}\n")
        f.write(f"{float(tau1)!r}\n")
        for i in range(4): f.write(f"{float(xlo[i])!r} {float(xhi[i])!r}\n")
        for i in range(4): f.write(" ".join(repr(float(Ct[i][j])) for j in range(4)) + "\n")
        for i in range(4): f.write(f"{float(r0lo[i])!r} {float(r0hi[i])!r}\n")
    return R

def run(lam_a, lam_b, path="seg_run.txt", order=25, box=None, eta=None):
    R = make(lam_a, lam_b, path, box=box, eta=eta)
    out = subprocess.run(["./evans_seg", path, str(order)], capture_output=True, text=True).stdout
    if "RESULT OK" not in out:
        return None, out, R
    vals = {}
    for l in out.splitlines():
        if l.startswith("TIP "):
            _, name, lo, hi = l.split(); vals[name] = (float(lo), float(hi))
    De = (vals["W13r"], vals["W13i"]); Do = (vals["W02r"], vals["W02i"])
    return dict(even=De, odd=Do, tip=vals), out, R

if __name__ == "__main__":
    import time
    la, lb = complex(sys.argv[1]), complex(sys.argv[2])
    t0 = time.time(); res, out, R = run(la, lb); t1 = time.time()
    print(f"eta={float(R['eta'].upper()):.4g} time {t1-t0:.1f}s")
    if res is None:
        print(out)
    else:
        for k in ("even", "odd"):
            (a, b), (c, d) = res[k]
            print(f"{k}: Re [{a:.6e}, {b:.6e}]  Im [{c:.6e}, {d:.6e}]   rel.width {max(b-a, d-c)/max(abs(a), abs(b), abs(c), abs(d)):.2e}")
