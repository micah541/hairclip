"""CAPD input for evans_cap2 on a lambda-box: W at sigma0 = 32 from the Lyapunov-Perron theorem (lintail2),
profile at sigma0 from the step-1 tail theorem, stage-C profile set from the step-1 shooting bracket."""
import sys, math
sys.path.insert(0, '/home/claude/s1spec/lin'); sys.path.insert(0, '/home/claude/s1_cap')
from flint import arb, acb
import lintail2 as LT
from gen_input import dn, up
from gen_evans import profile_part, EA, EB
S0 = 32

def make(box, path, R=None, eta_override=None):
    rl, rh, il, ih = box
    if R is None:
        R = LT.admissible_W(box)
    eta = R["eta"] if eta_override is None else eta_override
    prof = LT.profile_box_at(S0, EA, EB)
    W0 = R["W0"]
    mx = max(max(abs(float(w.real.mid())), abs(float(w.imag.mid()))) for w in W0)
    e2 = math.floor(math.log2(1e-8 / mx)); kap = 2.0 ** e2
    # rr for W: kappa * (eta * sum_k |dcol_k| + eta^2 |W2dir|)
    rrW = []
    for ij in range(6):
        t = arb(0)
        for dc in R["dcols"]:
            t += abs(dc[ij])
        t = eta * t + eta * eta * abs(R["W2dir"][ij])
        rrW.append(arb(kap) * t)
    X = []
    for v in prof:
        X.append((dn(v), up(v)))
    for ij in range(6):
        for part, rr in ((W0[ij].real, rrW[ij]), (W0[ij].imag, rrW[ij])):
            v = arb(kap) * part + arb(0, 1) * rr
            X.append((dn(v), up(v)))
    n = 16
    C = [[0.0] * n for _ in range(n)]
    R0 = [(0.0, 0.0)] * n
    ecol = 4
    eu = float(eta.upper()) * (1 + 1e-12)
    for col in R["cols"]:
        for mult in (1, 1j):
            vv = [z * mult for z in col]
            for ij in range(6):
                C[4 + 2 * ij][ecol] = kap * vv[ij].real
                C[4 + 2 * ij + 1][ecol] = kap * vv[ij].imag
            R0[ecol] = (-eu, eu)
            ecol += 1
    tau1, xlo, xhi, Ct, r0lo, r0hi = profile_part()
    with open(path, "w") as f:
        f.write(f"{float(rl)!r} {float(rh)!r} {float(il)!r} {float(ih)!r}\n")
        f.write(f"{float(S0 - 8)!r}\n")
        for a, b in X: f.write(f"{float(a)!r} {float(b)!r}\n")
        for i in range(n): f.write(" ".join(repr(float(C[i][j])) for j in range(n)) + "\n")
        for a, b in R0: f.write(f"{float(a)!r} {float(b)!r}\n")
        f.write(f"{float(tau1)!r}\n")
        for i in range(4): f.write(f"{float(xlo[i])!r} {float(xhi[i])!r}\n")
        for i in range(4): f.write(" ".join(repr(float(Ct[i][j])) for j in range(4)) + "\n")
        for i in range(4): f.write(f"{float(r0lo[i])!r} {float(r0hi[i])!r}\n")
    return R

if __name__ == "__main__":
    import subprocess, time
    box = tuple(float(v) for v in sys.argv[1:5])
    t0 = time.time(); R = make(box, "seg_test.txt"); t1 = time.time()
    out = subprocess.run(["./evans_cap2", "seg_test.txt", "25"], capture_output=True, text=True).stdout
    t2 = time.time()
    print(f"eta={float(R['eta'].upper()):.4g}  (tail {t1-t0:.1f}s, capd {t2-t1:.1f}s)")
    print("\n".join(l for l in out.splitlines() if l.startswith(("STAGE", "RETURN", "TIP W13", "TIP W02", "TIP p", "RESULT", "W8 W13"))))
