"""Bisection for eps* on the certified sign of M(eps) = m at the first return to q = 0, tight tail data."""
import subprocess, re, sys, json
from gen_input2 import build, write
EMAX = 1.62e-4 * (1 + 1e-12)
BIN = "/home/claude/s1/s1proof"
def run(elo, ehi, path="in_bis2.txt"):
    t, xlo, xhi, Ct, r0lo, r0hi, info = build("8", repr(elo), repr(ehi), 8.5, emax=EMAX)
    write(path, t, xlo, xhi, Ct, r0lo, r0hi)
    return subprocess.run([BIN, path, "25"], capture_output=True, text=True).stdout
def M_sign(eps):
    out = run(eps, eps)
    mline = [l for l in out.splitlines() if l.strip().startswith("m =")][0]
    lo, hi = map(float, re.findall(r"[-+0-9.eE]+", mline.split("=", 1)[1]))
    if "RESULT OK" not in out:
        return None, (lo, hi)
    return (1 if lo > 0 else (-1 if hi < 0 else 0)), (lo, hi)
if __name__ == "__main__":
    a, b = 1.45e-4, 1.62e-4
    sa, _ = M_sign(a); sb, _ = M_sign(b); assert sa == 1 and sb == -1
    for it in range(80):
        # regula-falsi-safeguarded bisection: try the secant point first when it is well inside
        c = 0.5 * (a + b)
        sc, (lo, hi) = M_sign(c)
        if sc == 1: a = c
        elif sc == -1: b = c
        else:
            print(f"undetermined at eps={c!r}: m in [{lo:.3e},{hi:.3e}]  (it {it})"); break
        if b - a < 1e-16: break
    print(f"bracket: eps* in [{a!r}, {b!r}]  width={b-a:.3e}  (relative {(b-a)/a:.2e})")
    ra, rb = M_sign(a), M_sign(b)
    print("M at a:", ra, "  M at b:", rb)
    json.dump(dict(eps_a=a, eps_b=b, M_a=ra[1], M_b=rb[1]), open("bracket2.json", "w"), indent=1)
