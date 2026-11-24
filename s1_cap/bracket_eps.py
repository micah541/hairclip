import subprocess, re
from gen_input import build, write
def M_sign(eps, tau1=8.5):
    t, xlo, xhi, Ct, r0lo, r0hi, info = build("8", repr(eps), repr(eps), tau1)
    write("in_bis.txt", t, xlo, xhi, Ct, r0lo, r0hi)
    out = subprocess.run(["./s1proof", "in_bis.txt", "25"], capture_output=True, text=True).stdout
    mline = [l for l in out.splitlines() if l.strip().startswith("m =")][0]
    lo, hi = map(float, re.findall(r"[-+0-9.eE]+", mline.split("=",1)[1]))
    if "RESULT OK" not in out: return None, (lo, hi)
    return (1 if lo > 0 else (-1 if hi < 0 else 0)), (lo, hi)
a, b = 1.45e-4, 1.62e-4
sa, _ = M_sign(a); sb, _ = M_sign(b); assert sa == 1 and sb == -1
for it in range(60):
    c = 0.5*(a+b); sc, (lo, hi) = M_sign(c)
    if sc == 1: a = c
    elif sc == -1: b = c
    else:
        print(f"undetermined at eps={c!r}: m in [{lo:.3e},{hi:.3e}]"); break
print(f"bracket: eps* in [{a!r}, {b!r}]  width={b-a:.3e}  (relative {(b-a)/a:.2e})")
print("sign at a:", M_sign(a)[1], " sign at b:", M_sign(b)[1])
