"""Assemble the 16-dim CAPD input for evans_cap: profile part from the step-1 tail theorem
(gen_input.build, bracket set) and a W part supplied by the caller (linear tail theorem)."""
import sys, math
sys.path.insert(0, '/home/claude/s1_cap')
from gen_input import build, dn, up
EA, EB = "0.00015345849609375", "0.00015346679687499998"
_PROF = None
def profile_part():
    global _PROF
    if _PROF is None:
        tau1, xlo, xhi, Ct, r0lo, r0hi, info = build("8", EA, EB, 8.5)
        _PROF = (tau1, xlo[:4], xhi[:4], [row[:4] for row in Ct[:4]], r0lo[:4], r0hi[:4])
    return _PROF
def write(path, lam_box, Wlo, Whi, Wcols, Wr0lo, Wr0hi, Wrrlo, Wrrhi):
    """lam_box = (lr_lo, lr_hi, li_lo, li_hi) doubles (already outward);
    Wlo/Whi: 12 doubles (centre interval of W), Wcols: list of 12-vectors (at most 12 columns, exact doubles),
    Wr0lo/hi: parameter ranges for those columns, Wrr: extra box (12)."""
    tau1, xlo, xhi, Ct, r0lo, r0hi = profile_part()
    n = 16
    X = [(xlo[i], xhi[i]) for i in range(4)] + [(Wlo[i], Whi[i]) for i in range(12)]
    C = [[0.0]*n for _ in range(n)]
    for i in range(4):
        for j in range(4):
            C[i][j] = Ct[i][j]
    assert len(Wcols) <= 12
    for j, col in enumerate(Wcols):
        for i in range(12):
            C[4 + i][4 + j] = col[i]
    R0 = [(r0lo[i], r0hi[i]) for i in range(4)] + [(Wr0lo[j], Wr0hi[j]) for j in range(len(Wcols))] + [(0.0, 0.0)]*(12 - len(Wcols))
    RR = [(0.0, 0.0)]*4 + [(Wrrlo[i], Wrrhi[i]) for i in range(12)]
    with open(path, "w") as f:
        f.write(" ".join(repr(float(v)) for v in lam_box) + "\n")
        f.write(repr(float(tau1)) + "\n")
        for a, b in X: f.write(f"{float(a)!r} {float(b)!r}\n")
        for i in range(n): f.write(" ".join(repr(float(C[i][j])) for j in range(n)) + "\n")
        for a, b in R0: f.write(f"{float(a)!r} {float(b)!r}\n")
        for a, b in RR: f.write(f"{float(a)!r} {float(b)!r}\n")
