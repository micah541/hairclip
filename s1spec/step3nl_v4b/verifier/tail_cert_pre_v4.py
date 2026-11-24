"""Machine certification of the STRUCTURAL facts used by the tail lemma / branch-and-bound (TAIL_LEMMA.tex, Lemmas 1-3), by running the
verbatim rhs code (rhs_iv.rhs_pn) on 'degree' objects instead of intervals.  Each Deg carries (lo, hi) = (valuation, top degree) w.r.t. a scaling
x -> lambda x (lambda in (0,1]):
   G_nu : f(lambda.) <= lambda^nu f     (nu = lo)      sum: min, product: add, sqrt: halve, max: min, constants: 0
   F_p  : f(lambda.) >= lambda^p  f     (p  = hi)      sum: max, product: add, sqrt: halve, max: max, constants: 0
(both are preserved by +,*,sqrt,max of nonnegative functions: proofs in TAIL_LEMMA.tex).  No sampling is involved.
Modes (which variables scale like lambda):
   radial : only D4, V4 (the unbounded directions)   -> need hi(H) <= 2 and hi(L) <= 2   (POS(r)/r^2 nonincreasing along rays in (D4,V4))
   s      : A,B,D4,V0,V2,V4 (fixed alpha)            -> need lo(H) >= 2 (H/s^2 nondecreasing) and hi(L) <= 2 (L/s^2 nonincreasing)
   joint  : alpha and all six x                      -> need lo(H) >= 2 and lo(L) >= 2 (QUAD: no floors) so that tot(s alpha, s x) <= s^2 tot(alpha,x)
NEG is a sum of n_k x_k^2 (exactly degree 2 in x), checked as (2,2) in mode s/joint."""
import os, sys, numpy as np
from ival import P, cc
import ival, rhs_iv as R, constants as K

class Deg(P):
    __slots__ = ("dl", "dh")
    def __init__(s, dl, dh): s.dl = dl; s.dh = dh; s.lo = np.array([1e-4]); s.hi = np.array([1e-4])      # lo/hi dummies (alpha level lookup)
    @staticmethod
    def of(b): return b if isinstance(b, Deg) else Deg(0, 0)
    def __add__(a, b): b = Deg.of(b); return Deg(min(a.dl, b.dl), max(a.dh, b.dh))
    __radd__ = __add__
    def __mul__(a, b): b = Deg.of(b); return Deg(a.dl + b.dl, a.dh + b.dh)
    __rmul__ = __mul__
    def __truediv__(a, b): return Deg(a.dl, a.dh)          # division by a (positive) CONSTANT only
    def _dsqrt(a): return Deg(a.dl / 2, a.dh / 2)
    def _dmx(a, b): b = Deg.of(b); return Deg(min(a.dl, b.dl), max(a.dh, b.dh))

def run(mode, floors=True):
    one = lambda: Deg(1, 1); zero = lambda: Deg(0, 0)
    sc = {"radial": dict(al=0, x=(0, 0, 1, 0, 0, 1)), "s": dict(al=0, x=(1,) * 6), "joint": dict(al=1, x=(1,) * 6)}[mode]
    al = Deg(sc["al"], sc["al"]); xs = [Deg(e, e) for e in sc["x"]]
    H, L, N = R.rhs_pn(al, *xs, floors=floors)
    return (H.dl, H.dh), (L.dl, L.dh), (N.dl, N.dh)

if __name__ == "__main__":
    out = []
    def log(s): print(s); out.append(s)
    log(f"ADD={sorted(K.ADD)} QUAD={K.QUAD} RAMP={K.RAMP}  (degrees: (lo=valuation, hi=top degree))")
    ok = True
    for mode in ("radial", "s", "joint"):
        H, L, N = run(mode)
        log(f"mode {mode:6s}: H={H}  L={L}  NEG={N}")
        if mode == "radial": c = H[1] <= 2 and L[1] <= 2
        elif mode == "s": c = H[0] >= 2 and L[1] <= 2 and N == (2, 2)
        else:
            H2, L2, N2 = run("joint", floors=False)          # tail lemma is applied to tot_s = tot without the degree-0 (floor) constants
            log(f"   joint, floors=False (tot_s): H={H2} L={L2}")
            c = H2[0] >= 2 and L2[0] >= 2 and N2 == (2, 2)
        log(f"   required condition for mode {mode}: {'OK' if c else 'FAIL'}"); ok &= c
    log("Split (version W): homogeneous part (floors=False, degrees >= 2 jointly in (alpha,x)) -> mode joint above, i.e. the scale-invariant tail lemma;\n"
        "  floor part = degree<2 terms (absolute strip, flux K0/K2, c0w V0, Phi_F-type constants): they sit in L only (checked: mode s, hi(L)<=2, and L-only in the lines above),\n"
        "  and enter the tail through the separate floor / L^1_loc argument (TAIL_LEMMA.tex), not through the scale-invariance condition.")
    log("TAIL_CERT " + ("PASS" if ok else "FAIL") + ": ALL STRUCTURAL CONDITIONS " + ("HOLD" if ok else "FAIL"))
    open("tail_cert_log.txt", "w").write("\n".join(out) + "\n")
    sys.exit(0 if ok else 1)
