#!/usr/bin/env python3
"""Task 4, U3 + U4: certified (arb) sweeps of Prop. U and the launch constants, weight Psi = eta^2 chi^2 W(sigma) D phi.

  PYTHONHASHSEED=0 nice -n 19 python3 -B u34_sweeps.py [--u1 out/u1_ghat.json] [--case primary] [--h 0.005] [--ny 400]  -> out/u34_<case>.json

W(sigma) = 1/(e^{-4(sigma+a)} + 1 + e^{4(sigma-l-a)}), a = 3/10, l = ln 1.5;  D = 1 + g y^-4 (g of the case, u1);  phi = exp(K_b/y), K_b = 1e-3.
Facts used [paper]: W'/W = 4(p - q), p - q has the sign of l/2 - sigma, so W increases on sigma <= l/2 and decreases on sigma >= l/2;
  for y <= 56 <= y_c and x >= 0 the argument ln(y/y_c) - x <= 0 < l/2, so every sweep integrand is nonincreasing in x (hence the upper Riemann sum at
  the left grid points is an upper bound) and nonincreasing in y_c (hence the sup over a y_c-cell [a, b] is the value at a).
  Tail x > X: W(sigma) <= e^{4(sigma + a)}, so int_X^inf sup_y f(y) W(ln(y/y_c) - x) dx <= sup_y [f(y) (y/y_c)^4] e^{4a} e^{-4X}/4.
Sweeps (U_ramp.tex Sec. 3, STATEMENTS_r2 Def. 7), for x = c^-(tau - tau_rho), c^- = 1/4 - abar:
  Sigma_in = (1/c^-) int_0^inf sup_{[28,56]} (eta^2)' D phi (y/28) W dx   (Prop. U uses c^+ Sigma_in, c^+ = 1/4 + abar + bbar/28)
  Sigma_E  = (1/c^-) int sup (28/y)^3 D phi W dx ;  S_X = sup W(ln(y/y_c)) D phi y/28 ;  Sigma_X = (1/c^-) int sup W D phi y/28 dx
  Sigma_rs^(j) = (1/c^-) int sup eps [Ghat_jj D - (1 - 4abar) g]_+ eta^2 phi (y/28) W dx       (zero for j = 0, 1 under (R<2): asserted)
  For each y-cell: f-factors enclosed on the cell (Taylor forms of ulib for eta), W bounded by its value at the upper end of the sigma-range (monotone).
Launch constants (launch_lemma.tex; STATEMENTS_r2 Remark 2):
  C_A^D = sup_{[28,350]} W(ln(y/56)) D phi y/28 ;  w_min = min_{[0,l]} W = W(0) = W(l) ;
  q = sum_{i>=0} W(sigma_0 + i l)/w_min, sigma_0 = ln(350/56) (W decreasing there; a window [rho', 1.5 rho'] of a launched neighbour has Psi >= w_min times the
      common factors eta^2 chi^2 D phi; the window ending at r_chi*1.5 lambda_1 is dominated by the next cell term), q^D = D(350) phi(350) q (as stated; >= q);
  Delta_in_j = p^{2j} c_dU^2 int_{ln28}^{ln350} W(s - ln56) D phi e^{2ps} ds ;
  c_2(y_c) = sup_{[28,60]} W(ln(y/y_c)) D phi y/28 ;  c_3(y_c) = sum_{k=0..4} sup_{[y_k, 1.5 y_k]} W(ln(y/y_c)) (1.5 y_k/28) D(y_k) phi(y_k), y_k = 50 * 1.5^k
      (present-centre launch data: E^rho_j(tau_1) <= c_2 H2_j + c_3 H3_j; on a y_c-cell [a, b]: c_2 at a, c_3 with sigma over [ln(y_k/b), ln(1.5 y_k/a)]);
  c_dU = 1.3056 (betabar + DL alpha_1 + p abar alpha_1)/(p (1/4 - abar)), DL = |lambda~ + p/4| (fix4_r3 R3; was without DL, i.e. lambda~ = -p/4)
      (launch_lemma.tex Lemma 3 with alpha_tau = lambda~ alpha + beta; betabar, abar from the verifier);
  Phi_j(y_c) = (1/c^-) [ bbar p U_max/y_c (D(28) phi(28) int W e^{-2 sigma})^{1/2} + c'_j U_max e^{2a} (D(28) phi(28))^{1/2}/(4 * 28^4) (28/y_c)^2 ],
      c'_j = (0.0520, 0.1446, 0.4324) (Prop. B' of v3.2, step3_proof_v32.tex l.245) [paper input];
  window Sobolev on an interval of s-length l: sup u^2 <= |u|^2/l + 2 |u| |u'| [paper]."""
import sys, os, json, math, time, argparse
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import t4lib as L
sys.path.insert(0, L.P2)
from flint import arb, ctx
import ulib as U
sys.path.insert(0, '/home/claude/s1spec/writeup/lemmaU/bulk/fix13_v37chain/r3'); import r3_forcing as R3   # fix4_r3 (R3, v36pre): DL = |lambda~ + p/4|; v36: gated by V36_R3
from t4lib import A, Q_, fup, flo, num
ctx.prec = 128
ap = argparse.ArgumentParser(); ap.add_argument('--u1', default=os.path.join(L.OUTD, 'u1_ghat.json')); ap.add_argument('--case', default='primary')
ap.add_argument('--h', default='0.005'); ap.add_argument('--X', default='8'); ap.add_argument('--ny', type=int, default=400); args = ap.parse_args()
T0 = time.time()
U1 = json.load(open(args.u1)); CASE = U1['cases'][args.case]
gq = Q_(CASE['design']['g']); g = A(gq)
GJJ = [arb(CASE['Ghat_jj'][j]['double']) for j in range(3)]
VS = L.verifier_sups(); ABAR, BBAR, BETAB = VS['abar'], VS['bbar'], VS['betabar']
CM = arb(1)/4 - ABAR; CP = arb(1)/4 + ABAR + BBAR/28
AW = A(Q_(3, 10)); LL = (arb(3)/2).log(); KB = A(Q_(1, 1000)); P = A(Q_(10246, 10000)); ALPHA1 = A(Q_(15, 100000)); UMAX = A(Q_(9, 100))
LN28, LN56, LN2 = arb(28).log(), arb(56).log(), arb(2).log()
def Wf(sig): return 1/((-4*(sig + AW)).exp() + 1 + (4*(sig - LL - AW)).exp())
def Dy(y): return 1 + g/(y**4)
def phy(y): return (KB/y).exp()
def up(x): return arb(x.upper())

# ---------------------------------------------------------------- y-cells on the ramp [28, 56]
NY = args.ny; CELLS = []
for i in range(NY):
    tlo, thi = arb(i)/NY, arb(i + 1)/NY; t = tlo.union(thi)
    s = LN28 + t*LN2; y = s.exp(); s_hi = LN28 + thi*LN2
    e2v, de2 = U.e2(t, 0), U.e2(t, 1)
    Dp = Dy(y)*phy(y); eps = 1/(y**4)
    br = GJJ[2]*Dy(y) - (1 - 4*ABAR)*g
    brp = arb(0) if br.upper() <= 0 else (br if br.lower() >= 0 else arb(0).union(br))
    for k in (0, 1):
        assert (GJJ[k]*Dy(arb(28)) - (1 - 4*ABAR)*g).upper() < 0           # (R<2): Sigma_rs^(0,1) = 0
    f = dict(inn=up(abs(de2)*Dp*y/28), E=up((28/y)**3*Dp), X=up(Dp*y/28), rs=up(eps*brp*e2v*phy(y)*y/28))
    CELLS.append(dict(s_hi=arb(s_hi.upper()), y_hi=up(y), f=f))
H = A(Q_(args.h)); XMAX = A(Q_(args.X)); NX = int(Q_(args.X)/Q_(args.h))

def sweeps(yc):
    lyc = arb(yc).log() if not isinstance(yc, arb) else yc.log()
    tot = {k: arb(0) for k in ('inn', 'E', 'X', 'rs')}; SX = arb(0)
    for ix in range(NX):
        x = H*ix; sup = {k: arb(0) for k in tot}
        for c in CELLS:
            w = up(Wf(c['s_hi'] - lyc - x))
            for k in tot:
                v = c['f'][k]*w
                if v.upper() > sup[k]: sup[k] = up(v)
        if ix == 0: SX = sup['X']
        for k in tot: tot[k] += H*sup[k]
    tail = {}
    for k in tot:
        m = arb(0)
        for c in CELLS:
            v = c['f'][k]*(c['y_hi']/arb(yc))**4
            if v.upper() > m: m = up(v)
        tail[k] = m*(4*AW).exp()*(-4*XMAX).exp()/4
        tot[k] = (tot[k] + tail[k])/CM
    return dict(Sigma_in=tot['inn'], Sigma_E=tot['E'], Sigma_X=tot['X'], Sigma_rs2=tot['rs'], S_X=SX, tail_rel=max(fup(tail[k]/CM/tot[k]) for k in tot if tot[k] > 0))

# ---------------------------------------------------------------- launch constants
def sup_cells(fun, s0, s1, n):
    best = arb(0)
    for i in range(n):
        s = (s0 + (s1 - s0)*arb(i)/n).union(s0 + (s1 - s0)*arb(i + 1)/n)
        v = fun(s)
        if v.upper() > best: best = up(v)
    return best
def int_cells(fun, s0, s1, n):
    tot = arb(0)
    for i in range(n):
        a_, b_ = s0 + (s1 - s0)*arb(i)/n, s0 + (s1 - s0)*arb(i + 1)/n
        tot += up(fun(a_.union(b_)))*(b_ - a_)
    return tot
wmin = Wf(arb(0)); wl = Wf(LL); assert abs(wmin - wl).upper() < 1e-30
LN350 = arb(350).log()
CAD = sup_cells(lambda s: Wf(s - LN56)*Dy(s.exp())*phy(s.exp())*s.exp()/28, LN28, LN350, 4000)
sig0 = LN350 - LN56
qsum = sum((Wf(sig0 + i*LL) for i in range(60)), arb(0)) + (-4*(sig0 + 60*LL - LL - AW)).exp()/(1 - (-4*LL).exp())
assert (sig0 - LL/2).lower() > 0
q = qsum/wmin; qD = Dy(arb(350))*phy(arb(350))*q
cdU = A(Q_(13056, 10000))*(BETAB + (A(R3.DL_Q)*ALPHA1 if L.CFG.R3ON else 0) + P*ABAR*ALPHA1)/(P*CM)   # R3: betabar' = betabar + |lambda~ + p/4| alpha_1 (gated)
assert cdU.upper() < A(Q_(2, 10**8))                                   # used as a bound in u1 (M1z check) and by task 5 ([hyp] there)
IDin = int_cells(lambda s: Wf(s - LN56)*Dy(s.exp())*phy(s.exp())*(2*P*s).exp(), LN28, LN350, 4000)
Din = [P**(2*j)*cdU**2*IDin for j in range(3)]
IWe = int_cells(lambda s: Wf(s)*(-2*s).exp(), arb(-12), arb(12), 4000) + (4*AW).exp()*arb(-24).exp()/2 + (4*(LL + AW)).exp()*arb(-72).exp()/6
CPR = [A(Q_('0.0520')), A(Q_('0.1446')), A(Q_('0.4324'))]
D28p = Dy(arb(28))*phy(arb(28))
def Phi(yc):
    yc = arb(yc)
    return [(BBAR*P*UMAX/yc*(D28p*IWe).sqrt() + CPR[j]*UMAX*(2*AW).exp()*D28p.sqrt()/(4*arb(28)**4)*(28/yc)**2)/CM for j in range(3)]
def Wsup_interval(s0, s1):
    """sup of W on [s0, s1] (arb points), by unimodality (max at l/2)"""
    if s1.upper() < (LL/2).lower(): return up(Wf(arb(s1.upper())))
    if s0.lower() > (LL/2).upper(): return up(Wf(arb(s0.lower())))
    return up(Wf(LL/2))
YK = [A(50)*(arb(3)/2)**k for k in range(5)]
def c23(a_, b_):
    """present-centre launch coefficients on the y_c-cell [a_, b_]"""
    la, lb = arb(a_).log(), arb(b_).log()
    c2 = sup_cells(lambda s: Wf(s - la)*Dy(s.exp())*phy(s.exp())*s.exp()/28, LN28, arb(60).log(), 800)
    assert (arb(60).log() - la - LL/2).upper() < 0                     # W increasing on the [28, 60] range: c_2 nonincreasing in y_c
    c3 = arb(0)
    for yk in YK:
        c3 += Wsup_interval(yk.log() - lb, (3*yk/2).log() - la)*(3*yk/2/28)*Dy(yk)*phy(yk)
    return c2, c3
print(f'case {args.case}: g {gq}, w_min {fup(wmin):.10f}, C_A^D {fup(CAD):.6f}, q {fup(q):.6f}, q^D {fup(qD):.6f}, c_dU {fup(cdU):.6e}, '
      f'Delta_in {[fup(x) for x in Din]}, int W e^-2s {fup(IWe):.6f}', flush=True)

YCG = [Q_(56), Q_(60), Q_(64), Q_(70), Q_(76), Q_(84), Q_(92), Q_(100), Q_(110), Q_(120), Q_(130), Q_(140), Q_(150), Q_(159), Q_(16867, 100), Q_(185), Q_(200), Q_(220), Q_(253)]
SW = {}
for yc in YCG:
    t1 = time.time(); r = sweeps(A(yc)); SW[str(yc)] = r
    print(f'  y_c {float(yc):7.2f}: Sigma_in {fup(r["Sigma_in"]):.6f} Sigma_E {fup(r["Sigma_E"]):.6f} S_X {fup(r["S_X"]):.6f} Sigma_X {fup(r["Sigma_X"]):.6f} '
          f'Sigma_rs2 {fup(r["Sigma_rs2"]):.4e}  ({time.time() - t1:.1f} s)', flush=True)
CELLS_P = [(YCG[i], YCG[i + 1]) for i in range(len(YCG) - 1)]
PRES = []
for a_, b_ in CELLS_P:
    c2, c3 = c23(A(a_), A(b_)); PRES.append(dict(yc=[str(a_), str(b_)], c2=num(c2), c3=num(c3), Phi=[num(x) for x in Phi(A(a_))],
                                                Ubar=num(A(Q_(13056, 10000))*ALPHA1*((3*A(b_)/2).log()*P).exp()), sweeps_at=str(a_)))
mono = all(fup(SW[str(YCG[i + 1])][k]) <= fup(SW[str(YCG[i])][k]) for i in range(len(YCG) - 1) for k in ('Sigma_in', 'Sigma_E', 'Sigma_X', 'S_X'))
out = dict(case=args.case, g=str(gq), weight=dict(mu=4, nu=4, a_W='3/10', l='ln 1.5', K_b='1/1000'),
           grid=dict(h=args.h, X=args.X, ny=NY), c_minus=num(CM), c_plus=num(CP),
           sweeps={str(k): {kk: (num(v) if isinstance(v, arb) else v) for kk, v in r.items()} for k, r in SW.items()},
           sweeps_monotone_check=mono,
           launch=dict(CA_D=num(CAD), q=num(q), qD=num(qD), wmin=num(wmin, d='lo'), c_dU=num(cdU), Delta_in=[num(x) for x in Din], int_W_em2s=num(IWe),
                       Phi_late=[num(x) for x in Phi(arb(56))], D350phi=num(Dy(arb(350))*phy(arb(350))),
                       sobolev=dict(rule='sup u^2 <= |u|^2/l + 2 |u| |u\'| on [rho, 1.5 rho] (s-length l = ln1.5); |u|^2 <= E/w_min', inv_l=num(1/LL), two=2, wmin=num(wmin, d='lo'))),
           present_cells=PRES)
out['args'] = L.args_block([args.u1, L.VERIF_ALPHA, os.path.join(L.P2, 'ulib.py'), os.path.abspath(__file__), os.path.join(HERE, 't4lib.py'), os.path.abspath(R3.__file__)])
L.dump(out, os.path.join(L.OUTD, f'u34_{args.case}.json'))
print('monotone in y_c (grid check):', mono, ' wall %.1f s' % (time.time() - T0))
