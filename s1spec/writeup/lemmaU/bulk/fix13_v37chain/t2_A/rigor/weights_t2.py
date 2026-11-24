#!/usr/bin/env python3
"""Task 2 (r2), stage 1: quadratic coefficient bounds of (A) with the GLOBAL identities localised (Lemma L of STATEMENTS_r2.tex), certified (arb).
Replaces rigor/weights_rep.py (kept for reference).  New files only.

Audit (A5).  Phase 2a / v3.1 take the z_{j+1}-interpolation, the L4 identity and the {2,3,3} integration by parts SEPARATELY on the ramp [14,28]
(box (0.05, 0.9)) and on the core [28,350] (box (0.11, 0.52)), each with its own constants, and drop the boundary terms at y = 28 (and at 350):
  interpolation: Y_ramp[j+1] != Y_core[j+1] (1156.7 vs 644.2 at level 2), so  -Theta Y z_j z_{j+1} + 1/2 Theta' Y z_j^2  jumps at 28 (not covered by L(i));
  L4: the per-region Young/AM-GM constants differ (M1 = 0.9 vs 0.52), i.e. no single mu (L(iii) violated), and [Theta z_1 z_2^3]_{28} is dropped.
L(ii) with Ytilde falling from Y_ramp to Y_core on [28,56] (the r1 suggestion) is impossible: the core Ghat_22 at y = 28 becomes 243e3 > 153664.
Repair used here:
  (a) ramp coefficient box := core box (0.11, 0.52), z0 <= 1/50, z1 <= 0.52.  Valid on [14,28]: the verifier certifies |f| <= 0.0071, |y f_y| <= 0.4396,
      |z_0| <= 0.00116, |z_1| <= 0.4328 there (assembly 'ramp M0/M1' margins).  Then the ramp's quadratic, L4 and S33 class weights are the core's.
  (b) one C^{1,1} majorant Ytilde_j(s) of the z_{j+1}^2 coefficient on [14,700]:  Y_c on [14,40]; quintic rise to Y_s on [40,60]; Y_s on [60,175];
      quintic rise to Y_o on [175,350]; Y_o on O = [350,700].  Cost of a rise: Lemma L(ii) bracket  |Y'| kappa + |Y''|/2  (Xi = eps*omega there).
  (c) the L4 and S33 identities are localised by a partition of unity chi = cos^2/sin^2((pi/2) S(t)): piece P1a on [14,60] (ramp + core; reference
      region = core box), P1b on [40,350] (segments; reference = segment box), P2 on [175,700] (O; reference = O box), transitions [40,60], [175,350].
      In each piece ONE constant per product (mu for AM-GM, delta for Young against z_n), fixed at the reference region (L(iii)); other regions of the
      piece get mu*W_R^2 etc.  Transition: X_chi <= M1 (3||sqrt(chi) z3|| + ||(kappa sqrt(chi) + q) z2||), q = |chi'|/sqrt(chi) <= pi S'/ell.
      S33:  int Theta chi w |z2| z3^2 <= 1/2 X_chi ( ||(kappa sqrt(chi) + q) w z3|| + ||sqrt(chi) w z4|| ).
  (d) ramp foot terms (eta_in in Theta') as in v3.1 (cutoff-free J2 data, 'footw', 'W33'), but their AM-GM halves now act on all of P1a (core included).
Segment boxes: core [28,60] (0.11, 0.52), U <= 1/20, z0 <= 1/50, z1 <= 0.52 (B2'); (60,350] (M0s, M1s) default (0.15, 0.52), U <= 0.08 (U(350) < 0.0792),
  z0 <= 0.161 (B3U(c) |v| <= 0.16 plus delta U), z1 <= M1s + p*0.08 (61/100; 69/100 for the fallback (0.20, 0.60));
  O: (1/2, 1), U <= 0.162, z0 <= Z0 = 331/500, z1 <= Z1 = 1166/1000 (T1).
No lam, no eps4 (R5': no cut at 350).  On O a reserve delta_o a_o int Theta z_n^2 is left for the cutoff terms (outer_L2.py).
Checks (cert): Ytilde >= coefficient of z_{j+1}^2 on every cell; sup_zone eps*Ghat_jk <= eps(28)*Ghat^core_jk for all j,k on every zone/cell;
  dissipation budget on every zone; retained >= r_j a_c (absolute) on [28,350] at levels 0,1.
Output: out/weights_t2[_<tag>].json.   Run: PYTHONHASHSEED=0 nice -n 19 python3 -B rigor/weights_t2.py [--seg-box 15/100 52/100] [--tag fb]"""
import os, sys, json, argparse, math
sys.dont_write_bytecode = True
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import arlib as L
from arlib import WX, A, Q_, arb
import weights_rep as W0
R2Q = Q_(L.CFG.R2_EXPORT) if L.CFG.R2ON else Q_(0)   # v36 (ROUTE_DECISION.md G1/F2'): level-2 retained dissipation r_2 = 1/1000

KAP = 5; LN2 = arb(2).log(); PI = arb.pi()
P_Q = Q_(10246, 10000)
DFRAC = Q_(1, 10)                     # Young fraction per product against z_n at the reference region (as v3.1)
DEL = {3: Q_(1, 20), 2: Q_(1, 2)}     # ramp-foot AM-GM parameters (z_3 and z_2 parts)
DD = Q_(1, 100)                       # ramp S33 foot: dissipation fraction (reference) for its z_4 part
FOOT_RES = Q_(1, 100)                 # O reserve fraction (units a_o) used by the S33 outer-foot z_4 half
NC = 64                               # cells per transition interval

def S_(t): return t**3*(10 - 15*t + 6*t*t)
def S1_(t): return 30*t*t*(1 - t)**2
def S2_(t): return 60*t*(1 - t)*(1 - 2*t)

class Reg:
    def __init__(self, name, M0q, M1q, Uq, z0q, z1q):
        self.name, self.M0q, self.M1q, self.Uq, self.z0q, self.z1q = name, Q_(M0q), Q_(M1q), Q_(Uq), Q_(z0q), Q_(z1q)
        self.amin = (1 + A(self.M0q)**2)**-2; self.At = L.at_eval(self.M0q, self.M1q, self.Uq); self.Ab = L.ab_honest(self.At)
        self.M1 = A(self.z1q)
    def level(self, T, j):
        """jet-free YC, pointwise P items (with z_n Young at dfrac*amin of THIS region), L4 groups W[r], S33 weight"""
        n = j + 2; YC = dict(WX.YC_base(T[j], j, self.Ab)); used = arb(0); L4 = {}; S33 = arb(0)
        for k in range(n + 1): YC.setdefault(k, arb(0))
        for it in WX.classify_x(T, j, self.z1q, self.At, M0zq=self.z0q):
            if it[0] == 'P':
                x, y, wt = it[1], it[2], it[3]
                if x > y: x, y = y, x
                if y > n: raise ValueError(('P beyond n', j, it))
                if x == y == n: used += wt
                elif y == n: d = A(DFRAC)*self.amin; YC[x] += wt**2/(4*d); used += d
                elif x == y: YC[x] += wt
                else: YC[x] += wt/2; YC[y] += wt/2
            elif it[0] == 'L4': L4[it[1]] = L4.get(it[1], arb(0)) + it[2]
            else:
                if n < 4: raise ValueError('S33 at n < 4')
                S33 += it[1]
        return dict(YC=YC, used=used, L4=L4, S33=S33)

def piece_constants(refs, j):
    """per-product constants of a piece, fixed ONCE at the first region of refs where the group occurs (L(iii)); keys: ('L4', r, .) and ('S33', .)"""
    if isinstance(refs, dict): refs = [refs]
    n = j + 2; c = {}
    for ref in reversed(refs): c.update(_pc(ref, j))
    return c

def _pc(ref, j):
    """reference weights (a_ref, b_ref, amin_ref) of every product of the piece; mu = b_ref/a_ref (AM-GM) or delta = DFRAC amin_ref/(n-side weight)^2 (Young)"""
    n = j + 2; c = {}; M1 = ref['M1']; amin = ref['amin']
    for r, W in ref['L4'].items():
        c[('L4', r, 3)] = (W, M1, amin); c[('L4', r, 'k')] = (W, KAP*M1, amin); c[('foot', r)] = W
    if ref['S33'] > 0:
        w = ref['S33']; g = KAP*w
        c[('S33', 33)] = (M1, g, amin); c[('S33', 34)] = (M1, w, amin); c[('S33', 23)] = (KAP*M1, g, amin); c[('S33', 24)] = (KAP*M1, w, amin)
        c[('S33', 'foot3')] = g; c[('S33', 'foot4')] = A(DD)*amin/w**2
    return c

def alloc(data, C, j, chi=arb(1), q=arb(0), foot=None, YC=None, used=None):
    """add the allocation of one piece (reference weights C) evaluated with region data 'data' at partition value chi, q = |chi'|/sqrt(chi).
    Every product c ||a z_x|| ||b z_y|| of a localised identity: if x or y = n (dissipation variable) Young  delta ||b z_n||^2 + c^2 ||a z_x||^2/(4 delta)
    with delta = DFRAC amin_ref/b_ref^2; otherwise AM-GM  c/2 (mu ||a z_x||^2 + ||b z_y||^2/mu), mu = b_ref/a_ref (one constant per piece, L(iii)).
    foot = (M1r, etap): ramp-foot halves (act on the whole piece P1).  Returns (YC, used_abs, fw)."""
    n = j + 2; M1 = data['M1']; fw = arb(0)
    sc = chi.sqrt() if chi.lower() > 0 else arb(0).union(arb(chi.upper()).sqrt())     # chi in [0,1]; ball arithmetic may put its lower end below 0
    YC = YC if YC is not None else {k: arb(0) for k in range(n + 1)}; used = [used if used is not None else arb(0)]
    def prod(c, a, x, b, y, key):
        ar, br, am = C[key]
        if x == n: a, x, b, y, ar, br = b, y, a, x, br, ar
        if y == n:
            d = A(DFRAC)*am/(br*br); used[0] += d*b*b; YC[x] += c*c*a*a/(4*d)
        else:
            mu = br/ar; YC[x] += c*mu*a*a/2; YC[y] += c*b*b/(2*mu)
    for r, W in data['L4'].items():
        Wc = W*sc
        prod(3, Wc, r, M1*sc, 3, ('L4', r, 3)); prod(1, Wc, r, M1*(KAP*sc + q), 2, ('L4', r, 'k'))
        if foot is not None:                 # c ||W z_r|| sqrt(Jf2), c = 2 M1r etap:  <= c (DEL/W_ref ||W z_r||^2 + W_ref Jf2/(4 DEL))
            assert r < n
            M1r, etap = foot; cf = 2*M1r*etap; Wr = C[('foot', r)]; D_ = A(DEL[3] if r == 3 else DEL[2])
            YC[r] += cf*D_*Wc*Wc/Wr; fw += cf*Wr/(4*D_)
    if data['S33'] > 0:
        w = data['S33']; g = (KAP*sc + q)*w + sc*data.get('S33p', arb(0)); wz = sc*w; X3 = M1*sc; X2 = M1*(KAP*sc + q)
        prod(Q_(3, 2), X3, 3, g, 3, ('S33', 33)); prod(Q_(3, 2), X3, 3, wz, 4, ('S33', 34))
        prod(Q_(1, 2), X2, 2, g, 3, ('S33', 23)); prod(Q_(1, 2), X2, 2, wz, 4, ('S33', 24))
        if foot is not None:                 # 1/2 F (||g z3|| + ||w z4||), F = 2 M1r etap sqrt(Jf2)
            M1r, etap = foot; cf = M1r*etap; gr = C[('S33', 'foot3')]; D_ = A(DEL[3])
            YC[3] += cf*D_*g*g/gr; fw += cf*gr/(4*D_)
            dl = C[('S33', 'foot4')]; used[0] += dl*wz*wz; fw += cf*cf/(4*dl)
    return YC, used[0], fw

def ghat(YCt, j, Yt, Yp, Ypp, bc, amin):
    G = {k: YCt[k] for k in range(j + 1)}
    G[j] = YCt[j] + Yt*(Yt/(4*bc*amin) + A(Q_(KAP**2, 2))) + KAP*Yp + Ypp/2
    return G

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--seg-box', nargs=2, default=['15/100', '52/100'])
    ap.add_argument('--delta-o', default='1/4', help='reserve on O (units a_o) for the cutoff terms (outer_L2.py)')
    ap.add_argument('--retain', default='5/4'); ap.add_argument('--tag', default='')
    ap.add_argument('--y-rise1', nargs=2, default=['40', '50']); ap.add_argument('--y-chi1', nargs=2, default=['50', '60'])
    ap.add_argument('--y-rise2', nargs=2, default=['175', '250']); ap.add_argument('--y-chi2', nargs=2, default=['250', '350'])
    ap.add_argument('--bc-ramp', default='137/100', help='ramp interpolation share at level 2 (levels 0,1: 3/2)')
    a = ap.parse_args()
    T = WX.build_tables()
    M0s, M1s = Q_(a.seg_box[0]), Q_(a.seg_box[1]); z1s = Q_(math.ceil(float(M1s + P_Q*Q_(8, 100))*100), 100)
    assert z1s >= M1s + P_Q*Q_(8, 100)
    assert (L.U_of_y(350) - A(Q_(8, 100))).upper() < 0 and (L.U_of_y(60) - A(Q_(1, 20))).upper() < 0 and (L.U_of_y(700) - A(Q_(162, 1000))).upper() < 0
    RC = Reg('core', '11/100', '52/100', '1/20', '1/50', '52/100')
    RS = Reg('seg', M0s, M1s, '8/100', '161/1000', z1s)
    RO = Reg('O', '1/2', '1', '162/1000', '331/500', '1166/1000')
    assert Q_(331, 500) >= Q_(1, 2) + Q_(162, 1000) and Q_(1166, 1000) >= 1 + P_Q*Q_(162, 1000)
    amin_c, amin_s, amin_o = RC.amin, RS.amin, RO.amin
    amin_r = (1 + A(Q_(5, 100))**2)**-2
    etap = A(Q_(15, 8))/LN2; M1r = A(Q_(52, 100))
    y1a, y1b = (A(Q_(x)) for x in a.y_rise1); y2a, y2b = (A(Q_(x)) for x in a.y_rise2)
    c1a, c1b = (A(Q_(x)) for x in a.y_chi1); c2a, c2b = (A(Q_(x)) for x in a.y_chi2)
    assert Q_(a.y_rise1[1]) == Q_(a.y_chi1[0]) and Q_(a.y_chi1[1]) == 60 and Q_(a.y_rise2[1]) == Q_(a.y_chi2[0]) and Q_(a.y_chi2[1]) == 350 and Q_(a.y_rise2[0]) >= Q_(a.y_chi1[1])
    eps = lambda y: y**-4; e28 = eps(arb(28))
    rtgt = A(Q_(a.retain)); dlo = A(Q_(a.delta_o))
    global HF, HS
    HF, HS = arb(0), arb(0)                                                # sup over O of h_f^2 = 4 eps om eta_o'^2 and h_s^2 (cell cover, arb)
    for i in range(1200):
        t = (arb(i)/1200).union(arb(i + 1)/1200); y = 350*(t*LN2).exp(); eo = y**-4*28/y
        e1 = -S1_(t)/LN2; e2_ = -S2_(t)/LN2**2
        HF = max(HF, (4*eo*e1*e1).upper())
        # ((eta^2)'' - kappa (eta^2)')/eta = 2 eta'^2/eta + 2 eta'' - 2 kappa eta';  eta'^2/eta in closed form: eta'^2/eta = (S'/ln2)^2/(1-S) = 900 t^4 (1-t)/((1+3t+6t^2) ln2^2) (closed form, no division by 0)
        u2 = 2*900*t**4*(1 - t)/((1 + 3*t + 6*t*t)*LN2**2) + 2*e2_ - 2*KAP*e1
        HS = max(HS, (eo*u2*u2).upper())
    HF, HS = arb(HF), arb(HS)
    out_sups = dict(HF=HF, HS=HS)
    out = dict(levels={}, boxes=dict(core=('11/100', '52/100', '1/20', '1/50', '52/100'), seg=(str(M0s), str(M1s), '8/100', '161/1000', str(z1s)),
               O=('1/2', '1', '162/1000', '331/500', '1166/1000')), params=dict(DFRAC=str(DFRAC), DEL={str(k): str(v) for k, v in DEL.items()}, DD=str(DD),
               delta_o=a.delta_o, retain=a.retain, rise1=a.y_rise1, chi1=a.y_chi1, rise2=a.y_rise2, chi2=a.y_chi2, bc_ramp=a.bc_ramp))
    allok = True
    for j in (0, 1, 2):
        n = j + 2
        dC, dS, dO = RC.level(T, j), RS.level(T, j), RO.level(T, j)
        for d_, R_ in ((dC, RC), (dS, RS), (dO, RO)): d_['M1'] = R_.M1; d_['amin'] = R_.amin
        CC, CO = piece_constants([dC, dS], j), piece_constants(dO, j)      # P1 = [14,350]: constants of the core (L(iii)); P2 = O side
        # ---------------- core reference (P1a, chi = 1): defines Ghat_core, gamma'
        YCc, uc, fw = alloc(dC, CC, j, foot=(M1r, etap))
        YCc = {k: dC['YC'][k] + YCc[k] for k in range(n + 1)}; used_c = dC['used'] + uc          # absolute (per Theta z_n^2)
        if j < 2:
            fr = 2 - rtgt - used_c/amin_c; bc_c = Q_(math.floor(float(fr.lower())*100), 100)
        else:   # v36 (r_2): level 2 retains r_2 = 1/1000 (units amin_c) when V36_R2 = 1; v35: r_2 = 0
            if L.CFG.R2ON:   # 1e-6 cushion so that the certified retained ball lies above r_2 (on the 1/1000 grid it equals r_2 up to rounding)
                fr = 2 - A(R2Q) - used_c/amin_c - A(Q_(1, 10**6)); bc_c = Q_(math.floor(float(fr.lower())*10**6), 10**6)
            else:
                fr = 2 - used_c/amin_c; bc_c = Q_(math.floor(float(fr.lower())*1000), 1000)
        assert bc_c > 0, ('core Bc', j)
        ret_c = 2 - A(bc_c) - used_c/amin_c                       # units amin_c
        Yc = arb(float(YCc[j + 1].upper())*(1 + 1e-12)); assert (YCc[j + 1] - Yc).upper() < 0    # Ytilde on [14,40]: an exact point above the coefficient
        Gc = ghat(YCc, j, Yc, arb(0), arb(0), A(bc_c), amin_c)
        # ---------------- segments (P1b, chi = 1), box S
        YCs, us_, _ = alloc(dS, CC, j, foot=(M1r, etap)); YCs = {k: dS['YC'][k] + YCs[k] for k in range(n + 1)}; used_s = dS['used'] + us_
        # ---------------- O (P2, chi = 1)
        YCo, uo, _ = alloc(dO, CO, j); YCo = {k: dO['YC'][k] + YCo[k] for k in range(n + 1)}; used_o = dO['used'] + uo
        # outer foot terms of the P2 identities (cutoff terms of R5'; eta_o in Theta'):  X_2 <= Z1(3||z3|| + kappa||z2|| + ||h_f z2||_ds),
        #   h_f^2 = (eps om (eta_o^2)')^2/Theta = 4 eps om eta_o'^2;  S33 foot (1/6) int |(eps om (eta_o^2)')'| w |z2|^3 <= (w/6)||h_s z2||_ds X_2.
        #   AM-GM: base halves -> YC_o (O rate), data halves -> outer source on Nbar_2^2;  z4 halves -> Young with the O reserve.
        ofoot = arb(0); ores = arb(0); Z1o = RO.M1
        for r, Wr in dO['L4'].items():
            mu = (dO['YC'][r] + arb(1))/(Z1o*Wr*Wr)                       # base add ~ the jet-free coefficient (rate room on O is x2.5e4)
            YCo[r] += Z1o/2*mu*Wr*Wr; ofoot += Z1o/2*HF/mu
        if dO['S33'] > 0:
            w = dO['S33']; g = KAP*w
            mu = (dO['YC'][3] + arb(1))/(Z1o*g*g); YCo[3] += Z1o/4*mu*g*g; ofoot += Z1o/4*HF/mu      # 1/2 Z1 ||h_f z2|| ||g z3||
            dl = A(FOOT_RES)*amin_o/(w*w); c = Z1o/2; ores += dl*w*w; ofoot += c*c*HF/(4*dl)    # 1/2 Z1 ||h_f z2|| ||w z4||
            c3 = w/6*3*Z1o; mu3 = (dO['YC'][3] + arb(1))/c3; YCo[3] += c3/2*mu3; ofoot += c3/2*HS/mu3   # (w/6) 3 Z1 ||h_s z2|| ||z3||_Theta
            c2 = w/6*KAP*Z1o; mu2 = (dO['YC'][2] + arb(1))/c2; YCo[2] += c2/2*mu2; ofoot += c2/2*HS/mu2
            ofoot += w/6*Z1o*(HS*HF).sqrt()
        used_o = used_o + ores*0                                          # the foot z4 Young is charged to the reserve (ores), not to the base budget
        # ---------------- zones of [28,700] (Lemma L): Ytilde rises where only one piece lives, the partition switches where Ytilde is flat.
        def cellvals(t, ell, parts):
            """coefficient enclosure on a cell: parts = list of (data or data(t, ell), consts, foot, side) side 'one' (chi = 1), 'a' (cos^2), 'b' (sin^2)"""
            th = PI/2*S_(t); qv = PI*abs(S1_(t))/ell              # |chi'|/sqrt(chi) <= 2 theta' = pi S'(t)/ell (both pieces)
            YC = {k: arb(0) for k in range(n + 1)}; used = arb(0)
            for (dd_, C_, ft, side) in parts:
                if callable(dd_): dd_ = dd_(t, ell)
                if side == 'one': chi, q = arb(1), arb(0)
                else:
                    c_ = th.cos()**2 if side == 'a' else th.sin()**2
                    chi = arb(0).union(arb(c_.upper())) if c_.lower() < 0 else c_; q = qv
                YC, used, _ = alloc(dd_, C_, j, chi=chi, q=q, foot=ft, YC=YC, used=used)
                for k in range(n + 1): YC[k] = YC[k] + chi*dd_['YC'][k]
                used = used + chi*dd_['used']
            return YC, used
        def zone(ya, yb, parts, Ya, Yb, amin_z, rreq, rise):
            """cells of [ya,yb]: Ytilde = Ya + (Yb - Ya) S(t) if rise else Yb (constant); checks Ytilde >= coefficient of z_{j+1}^2, budget, eps*Ghat <= eps(28)*Ghat_core"""
            ell = (yb/ya).log(); worst = {}; okz = True; fails = []; supY = arb(0)
            for i in range(NC):
                lo, hi = arb(i)/NC, arb(i + 1)/NC; t = lo.union(hi)
                YCt, used = cellvals(t, ell, parts)
                if rise: Yt = Ya + (Yb - Ya)*S_(t); Ylo = Ya + (Yb - Ya)*S_(lo); Yp = abs(Yb - Ya)*abs(S1_(t))/ell; Ypp = abs(Yb - Ya)*abs(S2_(t))/ell**2
                else: Yt = Ylo = Yb; Yp = Ypp = arb(0)
                supY = max(supY, YCt[j + 1].upper())
                ok_y = (YCt[j + 1] - Ylo).upper() <= 0
                avail = 2 - used/amin_z - rreq
                if avail.lower() <= 0: okz = False; fails.append(dict(i=i, why='budget')); continue
                bc = A(Q_(math.floor(float(avail.lower())*10**6), 10**6))
                y = ya*(t*ell).exp(); G = ghat(YCt, j, Yt, Yp, Ypp, bc, amin_z)
                ok_g = all((eps(y)*G[k] - e28*Gc[k]).upper() <= 0 for k in range(j + 1))
                okz = okz and ok_y and ok_g
                rw = max(float((eps(y)*G[k]/(e28*Gc[k])).upper()) for k in range(j + 1))
                if not worst or rw > worst['ratio']: worst = dict(i=i, ratio=rw, G_jj=L.fup(G[j]), bc=float(bc))
                if not (ok_y and ok_g): fails.append(dict(i=i, ok_y=ok_y, ok_g=ok_g, YC=L.fup(YCt[j + 1]), Ylo=L.flo(Ylo)))
            return dict(ok=okz, worst=worst, failures=fails[:6], interval=[L.flo(ya), L.fup(yb)], supY=supY)
        rr = rtgt*amin_c                                           # absolute retained requirement (levels 0,1)
        rreq_c = (rr/amin_c) if j < 2 else A(R2Q); rreq_s = (rr/amin_s) if j < 2 else A(R2Q)*amin_c/amin_s   # v36: absolute r_2 amin_c on [28, 350] at level 2
        wc, ws = dC['S33'], dS['S33']
        def dW(t, ell):                                            # [50,60]: core box, S33 weight rises C^1 from the core to the segment value
            d_ = dict(dC); d_['S33'] = wc + (ws - wc)*S_(t) if ws > wc else wc; d_['S33p'] = abs(ws - wc)*abs(S1_(t))/ell if ws > wc else arb(0); return d_
        PA, PS, P2 = (dC, CC, (M1r, etap)), (dS, CC, (M1r, etap)), (dO, CO, None)
        zW1 = zone(c1a, c1b, [(dW, CC, (M1r, etap), 'one')], arb(0), arb(0), amin_c, rreq_c, False)
        zB2 = zone(c2a, c2b, [PS + ('a',), P2 + ('b',)], arb(0), arb(0), amin_s, rreq_s, False)
        Ys = YCs[j + 1]; Yo = YCo[j + 1]
        Yt1 = arb(float(max(Ys.upper(), zW1['supY'], Yc.upper()))*(1 + 1e-12)); Yt2 = arb(float(max(Yo.upper(), zB2['supY'], Yt1.upper()))*(1 + 1e-12))
        zA1 = zone(y1a, y1b, [PA + ('one',)], Yc, Yt1, amin_c, rreq_c, True)            # [40,50]: core box, Ytilde Yc -> Yt1
        zB1 = zone(c1a, c1b, [(dW, CC, (M1r, etap), 'one')], Yt1, Yt1, amin_c, rreq_c, False)   # [50,60]: S33 weight rises, Ytilde = Yt1
        zA2 = zone(y2a, y2b, [PS + ('one',)], Yt1, Yt2, amin_s, rreq_s, True)            # [175,250]: segment box, Ytilde Yt1 -> Yt2
        zB2 = zone(c2a, c2b, [PS + ('a',), P2 + ('b',)], Yt2, Yt2, amin_s, rreq_s, False)  # [250,350]: partition P1 | P2, Ytilde = Yt2
        z1 = dict(rise=zA1, partition=zB1, ok=zA1['ok'] and zB1['ok']); z2 = dict(rise=zA2, partition=zB2, ok=zA2['ok'] and zB2['ok'])
        # plateau (60,175]: P1b only, Ytilde = Yt1
        avail_s = 2 - used_s/amin_s - rreq_s; bc_s = Q_(math.floor(float(avail_s.lower())*10**6), 10**6); assert bc_s > 0, ('seg Bc', j)
        Gs = ghat(YCs, j, Yt1, arb(0), arb(0), A(bc_s), amin_s)
        ok_s = all((eps(c1b)*Gs[k] - e28*Gc[k]).upper() <= 0 for k in range(j + 1)) and (YCs[j + 1] - Yt1).upper() <= 0
        seg = dict(ok=ok_s, rows={f'{j}{k}': (L.fup(eps(c1b)*Gs[k]), L.flo(e28*Gc[k])) for k in range(j + 1)}, bc=str(bc_s),
                   retained_abs_over_amin_c=L.flo((2 - A(bc_s) - used_s/amin_s)*amin_s/amin_c), Ys=L.fup(Ys), Ytilde=L.fup(Yt1))
        # O: Bc_o = 2 - used_o/amin_o - delta_o, Ytilde = Yt2
        avail_o = 2 - used_o/amin_o - dlo; bc_o = Q_(math.floor(float(avail_o.lower())*10**6), 10**6); assert bc_o > 0, ('O Bc', j)
        Go = ghat(YCo, j, Yt2, arb(0), arb(0), A(bc_o), amin_o)
        ok_o = all((eps(arb(350))*Go[k] - e28*Gc[k]).upper() <= 0 for k in range(j + 1)) and (YCo[j + 1] - Yt2).upper() <= 0
        # ramp (P1a with core data): YC = core YC, Ytilde = Yc; Bc_r; slack for the R1/R2 Young users (units amin_r)
        bc_r = Q_(a.bc_ramp) if j == 2 else Q_(3, 2)
        slack_r = 2 - A(bc_r) - used_c/amin_r
        assert slack_r.lower() >= 0, ('ramp slack', j)
        Gr = ghat(YCc, j, Yc, arb(0), arb(0), A(bc_r), amin_r)
        ok_lvl = z1['ok'] and ok_s and z2['ok'] and ok_o
        allok = allok and ok_lvl
        out['levels'][str(j)] = dict(
            core=dict(YC=YCc, Ghat=Gc, Bc=str(bc_c), used=used_c/amin_c, retained=ret_c, L4=dC['L4'], S33=dC['S33'], footw=fw),
            ramp=dict(YC=YCc, Ghat=Gr, footw=fw, W33=dC['S33'], slack=slack_r, Bc=str(bc_r), note='coefficient box = core box (0.11,0.52); slack in units a_min(0.05)'),
            seg=dict(YC=YCs, Ghat=Gs, used=used_s/amin_s, L4=dS['L4'], S33=dS['S33'], check=seg),
            O=dict(YC=YCo, Ghat=Go, outer_foot_N2=ofoot, outer_foot_reserve=ores/amin_o, used=used_o/amin_o, Bc=str(bc_o), delta_o=a.delta_o, L4=dO['L4'], S33=dO['S33'], Yo=Yo, a_o=amin_o, Z1=RO.M1,
                   check=dict(ok=ok_o, rows={f'{j}{k}': (L.fup(eps(arb(350))*Go[k]), L.flo(e28*Gc[k])) for k in range(j + 1)})),
            transition1=z1, transition2=z2, Ytilde=dict(Yc=Yc, Ys=Ys, Yo=Yo, plateau1=Yt1, plateau_O=Yt2), ok=ok_lvl)
        print('level %d: core Ghat_jj %.2f (Bc %s, used %.4f, retained %.4f) | ramp Ghat_jj %.1f slack %.4f | seg ok %s (Ghat_jj %.1f, Bc %s) | t1 %s worst %.4g | t2 %s worst %.3g | O ok %s (Ghat_jj %.4g, Bc %s) | Y c/s/o %.1f/%.1f/%.1f'
              % (j, L.fup(Gc[j]), bc_c, L.fup(used_c/amin_c), L.flo(ret_c), L.fup(Gr[j]), L.flo(slack_r), ok_s, L.fup(Gs[j]), bc_s, z1['ok'], max(z1['rise']['worst'].get('ratio', -1), z1['partition']['worst'].get('ratio', -1)),
                 z2['ok'], max(z2['rise']['worst'].get('ratio', -1), z2['partition']['worst'].get('ratio', -1)), ok_o, L.fup(Go[j]), bc_o, L.fup(Yc), L.fup(Ys), L.fup(Yo)), flush=True)
        for nm, z in (('rise1', z1['rise']), ('chi1', z1['partition']), ('rise2', z2['rise']), ('chi2', z2['partition'])):
            if not z['ok']: print('   %s failures' % nm, z['failures'])
    r_exp = [Q_(math.floor(float(out['levels'][str(j)]['core']['retained'].lower())*1000 - 1e-9), 1000) for j in range(2)] + [R2Q]
    if L.CFG.R2ON: assert A(R2Q) <= out['levels']['2']['core']['retained'].lower()   # v36 (only with r_2 > 0; v35 retains exactly 0): the exported r_2 is below the certified core retained fraction
    for j in range(2):     # segments, transitions keep the same absolute retained amount (built into rreq); core: exported rational below the certified value
        assert A(r_exp[j]) <= out['levels'][str(j)]['core']['retained'].lower()
    out['r_export'] = [str(x) for x in r_exp]; out['all_ok'] = allok; out['outer_sups'] = out_sups
    out['windows'] = W0.windows(T, Q_(15, 100000), Q_(15, 100), honest=True)
    # compatibility keys (same shape as weights_rep.json) for ramp_rep.py and assembly_rep.py
    out['core'] = {j: dict(Ghat=out['levels'][j]['core']['Ghat'], YC=out['levels'][j]['core']['YC'], Bc=out['levels'][j]['core']['Bc']) for j in out['levels']}
    out['ramp'] = {j: {k: out['levels'][j]['ramp'][k] for k in ('Ghat', 'YC', 'footw', 'W33', 'slack')} for j in out['levels']}
    inputs = [L.BASE + '/rigor/weights_exact.py', L.BASE + '/rigor/tables_exact.pkl', os.path.abspath(__file__), L.ROOT + '/rigor/arlib.py', L.ROOT + '/rigor/weights_rep.py']
    path = L.OUT + '/weights_t2%s.json' % (('_' + a.tag) if a.tag else '')
    L.dump(out, path, inputs); print('all checks:', allok, ' r_export', out['r_export'], ' wrote', path)

if __name__ == '__main__':
    main()
