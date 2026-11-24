#!/usr/bin/env python3
"""Task 2 / A1: H-J-W (STATEMENTS.tex Def. 3, Lemma 4) on the verifier result that feeds (A)'s J's.
Verifier file: step3nl/verifier/rate_table_fin11_a0.00015.json (the J's of base_ibp/lemmaS_F10p_M0.52.json args = its 'unc' block, rounded up;
asserted below).  Per march cell k: L_up_k = Ly_hi (barrier upper end on the cell), D_k.  The per-cell r1_k, Phi_k of closure.unconditional()
are not stored, so P_k = sum_{i<k} c_i D_i is bounded below by 0 (valid: c_i >= 0).
(W1) L_up_k <= Ly1 + 2 P_k : checked with P_k >= 0 only (sufficient, not necessary); fallback kappa_J from Lemma 4.
(W2) L_F <= L* and the tail barrier (verifier tail lemma, closure.tail(): rate_s <= -rho, alpha-condition) -- read from the stored tail block.
(W3) audit: closure.unconditional() runs the march with rhs_iv.set_mode(delta, dropL=True) (rhs_iv.py:198 multiplies L by 0) and
     Phi = unc.phi_val (degree-1 forcing only, unc.py:17-27); the degree-0 flux floor F0 = MU*FLUX0 + MU*EPSI*FLUX2 (rhs_iv.py:172,185) is in L,
     so it is NOT accounted on march cells (closure.py:182 adds (r1^+ Lup + Phi) D only; F0 enters only the tail slope X_s, closure.py:197).
     Repair: add F0 to c_k, i.e. + F0 * tau_t to X_c.
Result: J_m -> J_m (1 + kappa_J) for m = 0,1,2,4 (all homogeneous of degree 1 in X_c), slopes unchanged (X_s unchanged), with
     kappa_J = [ max_k (L_up_k - Ly1)^+ / 2 + F0 * sum_k D_k ] / X_c      (arb upper; X_c the stored float, which is what the J's divide).
Output: out/hjw_rep.json.   Run: PYTHONHASHSEED=0 nice -n 19 python3 -B rigor/hjw_rep.py [--rt PATH]"""
import os, sys, json, argparse, math
sys.dont_write_bytecode = True
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import arlib as L
from arlib import A, arb
RT_DEF = '/home/claude/s1spec/step3nl/verifier/rate_table_fin11_a0.00015.json'
if __name__ == '__main__':
    ap = argparse.ArgumentParser(); ap.add_argument('--rt', default=RT_DEF); a = ap.parse_args()
    rt = json.load(open(a.rt)); u = rt['unc']; tl = rt['tail']; B = json.load(open(L.BASE + '/lemmaS_F10p_M0.52.json'))['args']
    # the J's used by (A) are the verifier's, rounded up
    for m, nm in enumerate(('V0sq', 'V1sq', 'V2sq')):
        assert A(repr(B['J'][m])) >= A(repr(u[nm]['const'])) and A(repr(B['slopes'][m])) >= A(repr(u[nm]['slope_per_tau'])), nm
    assert A(repr(B['J4'])) >= A(repr(u['V4sq']['const'])) and A(repr(B['sJ4'])) >= A(repr(u['V4sq']['slope_per_tau']))
    X = lambda v: arb(repr(v))
    Ly1 = X(rt['Ly1']); Xc = X(u['X_const']); F0 = X(tl['floor_terms']['F0']); LF = X(tl['sqrt_LF'])**2; Ls = X(tl['L_star'])
    cells = []; Dsum = arb(0); worst = None; exc_max = arb(0)
    for st in rt['steps']:
        Lup = X(st['Ly_hi']); exc = Lup - Ly1; Dsum += X(st['D'])
        w1 = exc.upper() <= 0
        cells.append(dict(k=st['k'], D=st['D'], L_up=st['Ly_hi'], excess_over_Ly1=L.fup(exc) if exc.upper() > 0 else 0.0, W1_with_P_ge_0=w1))
        if exc.upper() > exc_max.upper(): exc_max = exc.upper(); worst = st['k']
    W1_ok = all(c['W1_with_P_ge_0'] for c in cells)
    W2_ok = bool((LF - Ls).upper() < 0 and tl['rho'] > 0 and tl['rate_s_U'] < 0)
    kJ_W1 = (exc_max/2)/Xc; kJ_W3 = F0*Dsum/Xc; kJ = kJ_W1 + kJ_W3
    res = dict(verifier_file=a.rt, verifier_sha256=L.sha(a.rt), Ly1=rt['Ly1'], X_const=u['X_const'], X_slope=u['X_slope'], F0=tl['floor_terms']['F0'],
               tau_march=L.fup(Dsum), Lup_max=max(c['L_up'] for c in cells), cells=cells,
               W1_ok=W1_ok, W1_note='checked with P_k >= 0 (per-cell r1_k, Phi_k not stored); fails at k = 0 already (L_up_0 = Ly_hi_0 > Ly1 by construction of the barrier), so Lemma 4 fallback',
               W2_ok=W2_ok, W2_note='L_F = %.3e <= L* = %.4e; tail barrier = verifier tail lemma (closure.tail, rate_s_U = %.5f <= -rho), not re-run here' % (float(LF.mid()), float(Ls.mid()), tl['rate_s_U']),
               W3_ok=False, W3_note='F0 (degree-0 flux floor, %.4e per unit tau) is NOT inside the march rate: unconditional() uses dropL=True (rhs_iv.py:198) and Phi has no constant term (unc.py:17-27); closure.py:182 omits F0*D. Repaired by +F0*tau_march in X_c (kappa_W3 below). Same omission in the tau_a = tau_1 J of v3.1 (relative 5e-7).' % tl['floor_terms']['F0'],
               kappa_J_W1=L.fup(kJ_W1), kappa_J_W3=L.fup(kJ_W3), kappa_J=L.fup(kJ), worst_cell=worst, label='cert',
               statement='int_{tau_a}^{tau_a+T} V_m^2 <= J_m (1 + kappa_J) + s_m T for all windows, m = 0,1,2,4 (J_m, s_m of lemmaS_F10p_M0.52 args)')
    print('H-J-W: W1 (P>=0) ok=%s worst cell %s; W2 ok=%s; W3 ok=False (F0 omitted on march, repaired); kappa_J = %.6g (W1 part %.6g, W3 part %.3g)' % (W1_ok, worst, W2_ok, L.fup(kJ), L.fup(kJ_W1), L.fup(kJ_W3)))
    L.dump(res, L.OUT + '/hjw_rep.json', [a.rt, L.BASE + '/lemmaS_F10p_M0.52.json', os.path.abspath(__file__), L.ROOT + '/rigor/arlib.py'])
