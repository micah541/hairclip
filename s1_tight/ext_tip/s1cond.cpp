// Tight S1 profile from the step-1 shooting condition (conditioned enclosure).
// The S1 state at sigma = 8 is a point y* of the certified set Y = xm + C r0 + rr (in_br2.txt, as s1prof) with
// M(y*) = 0, where M = m at the first return to q = 0 after tau1 (MinusPlus): the function whose zero defines
// eps* in step 1 (s1proof.cpp / bracket2.py).  Write y* = xm + C d + e, d in r0, e in rr.  By the mean value
// theorem on the segment [xm, y*] (inside the convex set Y), with A = DM(eta) C, G = DPhi_t(xi) C:
//    0 = M(xm) + A d + DM(eta) e  =>  d_0 = -(M(xm) + sum_{j>=1} A_j d_j + DM(eta) e)/A_0,
//    Phi_t(y*) = Phi_t(xm) - (G_0/A_0)(M(xm) + DM e) + sum_{j>=1} (G_j - G_0 A_j/A_0) d_j + DPhi_t e.
// DM(eta) and DPhi_t(xi) are enclosed over all of Y by CAPD C1 computations (C1 Poincare map + computeDP,
// C1 solution curve).  Phi_t(xm), M(xm) come from point runs.  This holds for EVERY zero of M in Y (no
// uniqueness is used).  The result is intersected with the plain box enclosure (s1prof), which is also valid.
// Reversed time tau = 8 - sigma, field (p,q,k,m,th)' = (k q, -1 - k p, -m, -p/4, -k), as s1prof.cpp.
// Usage: s1cond input Tend N [order [tipfile]]
// stdout: N rows in the prof_tight format (ta tb tc, centre lo/hi x5, hull lo/hi x5), same grid as s1prof;
//         centre = conditioned enclosure intersected with the box centre, hull = box hull (= prof_tight hull).
// tipfile: conditioned tip data (return time, tip state), A, and diagnostics.
#include <iostream>
#include <fstream>
#include "capd/capdlib.h"
using namespace std;
using namespace capd;
static const char* VF = "var:p,q,k,m,th;fun:k*q,-1-k*p,-m,-0.25*p,-k;";
static interval meet(const interval& a, const interval& b, bool& ok) {
  double lo = std::max<double>(a.leftBound(), b.leftBound()), hi = std::min<double>(a.rightBound(), b.rightBound());
  if (lo > hi) { ok = false; return a; }
  return interval(lo, hi);
}
int main(int argc, char* argv[]) {
  cout.precision(17);
  if (argc < 4) { cerr << "usage: s1cond input Tend N [order [tipfile]]\n"; return 2; }
  double Tend = atof(argv[2]); int N = atoi(argv[3]);
  int order = (argc > 4) ? atoi(argv[4]) : 25;
  ifstream in(argv[1]);
  double tau1; in >> tau1;
  IVector x(5), r0(5); IMatrix C(5, 5);
  for (int i = 0; i < 5; ++i) { double lo, hi; in >> lo >> hi; x[i] = interval(lo, hi); }
  for (int i = 0; i < 5; ++i) for (int j = 0; j < 5; ++j) { double c; in >> c; C[i][j] = interval(c); }
  for (int i = 0; i < 5; ++i) { double lo, hi; in >> lo >> hi; r0[i] = interval(lo, hi); }
  if (!in) { cerr << "bad input\n"; return 2; }
  try {
    IMap vf(VF);
    IVector xm(5), rr(5);
    for (int i = 0; i < 5; ++i) { xm[i] = interval(x[i].mid().leftBound()); rr[i] = x[i] - xm[i]; }
    for (int i = 0; i < 5; ++i)
      if (!(r0[i].contains(0.0) && rr[i].contains(0.0))) { cerr << "xm not in Y\n"; return 1; }
    ICoordinateSection section(5, 1);   // q = 0
    // ---- (a) C1 Poincare map over Y: enclosures of DP (all y in Y) and of the return-time derivative dT
    IMatrix DP(5, 5); IVector dT(5), PY(5); interval rtY;
    {
      IOdeSolver solver(vf, order); ITimeMap tm(solver);
      C1Rect2Set s(xm, C, r0, rr);
      IVector y1 = tm(tau1, s);
      if (!(y1[1].rightBound() < 0.0)) { cerr << "C1: q(tau1) not negative\n"; return 1; }
      IPoincareMap pm(solver, section, poincare::MinusPlus);
      IMatrix dF(5, 5);
      PY = pm(s, dF, rtY);
      DP = pm.computeDP(PY, dF, dT, rtY);
    }
    IVector A(5); interval DMe(0.0), dTe(0.0); IVector dTC(5);
    for (int j = 0; j < 5; ++j) {
      A[j] = 0.0; dTC[j] = 0.0;
      for (int i = 0; i < 5; ++i) { A[j] += DP[3][i] * C[i][j]; dTC[j] += dT[i] * C[i][j]; }
    }
    for (int i = 0; i < 5; ++i) { DMe += DP[3][i] * rr[i]; dTe += dT[i] * rr[i]; }
    if (!(A[0].rightBound() < 0.0)) { cerr << "A0 not negative\n"; return 1; }
    // ---- (b) point runs from xm: M(xm), return time, tip state; trajectory
    IVector P0(5); interval T0;
    {
      IOdeSolver solver(vf, order); ITimeMap tm(solver);
      C0Rect2Set s(xm);
      IVector y1 = tm(tau1, s);
      if (!(y1[1].rightBound() < 0.0)) { cerr << "point: q(tau1) not negative\n"; return 1; }
      IPoincareMap pm(solver, section, poincare::MinusPlus);
      P0 = pm(s, T0);
    }
    interval M0 = P0[3];
    IOdeSolver solP(vf, order); ITimeMap tmP(solP);
    C0Rect2Set sP(xm); ITimeMap::SolutionCurve psol(0.0); tmP(Tend, sP, psol);
    // ---- (c) C1 solution curve over Y
    IOdeSolver sol1(vf, order); ITimeMap tm1(sol1);
    C1Rect2Set s1(xm, C, r0, rr); ITimeMap::SolutionCurve c1sol(0.0); tm1(Tend, s1, c1sol);
    // ---- (d) plain box run, exactly as s1prof.cpp
    IOdeSolver solB(vf, order); ITimeMap tmB(solB);
    C0Rect2Set sB(xm, C, r0, rr); ITimeMap::SolutionCurve bsol(0.0); tmB(Tend, sB, bsol);
    // ---- conditioned tip data
    IVector Ptip(5); IMatrix DPC = DP * C; interval Ttip;
    for (int c = 0; c < 5; ++c) {
      interval g0 = DPC[c][0];
      interval v = P0[c] - (g0 / A[0]) * (M0 + DMe);
      for (int j = 1; j < 5; ++j) v += (DPC[c][j] - g0 * A[j] / A[0]) * r0[j];
      for (int k = 0; k < 5; ++k) v += DP[c][k] * rr[k];
      Ptip[c] = v;
    }
    {
      interval v = T0 - (dTC[0] / A[0]) * (M0 + DMe);
      for (int j = 1; j < 5; ++j) v += (dTC[j] - dTC[0] * A[j] / A[0]) * r0[j];
      Ttip = v + dTe;
    }
    bool ok = true; double wmax[5] = {0, 0, 0, 0, 0}, wbox[5] = {0, 0, 0, 0, 0};
    for (int i = 0; i < N; ++i) {
      double ta = Tend * i / N, tb = (i == N - 1) ? Tend : Tend * (i + 1) / N;
      double tc = 0.5 * (ta + tb);
      IVector Fc = psol(interval(tc));
      IMatrix DF = c1sol.derivative(interval(tc));
      IVector bc = bsol(interval(tc));
      IVector vh = bsol(interval(ta, tb));
      IMatrix G = DF * C;
      cout << ta << " " << tb << " " << tc;
      for (int c = 0; c < 5; ++c) {
        interval g0 = G[c][0];
        interval v = Fc[c] - (g0 / A[0]) * (M0 + DMe);
        for (int j = 1; j < 5; ++j) v += (G[c][j] - g0 * A[j] / A[0]) * r0[j];
        for (int k = 0; k < 5; ++k) v += DF[c][k] * rr[k];
        bool okc = true; interval cen = meet(v, bc[c], okc);
        if (!okc) { ok = false; cerr << "EMPTY intersection row " << i << " comp " << c << "\n"; }
        wmax[c] = std::max<double>(wmax[c], cen.rightBound() - cen.leftBound());
        wbox[c] = std::max<double>(wbox[c], bc[c].rightBound() - bc[c].leftBound());
        cout << " " << cen.leftBound() << " " << cen.rightBound();
      }
      for (int c = 0; c < 5; ++c) cout << " " << vh[c].leftBound() << " " << vh[c].rightBound();
      cout << "\n";
    }
    cout << "END\n";
    if (argc > 5) {
      ofstream tf(argv[5]); tf.precision(17);
      tf << "tau_tip " << Ttip.leftBound() << " " << Ttip.rightBound() << "\n";
      const char* nm[5] = {"p", "q", "k", "m", "th"};
      for (int c = 0; c < 5; ++c) tf << "tip_" << nm[c] << " " << Ptip[c].leftBound() << " " << Ptip[c].rightBound() << "\n";
      tf << "box_tau_tip " << rtY.leftBound() << " " << rtY.rightBound() << "\n";
      for (int c = 0; c < 5; ++c) tf << "box_tip_" << nm[c] << " " << PY[c].leftBound() << " " << PY[c].rightBound() << "\n";
      tf << "M_xm " << M0.leftBound() << " " << M0.rightBound() << "\n";
      tf << "T_xm " << T0.leftBound() << " " << T0.rightBound() << "\n";
      for (int j = 0; j < 5; ++j) tf << "A" << j << " " << A[j].leftBound() << " " << A[j].rightBound() << "\n";
      for (int j = 0; j < 5; ++j) tf << "dTC" << j << " " << dTC[j].leftBound() << " " << dTC[j].rightBound() << "\n";
      tf << "DMe " << DMe.leftBound() << " " << DMe.rightBound() << "\n";
      for (int c = 0; c < 5; ++c) tf << "maxwidth_" << nm[c] << " cond " << wmax[c] << " box " << wbox[c] << "\n";
      tf << "intersections_nonempty " << (ok ? 1 : 0) << "\n";
    }
    if (!ok) return 1;
  } catch (std::exception& e) { cout << "EXCEPTION: " << e.what() << "\n"; return 1; }
  return 0;
}
