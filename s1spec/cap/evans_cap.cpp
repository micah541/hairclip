// Rigorous Evans function enclosure for the linearisation at S1 (computer-assisted part of step 2).
//
// State (16 real):  y = (p, q, k, m,  W01r,W01i, W02r,W02i, W03r,W03i, W12r,W12i, W13r,W13i, W23r,W23i)
// Reversed time tau = 8 - sigma.  Profile (reversed):  p' = k q, q' = -1 - k p, k' = -m, m' = -p/4.
// Compound (second exterior power) of  u'''' = a2 u'' + a1 u' + a0 u, reversed sign:
//   a0 = 1/4 - 3 m^2 - k p - lam,  a1 = -(5 k m + q/4),  a2 = -k^2,   lam = lr + i li  (interval parameters)
//   dW01 = W02, dW02 = W03 + W12, dW03 = a1 W01 + a2 W02 + W13, dW12 = W13,
//   dW13 = -a0 W01 + a2 W12 + W23, dW23 = -a0 W02 - a1 W12      (d/dsigma; reversed: minus all)
// Input file:  lr_lo lr_hi li_lo li_hi / tau1 / 16 lines x_lo x_hi / 256 numbers C (row-major) / 16 lines r0_lo r0_hi
//              / 16 lines rr_lo rr_hi
// Initial set = x_mid + C r0 + rr  (CAPD C0Rect2Set).  Output: W at the first return to q = 0 (tip).
#include <iostream>
#include <fstream>
#include <iomanip>
#include "capd/capdlib.h"
using namespace std;
using namespace capd;

int main(int argc, char* argv[]) {
  cout.precision(17);
  if (argc < 2) { cerr << "usage: evans_cap input [order]\n"; return 2; }
  int order = (argc > 2) ? atoi(argv[2]) : 20;
  const int n = 16;
  ifstream in(argv[1]);
  double a, b, c, d; in >> a >> b >> c >> d;
  interval lr(a, b), li(c, d);
  double tau1; in >> tau1;
  IVector x(n), r0(n), rr(n); IMatrix C(n, n);
  for (int i = 0; i < n; ++i) { double lo, hi; in >> lo >> hi; x[i] = interval(lo, hi); }
  for (int i = 0; i < n; ++i) for (int j = 0; j < n; ++j) { double v; in >> v; C[i][j] = interval(v); }
  for (int i = 0; i < n; ++i) { double lo, hi; in >> lo >> hi; r0[i] = interval(lo, hi); }
  for (int i = 0; i < n; ++i) { double lo, hi; in >> lo >> hi; rr[i] = interval(lo, hi); }
  if (!in) { cerr << "bad input\n"; return 2; }
  try {
    const char* f =
      "par:lr,li;"
      "var:p,q,k,m,a1r,a1i,a2r,a2i,a3r,a3i,b2r,b2i,b3r,b3i,c3r,c3i;"
      "fun:k*q, -1-k*p, -m, -0.25*p,"
      // W01' = W02
      "-a2r, -a2i,"
      // W02' = W03 + W12
      "-(a3r+b2r), -(a3i+b2i),"
      // W03' = a1 W01 + a2 W02 + W13   (a1 = -(5km+q/4), a2 = -k^2 real)
      "-(-(5*k*m+0.25*q)*a1r - k*k*a2r + b3r), -(-(5*k*m+0.25*q)*a1i - k*k*a2i + b3i),"
      // W12' = W13
      "-b3r, -b3i,"
      // W13' = -a0 W01 + a2 W12 + W23,  a0 = A - lam,  A = 1/4 - 3m^2 - kp
      "-(-(0.25-3*m*m-k*p)*a1r + (lr*a1r - li*a1i) - k*k*b2r + c3r),"
      "-(-(0.25-3*m*m-k*p)*a1i + (lr*a1i + li*a1r) - k*k*b2i + c3i),"
      // W23' = -a0 W02 - a1 W12
      "-(-(0.25-3*m*m-k*p)*a2r + (lr*a2r - li*a2i) + (5*k*m+0.25*q)*b2r),"
      "-(-(0.25-3*m*m-k*p)*a2i + (lr*a2i + li*a2r) + (5*k*m+0.25*q)*b2i);";
    IMap vf(f);
    vf.setParameter("lr", lr); vf.setParameter("li", li);
    IOdeSolver solver(vf, order);
    ITimeMap tm(solver);
    IVector xm(n), rx(n);
    for (int i = 0; i < n; ++i) { xm[i] = interval(x[i].mid().leftBound()); rx[i] = (x[i] - xm[i]) + rr[i]; }
    C0Rect2Set set(xm, C, r0, rx);
    IVector y1 = tm(tau1, set);
    bool qneg = (y1[1].rightBound() < 0.0);
    cout << "STAGE1 q(tau1) = " << y1[1] << (qneg ? "  (<0 ok)" : "  NOT NEGATIVE") << "\n";
    if (!qneg) { cout << "RESULT FAIL\n"; return 1; }
    ICoordinateSection section(n, 1);
    IPoincareMap pm(solver, section, poincare::MinusPlus);
    interval rt;
    IVector y2 = pm(set, rt);
    cout << "RETURN " << rt << "\n";
    const char* nm[16] = {"p","q","k","m","W01r","W01i","W02r","W02i","W03r","W03i","W12r","W12i","W13r","W13i","W23r","W23i"};
    for (int i = 0; i < n; ++i) cout << nm[i] << " " << y2[i].leftBound() << " " << y2[i].rightBound() << "\n";
    cout << "RESULT OK\n";
  } catch (std::exception& e) {
    cout << "RESULT EXCEPTION: " << e.what() << "\n";
    return 1;
  }
  return 0;
}
