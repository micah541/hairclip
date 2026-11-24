// Rigorous backward integration for the CDF shrinker S1 (computer-assisted part).
//
// Reversed-time vector field (tau = S - sigma) for y = (p, q, k, m, th):
//     p' = k q,  q' = -1 - k p,  k' = -m,  m' = -p/4,  th' = -k.
//
// Input file (all numbers are lower/upper bounds, already outward-rounded to doubles):
//   tau1
//   5 lines: x_i_lo x_i_hi            (center, as an interval vector)
//   25 numbers: C (row-major), exact doubles
//   5 lines: r0_i_lo r0_i_hi
// The initial set is  x + C*r0  (CAPD C0Rect2Set / Lohner doubleton).
//
// Output: enclosure of the state at tau1 (to check q<0), then the first return to q=0
// with crossing direction MinusPlus (in tau-time), the return time, and the enclosure of
// the crossing point.  Any failure of CAPD's validation throws and the program reports it.
#include <iostream>
#include <fstream>
#include <iomanip>
#include "capd/capdlib.h"
using namespace std;
using namespace capd;

int main(int argc, char* argv[]) {
  cout.precision(17);
  if (argc < 2) { cerr << "usage: s1proof input [order]\n"; return 2; }
  int order = (argc > 2) ? atoi(argv[2]) : 25;
  ifstream in(argv[1]);
  double tau1; in >> tau1;
  IVector x(5), r0(5); IMatrix C(5, 5);
  for (int i = 0; i < 5; ++i) { double lo, hi; in >> lo >> hi; x[i] = interval(lo, hi); }
  for (int i = 0; i < 5; ++i) for (int j = 0; j < 5; ++j) { double c; in >> c; C[i][j] = interval(c); }
  for (int i = 0; i < 5; ++i) { double lo, hi; in >> lo >> hi; r0[i] = interval(lo, hi); }

  try {
    IMap vf("var:p,q,k,m,th;fun:k*q,-1-k*p,-m,-0.25*p,-k;");
    IOdeSolver solver(vf, order);
    ITimeMap tm(solver);
    IVector xm(5), rr(5);
    for (int i = 0; i < 5; ++i) { xm[i] = interval(x[i].mid().leftBound()); rr[i] = x[i] - xm[i]; }
    C0Rect2Set set(xm, C, r0, rr);   // point centre + C*r0 + r (r encloses x - xm)

    // Stage 1: fixed-time map to tau1, enclosure of the trajectory segment is checked for q
    IVector y1 = tm(tau1, set);
    cout << "STAGE1 tau1 = " << tau1 << "\n";
    for (int i = 0; i < 5; ++i) cout << "  y1[" << i << "] = " << y1[i] << "\n";
    bool qneg = (y1[1].rightBound() < 0.0);
    cout << "  q(tau1) < 0 : " << (qneg ? "YES" : "NO") << "\n";
    if (!qneg) { cout << "RESULT FAIL q not negative at tau1\n"; return 1; }

    // Stage 2: first return to the section q = 0 from q < 0 (MinusPlus)
    ICoordinateSection section(5, 1);           // q = y[1] = 0
    IPoincareMap pm(solver, section, poincare::MinusPlus);
    interval rt;
    IVector y2 = pm(set, rt);
    cout << "STAGE2 return time (after tau1) = " << rt << "\n";
    const char* nm[5] = {"p", "q", "k", "m", "th"};
    for (int i = 0; i < 5; ++i) cout << "  " << nm[i] << " = " << y2[i] << "\n";
    cout << "  m sign: " << (y2[3].leftBound() > 0 ? "POSITIVE" : (y2[3].rightBound() < 0 ? "NEGATIVE" : "UNDETERMINED")) << "\n";
    cout << "RESULT OK\n";
  } catch (std::exception& e) {
    cout << "RESULT EXCEPTION: " << e.what() << "\n";
    return 1;
  }
  return 0;
}
