// Rigorous enclosures of the S1 half-branch trajectory (reversed time tau = S - sigma).
// Same input format as s1proof.cpp.  Usage: s1curve input Tend Nsub [order]
// Prints:
//   TIP <return time lo hi> then 5 lines "lo hi" for (p,q,k,m,th) at the q=0 crossing,
//   then Nsub lines: "tlo thi  p_lo p_hi  q_lo q_hi  k_lo k_hi  m_lo m_hi  th_lo th_hi"
//   enclosing the trajectory for tau in [tlo, thi].
#include <iostream>
#include <fstream>
#include <iomanip>
#include "capd/capdlib.h"
using namespace std;
using namespace capd;

int main(int argc, char* argv[]) {
  cout.precision(17);
  if (argc < 4) { cerr << "usage: s1curve input Tend Nsub [order]\n"; return 2; }
  double Tend = atof(argv[2]); int Nsub = atoi(argv[3]);
  int order = (argc > 4) ? atoi(argv[4]) : 25;
  ifstream in(argv[1]);
  double tau1; in >> tau1;
  IVector x(5), r0(5); IMatrix C(5, 5);
  for (int i = 0; i < 5; ++i) { double lo, hi; in >> lo >> hi; x[i] = interval(lo, hi); }
  for (int i = 0; i < 5; ++i) for (int j = 0; j < 5; ++j) { double c; in >> c; C[i][j] = interval(c); }
  for (int i = 0; i < 5; ++i) { double lo, hi; in >> lo >> hi; r0[i] = interval(lo, hi); }
  try {
    IMap vf("var:p,q,k,m,th;fun:k*q,-1-k*p,-m,-0.25*p,-k;");
    {   // tip via Poincare map (same two-stage procedure as s1proof)
      IOdeSolver solver(vf, order);
      ITimeMap tm(solver);
      IVector xm(5), rr(5);
    for (int i = 0; i < 5; ++i) { xm[i] = interval(x[i].mid().leftBound()); rr[i] = x[i] - xm[i]; }
    C0Rect2Set set(xm, C, r0, rr);   // point centre + C*r0 + r (r encloses x - xm)
      IVector y1 = tm(tau1, set);
      if (!(y1[1].rightBound() < 0.0)) { cout << "FAIL q(tau1) not negative\n"; return 1; }
      ICoordinateSection section(5, 1);
      IPoincareMap pm(solver, section, poincare::MinusPlus);
      interval rt; IVector y2 = pm(set, rt);
      cout << "TIP " << rt.leftBound() << " " << rt.rightBound() << "\n";
      for (int i = 0; i < 5; ++i) cout << y2[i].leftBound() << " " << y2[i].rightBound() << "\n";
    }
    IOdeSolver solver(vf, order);
    ITimeMap tm(solver);
    IVector xm(5), rr(5);
    for (int i = 0; i < 5; ++i) { xm[i] = interval(x[i].mid().leftBound()); rr[i] = x[i] - xm[i]; }
    C0Rect2Set set(xm, C, r0, rr);   // point centre + C*r0 + r (r encloses x - xm)
    ITimeMap::SolutionCurve sol(0.0);
    tm(Tend, set, sol);
    for (int i = 0; i < Nsub; ++i) {
      interval t(Tend * i / Nsub, Tend * (i + 1) / Nsub);
      if (i == Nsub - 1) t = interval(Tend * i / Nsub, Tend);
      IVector v = sol(t);
      cout << t.leftBound() << " " << t.rightBound();
      for (int j = 0; j < 5; ++j) cout << " " << v[j].leftBound() << " " << v[j].rightBound();
      cout << "\n";
    }
    cout << "END\n";
  } catch (std::exception& e) {
    cout << "EXCEPTION: " << e.what() << "\n"; return 1;
  }
  return 0;
}
