// Tight S1 profile data for downstream certificates.  Same input format as s1proof.cpp.
// Usage: s1prof input Tend N [order]
// Output: N lines  "ta tb  (p q k m th at tc: lo hi)x5  (p q k m th hull over [ta,tb]: lo hi)x5"
// with tc = 0.5*(ta+tb) computed in double and printed; reversed time tau = 8 - sigma.
#include <iostream>
#include <fstream>
#include "capd/capdlib.h"
using namespace std;
using namespace capd;
int main(int argc, char* argv[]) {
  cout.precision(17);
  if (argc < 4) { cerr << "usage: s1prof input Tend N [order]\n"; return 2; }
  double Tend = atof(argv[2]); int N = atoi(argv[3]);
  int order = (argc > 4) ? atoi(argv[4]) : 25;
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
    C0Rect2Set set(xm, C, r0, rr);
    ITimeMap::SolutionCurve sol(0.0);
    tm(Tend, set, sol);
    for (int i = 0; i < N; ++i) {
      double ta = Tend * i / N, tb = (i == N - 1) ? Tend : Tend * (i + 1) / N;
      double tc = 0.5 * (ta + tb);
      IVector vc = sol(interval(tc));
      IVector vh = sol(interval(ta, tb));
      cout << ta << " " << tb << " " << tc;
      for (int j = 0; j < 5; ++j) cout << " " << vc[j].leftBound() << " " << vc[j].rightBound();
      for (int j = 0; j < 5; ++j) cout << " " << vh[j].leftBound() << " " << vh[j].rightBound();
      cout << "\n";
    }
    cout << "END\n";
  } catch (std::exception& e) { cout << "EXCEPTION: " << e.what() << "\n"; return 1; }
  return 0;
}
