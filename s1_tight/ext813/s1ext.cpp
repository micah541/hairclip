// Outward extension of the tight S1 profile (prof/s1prof.cpp) to sigma in [8, 8 + Send].
// Same input file as s1prof (in_br2.txt: the certified parallelepiped x + C r0 at sigma = 8, all eps in the
// tight bracket) and the same state (p, q, k, m, th).  Integrates the profile ODE forward in s = sigma - 8:
//     p' = -k q,  q' = 1 + k p,  k' = m,  m' = p/4,  th' = k        (' = d/dsigma = -d/dtau)
// i.e. exactly the s1prof vector field with the time reversed (tau = 8 - sigma = -s).
// Usage: s1ext input Send N [order [pointsfile NP]]
// stdout: N rows in the prof_tight format, in tau = 8 - sigma = -s:
//     ta tb tc  (p q k m th at tc: lo hi)x5  (p q k m th hull over [ta,tb]: lo hi)x5
//   s_a = Send*i/N, s_b = Send*(i+1)/N (s_b = Send on the last row), s_c = 0.5*(s_a+s_b), all in double;
//   printed ta = -s_b, tb = -s_a, tc = -s_c (exact negations), 17 significant digits; then "END".
// pointsfile: enclosures at the exact points s_j = Send*j/NP, j = 0..NP:  "s_j (p q k m th: lo hi)x5".
#include <iostream>
#include <fstream>
#include "capd/capdlib.h"
using namespace std;
using namespace capd;
int main(int argc, char* argv[]) {
  cout.precision(17);
  if (argc < 4) { cerr << "usage: s1ext input Send N [order [pointsfile NP]]\n"; return 2; }
  double Send = atof(argv[2]); int N = atoi(argv[3]);
  int order = (argc > 4) ? atoi(argv[4]) : 25;
  ifstream in(argv[1]);
  double tau1; in >> tau1;
  IVector x(5), r0(5); IMatrix C(5, 5);
  for (int i = 0; i < 5; ++i) { double lo, hi; in >> lo >> hi; x[i] = interval(lo, hi); }
  for (int i = 0; i < 5; ++i) for (int j = 0; j < 5; ++j) { double c; in >> c; C[i][j] = interval(c); }
  for (int i = 0; i < 5; ++i) { double lo, hi; in >> lo >> hi; r0[i] = interval(lo, hi); }
  if (!in) { cerr << "bad input file\n"; return 2; }
  try {
    IMap vf("var:p,q,k,m,th;fun:-k*q,1+k*p,m,0.25*p,k;");
    IOdeSolver solver(vf, order);
    ITimeMap tm(solver);
    IVector xm(5), rr(5);
    for (int i = 0; i < 5; ++i) { xm[i] = interval(x[i].mid().leftBound()); rr[i] = x[i] - xm[i]; }
    C0Rect2Set set(xm, C, r0, rr);
    ITimeMap::SolutionCurve sol(0.0);
    tm(Send, set, sol);
    for (int i = 0; i < N; ++i) {
      double sa = Send * i / N, sb = (i == N - 1) ? Send : Send * (i + 1) / N;
      double sc = 0.5 * (sa + sb);
      IVector vc = sol(interval(sc));
      IVector vh = sol(interval(sa, sb));
      cout << -sb << " " << (i == 0 ? 0.0 : -sa) << " " << -sc;
      for (int j = 0; j < 5; ++j) cout << " " << vc[j].leftBound() << " " << vc[j].rightBound();
      for (int j = 0; j < 5; ++j) cout << " " << vh[j].leftBound() << " " << vh[j].rightBound();
      cout << "\n";
    }
    cout << "END\n";
    if (argc > 6) {
      ofstream pf(argv[5]); pf.precision(17);
      int NP = atoi(argv[6]);
      for (int j = 0; j <= NP; ++j) {
        double s = (j == NP) ? Send : Send * j / NP;
        IVector v = sol(interval(s));
        pf << s;
        for (int c = 0; c < 5; ++c) pf << " " << v[c].leftBound() << " " << v[c].rightBound();
        pf << "\n";
      }
    }
  } catch (std::exception& e) { cout << "EXCEPTION: " << e.what() << "\n"; return 1; }
  return 0;
}
