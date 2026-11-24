// Rigorous Evans function enclosure (step 2), two stages.
//  Stage B: sigma0 -> 8 (fixed reversed time tauB = sigma0 - 8): profile box from the step-1 tail theorem at sigma0,
//           W from the Lyapunov-Perron admissible plane at sigma0 (lintail2.py).
//  Stage C: restart at sigma = 8 with the tight step-1 profile set (shooting bracket) and W(8) from stage B (rescaled
//           by a positive factor), integrate to tau1, then Poincare map to the tip section q = 0 (MinusPlus).
// State: (p,q,k,m, W01r,W01i,W02r,W02i,W03r,W03i,W12r,W12i,W13r,W13i,W23r,W23i); reversed-time compound field.
// Input:  lr_lo lr_hi li_lo li_hi
//         tauB
//         16 lines xB_lo xB_hi ; 256 numbers CB ; 16 lines r0B_lo r0B_hi        (stage-B initial set xB + CB r0B)
//         tau1
//         4 lines xC_lo xC_hi ; 16 numbers CC ; 4 lines r0C_lo r0C_hi           (stage-C profile set)
// Output: enclosure of W at sigma = 8 (after rescaling) and at the tip.
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include "capd/capdlib.h"
using namespace std;
using namespace capd;

static const char* VF =
  "par:lr,li;"
  "var:p,q,k,m,a1r,a1i,a2r,a2i,a3r,a3i,b2r,b2i,b3r,b3i,c3r,c3i;"
  "fun:k*q, -1-k*p, -m, -0.25*p,"
  "-a2r, -a2i,"
  "-(a3r+b2r), -(a3i+b2i),"
  "-(-(5*k*m+0.25*q)*a1r - k*k*a2r + b3r), -(-(5*k*m+0.25*q)*a1i - k*k*a2i + b3i),"
  "-b3r, -b3i,"
  "-(-(0.25-3*m*m-k*p)*a1r + (lr*a1r - li*a1i) - k*k*b2r + c3r),"
  "-(-(0.25-3*m*m-k*p)*a1i + (lr*a1i + li*a1r) - k*k*b2i + c3i),"
  "-(-(0.25-3*m*m-k*p)*a2r + (lr*a2r - li*a2i) + (5*k*m+0.25*q)*b2r),"
  "-(-(0.25-3*m*m-k*p)*a2i + (lr*a2i + li*a2r) + (5*k*m+0.25*q)*b2i);";

// Stage B uses the positively rescaled field W~ = exp(-int (q/4)^{1/3}) W (same phase; keeps |W~| ~ 1)
static const char* VFB = "par:lr,li;var:p,q,k,m,a1r,a1i,a2r,a2i,a3r,a3i,b2r,b2i,b3r,b3i,c3r,c3i;fun:k*q,-1-k*p,-m,-0.25*p,(-a2r) - (0.25*q)^(1/3)*a1r,(-a2i) - (0.25*q)^(1/3)*a1i,(-(a3r+b2r)) - (0.25*q)^(1/3)*a2r,(-(a3i+b2i)) - (0.25*q)^(1/3)*a2i,(-(-(5*k*m+0.25*q)*a1r - k*k*a2r + b3r)) - (0.25*q)^(1/3)*a3r,(-(-(5*k*m+0.25*q)*a1i - k*k*a2i + b3i)) - (0.25*q)^(1/3)*a3i,(-b3r) - (0.25*q)^(1/3)*b2r,(-b3i) - (0.25*q)^(1/3)*b2i,(-(-(0.25-3*m*m-k*p)*a1r + (lr*a1r - li*a1i) - k*k*b2r + c3r)) - (0.25*q)^(1/3)*b3r,(-(-(0.25-3*m*m-k*p)*a1i + (lr*a1i + li*a1r) - k*k*b2i + c3i)) - (0.25*q)^(1/3)*b3i,(-(-(0.25-3*m*m-k*p)*a2r + (lr*a2r - li*a2i) + (5*k*m+0.25*q)*b2r)) - (0.25*q)^(1/3)*c3r,(-(-(0.25-3*m*m-k*p)*a2i + (lr*a2i + li*a2r) + (5*k*m+0.25*q)*b2i)) - (0.25*q)^(1/3)*c3i;";

int main(int argc, char* argv[]) {
  cout.precision(17);
  if (argc < 2) { cerr << "usage: evans_cap2 input [order]\n"; return 2; }
  int order = (argc > 2) ? atoi(argv[2]) : 25;
  const int n = 16;
  ifstream in(argv[1]);
  double a, b, c, d; in >> a >> b >> c >> d;
  interval lr(a, b), li(c, d);
  double tauB; in >> tauB;
  IVector xB(n), r0B(n); IMatrix CB(n, n);
  for (int i = 0; i < n; ++i) { double lo, hi; in >> lo >> hi; xB[i] = interval(lo, hi); }
  for (int i = 0; i < n; ++i) for (int j = 0; j < n; ++j) { double v; in >> v; CB[i][j] = interval(v); }
  for (int i = 0; i < n; ++i) { double lo, hi; in >> lo >> hi; r0B[i] = interval(lo, hi); }
  double tau1; in >> tau1;
  IVector xC(4), r0C(4); IMatrix CC(4, 4);
  for (int i = 0; i < 4; ++i) { double lo, hi; in >> lo >> hi; xC[i] = interval(lo, hi); }
  for (int i = 0; i < 4; ++i) for (int j = 0; j < 4; ++j) { double v; in >> v; CC[i][j] = interval(v); }
  for (int i = 0; i < 4; ++i) { double lo, hi; in >> lo >> hi; r0C[i] = interval(lo, hi); }
  if (!in) { cerr << "bad input\n"; return 2; }
  const char* nm[16] = {"p","q","k","m","W01r","W01i","W02r","W02i","W03r","W03i","W12r","W12i","W13r","W13i","W23r","W23i"};
  try {
    IMap vf(VF), vfb(VFB);
    vf.setParameter("lr", lr); vf.setParameter("li", li);
    vfb.setParameter("lr", lr); vfb.setParameter("li", li);
    IOdeSolver solverB(vfb, order);
    ITimeMap tmB(solverB);
    IOdeSolver solver(vf, order);
    double hfix = (argc > 3) ? atof(argv[3]) : 0.0;
    if (hfix > 0) { solver.setStep(hfix); solver.turnOffStepControl(); }
    ITimeMap tm(solver);
    // ---- stage B
    IVector xm(n), rx(n);
    for (int i = 0; i < n; ++i) { xm[i] = interval(xB[i].mid().leftBound()); rx[i] = xB[i] - xm[i]; }
    C0Rect2Set setB(xm, CB, r0B, rx);
    IVector y8 = tmB(interval(tauB), setB);
    // rescale W by a positive power of two so that max|W| ~ 1e-8 (keeps CAPD enclosures tight)
    double mx = 0;
    for (int i = 4; i < n; ++i) mx = std::max(mx, std::max(std::fabs(y8[i].leftBound()), std::fabs(y8[i].rightBound())));
    if (!(mx > 0) || !std::isfinite(mx)) { cout << "RESULT FAIL stageB W not finite\n"; return 1; }
    int e2 = (int)floor(log2(1e-8 / mx));
    interval sc = interval(ldexp(1.0, e2));
    cout << "STAGEB q(8) = " << y8[1] << "  (rescale 2^" << e2 << ")\n";
    for (int i = 4; i < n; ++i) cout << "W8 " << nm[i] << " " << (sc * y8[i]).leftBound() << " " << (sc * y8[i]).rightBound() << "\n";
    // ---- stage C
    IVector xc(n), rc(n), r0(n); IMatrix C(n, n);
    for (int i = 0; i < 4; ++i) {
      xc[i] = interval(xC[i].mid().leftBound()); rc[i] = xC[i] - xc[i]; r0[i] = r0C[i];
      for (int j = 0; j < 4; ++j) C[i][j] = CC[i][j];
    }
    for (int i = 4; i < n; ++i) {
      interval w = sc * y8[i];
      xc[i] = interval(w.mid().leftBound()); rc[i] = w - xc[i]; r0[i] = interval(0.0); C[i][i] = interval(1.0);
    }
    C0Rect2Set setC(xc, C, r0, rc);
    IVector y1 = tm(interval(tau1), setC);
    bool qneg = (y1[1].rightBound() < 0.0);
    cout << "STAGEC q(tau1) = " << y1[1] << (qneg ? "  (<0 ok)" : "  NOT NEGATIVE") << "\n";
    if (!qneg) { cout << "RESULT FAIL\n"; return 1; }
    ICoordinateSection section(n, 1);
    if (hfix > 0) solver.turnOnStepControl();
    IPoincareMap pm(solver, section, poincare::MinusPlus);
    interval rt;
    IVector y2 = pm(setC, rt);
    cout << "RETURN " << rt << "\n";
    for (int i = 0; i < n; ++i) cout << "TIP " << nm[i] << " " << y2[i].leftBound() << " " << y2[i].rightBound() << "\n";
    cout << "RESULT OK\n";
  } catch (std::exception& e) {
    cout << "RESULT EXCEPTION: " << e.what() << "\n";
    return 1;
  }
  return 0;
}
