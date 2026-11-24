// Rigorous Evans-function enclosure on a lambda-segment  lam(t) = lam_a + t (lam_b - lam_a),  t in [0,1].
// State (17): (p,q,k,m, W01r,W01i,W02r,W02i,W03r,W03i,W12r,W12i,W13r,W13i,W23r,W23i, t), reversed time tau = 8 - sigma.
// Profile (reversed): p' = k q, q' = -1 - k p, k' = -m, m' = -p/4.   Second compound of the eigen-ODE
//   u^(4) = a2 u'' + a1 u' + a0 u,  a0 = 1/4 - 3m^2 - kp - lam,  a1 = -(5km + q/4),  a2 = -k^2  (minus sign: reversed time).
// Stage B (sigma0 -> 8): W rescaled by exp(-int (q/4)^{1/3}) (a positive factor), profile box from the tail theorem.
// Stage C (8 -> tip): tight step-1 profile set; the W and t rows of the stage-B Lohner set are carried over (so the
//   correlation with t is kept), W multiplied by a power of two.  Then time tau1 and the Poincare map to q = 0 (MinusPlus).
// Input:  lra lia dr di / tauB / 17 lines xB / 17x17 CB / 17 lines r0B / tau1 / 4 lines xC / 4x4 CC / 4 lines r0C
#include <iostream>
#include <fstream>
#include <cmath>
#include "capd/capdlib.h"
using namespace std;
using namespace capd;
#ifndef SETB
#define SETB C0Rect2Set
#endif
static const char* VF = "par:lra,lia,dr,di;var:p,q,k,m,a1r,a1i,a2r,a2i,a3r,a3i,b2r,b2i,b3r,b3i,c3r,c3i,t;fun:k*q,-1-k*p,-m,-0.25*p,-a2r,-a2i,-(a3r+b2r),-(a3i+b2i),-(-(5*k*m+0.25*q)*a1r - k*k*a2r + b3r),-(-(5*k*m+0.25*q)*a1i - k*k*a2i + b3i),-b3r,-b3i,-(-(0.25-3*m*m-k*p)*a1r + ((lra+t*dr)*a1r - (lia+t*di)*a1i) - k*k*b2r + c3r),-(-(0.25-3*m*m-k*p)*a1i + ((lra+t*dr)*a1i + (lia+t*di)*a1r) - k*k*b2i + c3i),-(-(0.25-3*m*m-k*p)*a2r + ((lra+t*dr)*a2r - (lia+t*di)*a2i) + (5*k*m+0.25*q)*b2r),-(-(0.25-3*m*m-k*p)*a2i + ((lra+t*dr)*a2i + (lia+t*di)*a2r) + (5*k*m+0.25*q)*b2i),0;";
static const char* VFB = "par:lra,lia,dr,di;var:p,q,k,m,a1r,a1i,a2r,a2i,a3r,a3i,b2r,b2i,b3r,b3i,c3r,c3i,t;fun:k*q,-1-k*p,-m,-0.25*p,(-a2r) - (0.25*q)^(1/3)*a1r,(-a2i) - (0.25*q)^(1/3)*a1i,(-(a3r+b2r)) - (0.25*q)^(1/3)*a2r,(-(a3i+b2i)) - (0.25*q)^(1/3)*a2i,(-(-(5*k*m+0.25*q)*a1r - k*k*a2r + b3r)) - (0.25*q)^(1/3)*a3r,(-(-(5*k*m+0.25*q)*a1i - k*k*a2i + b3i)) - (0.25*q)^(1/3)*a3i,(-b3r) - (0.25*q)^(1/3)*b2r,(-b3i) - (0.25*q)^(1/3)*b2i,(-(-(0.25-3*m*m-k*p)*a1r + ((lra+t*dr)*a1r - (lia+t*di)*a1i) - k*k*b2r + c3r)) - (0.25*q)^(1/3)*b3r,(-(-(0.25-3*m*m-k*p)*a1i + ((lra+t*dr)*a1i + (lia+t*di)*a1r) - k*k*b2i + c3i)) - (0.25*q)^(1/3)*b3i,(-(-(0.25-3*m*m-k*p)*a2r + ((lra+t*dr)*a2r - (lia+t*di)*a2i) + (5*k*m+0.25*q)*b2r)) - (0.25*q)^(1/3)*c3r,(-(-(0.25-3*m*m-k*p)*a2i + ((lra+t*dr)*a2i + (lia+t*di)*a2r) + (5*k*m+0.25*q)*b2i)) - (0.25*q)^(1/3)*c3i,0;";

int main(int argc, char* argv[]) {
  cout.precision(17);
  if (argc < 2) { cerr << "usage: evans_seg input [order]\n"; return 2; }
  int order = (argc > 2) ? atoi(argv[2]) : 25;
  const int n = 17;
  ifstream in(argv[1]);
  double lra, lia, dr, di; in >> lra >> lia >> dr >> di;
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
  const char* nm[17] = {"p","q","k","m","W01r","W01i","W02r","W02i","W03r","W03i","W12r","W12i","W13r","W13i","W23r","W23i","t"};
  try {
    IMap vf(VF), vfb(VFB);
    IMap* fs[2] = {&vf, &vfb};
    for (int u = 0; u < 2; ++u) { fs[u]->setParameter("lra", interval(lra)); fs[u]->setParameter("lia", interval(lia));
                                  fs[u]->setParameter("dr", interval(dr)); fs[u]->setParameter("di", interval(di)); }
    IOdeSolver solverB(vfb, order), solver(vf, order);
    ITimeMap tmB(solverB), tm(solver);
    // ---------------- stage B
    IVector xm(n), rx(n);
    for (int i = 0; i < n; ++i) { xm[i] = interval(xB[i].mid().leftBound()); rx[i] = xB[i] - xm[i]; }
    SETB setB(xm, CB, r0B, rx);
    tmB(interval(tauB), setB);
    IVector x = setB.get_x(), r0 = setB.get_r0(), r = setB.get_r();
    IMatrix C = setB.get_C(), B = setB.get_B();
    IVector Br = B * r;
    IVector hullB = IVector(setB);
    cout << "STAGEB q(8) = " << hullB[1] << "  t = " << hullB[16] << "\n";
    if (getenv("EVDEBUG")) { cout << "r0 = " << r0 << "\nr = " << r << "\nx = " << x << "\n"; for (int i=0;i<n;++i){ cout << "C" << i << ": "; for(int j=0;j<n;++j) cout << C[i][j].mid().leftBound() << " "; cout << "\n";} }
    double mx = 0;
    for (int i = 4; i < 16; ++i) mx = std::max(mx, std::max(std::fabs(hullB[i].leftBound()), std::fabs(hullB[i].rightBound())));
    if (!(mx > 0) || !std::isfinite(mx)) { cout << "RESULT FAIL stageB\n"; return 1; }
    int e2 = (int)std::floor(std::log2(1e-8 / mx));
    interval sc(std::ldexp(1.0, e2));
    // ---------------- stage C set
    IVector xn(n), r0n(n), rn(n); IMatrix Cn(n, n);
    for (int i = 0; i < 4; ++i) {
      xn[i] = interval(xC[i].mid().leftBound()); rn[i] = xC[i] - xn[i];
      for (int j = 0; j < 4; ++j) Cn[i][j] = CC[i][j];
      r0n[i] = r0C[i];
    }
    int col = 4;
    int nused = 0; for (int j = 0; j < n; ++j) if (!(r0[j].leftBound() == 0.0 && r0[j].rightBound() == 0.0)) ++nused;
    if (nused + 4 > n) { cout << "RESULT FAIL columns " << nused << "\n"; return 1; }
    for (int j = 0; j < n; ++j) {
      bool used = !(r0[j].leftBound() == 0.0 && r0[j].rightBound() == 0.0);
      if (!used) continue;
      if (col >= n) { cout << "RESULT FAIL too many columns\n"; return 1; }
      for (int i = 4; i < n; ++i) Cn[i][col] = ((i < 16) ? sc : interval(1.0)) * C[i][j];
      r0n[col] = r0[j];
      ++col;
    }
    for (int i = 4; i < n; ++i) {
      interval s = (i < 16) ? sc : interval(1.0);
      interval xi = s * x[i];
      xn[i] = interval(xi.mid().leftBound());
      rn[i] = (xi - xn[i]) + s * Br[i];
    }
    C0Rect2Set setC(xn, Cn, r0n, rn);
    IVector y1 = tm(interval(tau1), setC);
    bool qneg = (y1[1].rightBound() < 0.0);
    cout << "STAGEC q(tau1) = " << y1[1] << (qneg ? "  (<0 ok)" : "  NOT NEGATIVE") << "\n";
    if (!qneg) { cout << "RESULT FAIL\n"; return 1; }
    ICoordinateSection section(n, 1);
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
