// experiment: stage B with exact cone coefficients (k=m=p=0), state (W 12, q, t)
#include <iostream>
#include <fstream>
#include <cmath>
#include "capd/capdlib.h"
using namespace std; using namespace capd;
int main(int argc, char* argv[]) {
  const char* W[12] = {"a1r","a1i","a2r","a2i","a3r","a3i","b2r","b2i","b3r","b3i","c3r","c3i"};
  int damp = atoi(argv[2]);
  string LR = "(lra+t*dr)", LI = "(lia+t*di)";
  string f[14] = {"-a2r", "-a2i", "-(a3r+b2r)", "-(a3i+b2i)", "-(-(0.25*q)*a1r + b3r)", "-(-(0.25*q)*a1i + b3i)", "-b3r", "-b3i",
    "-(-0.25*a1r + ("+LR+"*a1r - "+LI+"*a1i) + c3r)", "-(-0.25*a1i + ("+LR+"*a1i + "+LI+"*a1r) + c3i)",
    "-(-0.25*a2r + ("+LR+"*a2r - "+LI+"*a2i) + (0.25*q)*b2r)", "-(-0.25*a2i + ("+LR+"*a2i + "+LI+"*a2r) + (0.25*q)*b2i)", "-1", "0"};
  string s = "par:lra,lia,dr,di;var:"; for (int i=0;i<12;++i) s += string(W[i]) + ","; s += "q,t;fun:";
  for (int i=0;i<14;++i) { string g = f[i]; if (damp && i < 12) g = "(" + g + ") - " + string(damp == 2 ? "1.6" : "(0.25*q)^(1/3)") + "*" + W[i]; s += g + (i<13 ? "," : ";"); }
  // read the evans_seg input and map it
  ifstream in(argv[1]);
  double l0,l1,i0,i1,dr,di; in >> l0 >> l1 >> i0 >> i1 >> dr >> di; double tauB; in >> tauB;
  const int n0 = 17; IVector xB(n0), r0B(n0); IMatrix CB(n0,n0);
  for (int i=0;i<n0;++i){double lo,hi; in>>lo>>hi; xB[i]=interval(lo,hi);}
  for (int i=0;i<n0;++i) for(int j=0;j<n0;++j){double v; in>>v; CB[i][j]=interval(v);}
  for (int i=0;i<n0;++i){double lo,hi; in>>lo>>hi; r0B[i]=interval(lo,hi);}
  const int n = 14; IVector x(n), r0(n); IMatrix C(n,n);
  // W rows 4..15 -> 0..11, q row 1 -> 12, t row 16 -> 13 ; columns 4..12 (Gamma 8 + t) -> 0..8
  for (int i=0;i<12;++i) x[i] = xB[4+i]; x[12] = interval(32.0); x[13] = xB[16];
  for (int i=0;i<12;++i) for (int j=0;j<9;++j) C[i][j] = CB[4+i][4+j];
  for (int j=0;j<9;++j) { C[13][j] = CB[16][4+j]; r0[j] = r0B[4+j]; }
  IMap vf(s); vf.setParameter("lra", interval(l0,l1)); vf.setParameter("lia", interval(i0,i1)); vf.setParameter("dr", interval(dr)); vf.setParameter("di", interval(di));
  IOdeSolver solver(vf, 25); if (argc > 3) { solver.setStep(atof(argv[3])); solver.turnOffStepControl(); } ITimeMap tm(solver);
  IVector xm(n), rx(n); for (int i=0;i<n;++i){ xm[i]=interval(x[i].mid().leftBound()); rx[i]=x[i]-xm[i]; }
  C0Rect2Set set(xm, C, r0, rx);
  IVector y = tm(interval(tauB), set);
  double mw=0, md=0; for (int i=0;i<12;++i){ mw = std::max(mw, std::max(std::fabs(y[i].leftBound()), std::fabs(y[i].rightBound()))); md = std::max(md, (double)diam(y[i]).rightBound()); }
  cout << "damp " << damp << " q=" << y[12] << " rel width " << md/mw << "\n";
}
