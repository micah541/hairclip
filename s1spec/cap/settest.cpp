// compare set types for stage B (sigma0 -> 8)
#include <iostream>
#include <fstream>
#include <cmath>
#include "capd/capdlib.h"
using namespace std; using namespace capd;
extern const char* VF;
int main(int argc, char* argv[]) {
  cout.precision(10);
  const int n = 16; int order = atoi(argv[2]); double h = atof(argv[3]); int kind = atoi(argv[4]);
  ifstream in(argv[1]); double a,b,c,d; in >> a >> b >> c >> d; double tauB; in >> tauB;
  IVector xB(n), r0B(n); IMatrix CB(n,n);
  for (int i=0;i<n;++i){double lo,hi; in>>lo>>hi; xB[i]=interval(lo,hi);}
  for (int i=0;i<n;++i) for(int j=0;j<n;++j){double v; in>>v; CB[i][j]=interval(v);}
  for (int i=0;i<n;++i){double lo,hi; in>>lo>>hi; r0B[i]=interval(lo,hi);}
  IMap vf(VF); vf.setParameter("lr", interval(a,b)); vf.setParameter("li", interval(c,d));
  IOdeSolver solver(vf, order); if (h>0){ solver.setStep(h); solver.turnOffStepControl(); }
  ITimeMap tm(solver);
  IVector xm(n), rx(n); for(int i=0;i<n;++i){ xm[i]=interval(xB[i].mid().leftBound()); rx[i]=xB[i]-xm[i]; }
  IVector y;
  if (kind==0){ C0Rect2Set s(xm,CB,r0B,rx); y = tm(interval(tauB), s); }
  if (kind==1){ C0TripletonSet s(xm,CB,r0B,rx); y = tm(interval(tauB), s); }
  if (kind==2){ C0HORect2Set s(xm,CB,r0B,rx); y = tm(interval(tauB), s); }
  if (kind==3){ C0HOTripletonSet s(xm,CB,r0B,rx); y = tm(interval(tauB), s); }
  if (kind==4){ C1Rect2Set s(xm,CB,r0B,rx); y = tm(interval(tauB), s); }
  double w = 0, m = 0; for (int i=4;i<n;++i){ w = std::max(w, (double)diam(y[i]).rightBound()); m = std::max(m, std::max(std::fabs(y[i].leftBound()), std::fabs(y[i].rightBound()))); }
  cout << "kind " << kind << " order " << order << " h " << h << ": max diam/max|W| = " << w/m << "\n";
}
