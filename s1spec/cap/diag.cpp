#include <iostream>
#include <fstream>
#include <cmath>
#include "capd/capdlib.h"
using namespace std; using namespace capd;
extern const char* VF_D;
int main(int argc, char* argv[]) {
  cout.precision(6);
  const int n = 16; int order = 20; double h = atof(argv[2]);
  ifstream in(argv[1]); double a,b,c,d; in >> a >> b >> c >> d; double tauB; in >> tauB;
  IVector xB(n), r0B(n); IMatrix CB(n,n);
  for (int i=0;i<n;++i){double lo,hi; in>>lo>>hi; xB[i]=interval(lo,hi);}
  for (int i=0;i<n;++i) for(int j=0;j<n;++j){double v; in>>v; CB[i][j]=interval(v);}
  for (int i=0;i<n;++i){double lo,hi; in>>lo>>hi; r0B[i]=interval(lo,hi);}
  IMap vf(VF_D); vf.setParameter("lr", interval(a,b)); vf.setParameter("li", interval(c,d));
  IOdeSolver solver(vf, order); if (h>0){ solver.setStep(h); solver.turnOffStepControl(); }
  IVector xm(n), rx(n); for(int i=0;i<n;++i){ xm[i]=interval(xB[i].mid().leftBound()); rx[i]=xB[i]-xm[i]; }
  C0Rect2Set s(xm,CB,r0B,rx);
  double T = 0; int steps = 0;
  while (T < tauB - 1e-12) {
    s.move(solver); steps++; T = s.getCurrentTime().rightBound();
    if (steps % 40 == 0 || T >= tauB - 1e-12) {
      IVector x = s.get_x(); IMatrix C = s.get_C(); IVector r0 = s.get_r0(); IMatrix B = s.get_B(); IVector r = s.get_r();
      IVector Cr = C*r0, Br = B*r; IVector hull = IVector(s);
      double m = 0, wc = 0, wb = 0;
      for (int i=4;i<n;++i){ m = std::max(m, std::fabs(x[i].leftBound())); wc = std::max(wc, (double)rightBound(abs(Cr[i]))); wb = std::max(wb, (double)rightBound(abs(Br[i]))); }
      cout << "t=" << T << " |W|~" << m << "  |C r0|/|W|=" << wc/m << "  |B r|/|W|=" << wb/m << "\n";
    }
  }
}
