#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include "capd/capdlib.h"
using namespace std; using namespace capd;
int main(int argc, char* argv[]) {
  cout.precision(10);
  ifstream in(argv[1]); int order = atoi(argv[2]);
  int n, np; in >> n >> np; vector<string> pn(np); vector<interval> pv(np);
  for (int i = 0; i < np; ++i) { double lo, hi; in >> pn[i] >> lo >> hi; pv[i] = interval(lo, hi); }
  string line; getline(in, line); string VFB, VFC; getline(in, VFB); getline(in, VFC);
  double tauB; in >> tauB; int nch, npc; in >> nch >> npc; vector<string> pcn(npc); for (auto& s : pcn) in >> s;
  vector<double> L(nch), Hh(nch); vector<vector<interval>> P(nch, vector<interval>(npc));
  for (int c = 0; c < nch; ++c) { in >> L[c] >> Hh[c]; for (int i = 0; i < npc; ++i) { double a, b; in >> a >> b; P[c][i] = interval(a, b); } }
  IVector x(n), box(n), r0(n); IMatrix C(n, n);
  for (int i = 0; i < n; ++i) { double a, b; in >> a >> b; x[i] = interval(a, b); }
  int ncol; in >> ncol;
  for (int i = 0; i < n; ++i) for (int j = 0; j < ncol; ++j) { double v; in >> v; C[i][j] = interval(v); }
  for (int j = 0; j < ncol; ++j) { double a, b; in >> a >> b; r0[j] = interval(a, b); }
  for (int i = 0; i < n; ++i) { double a, b; in >> a >> b; box[i] = interval(a, b); }
  IMap vf(VFB.c_str()); for (int i = 0; i < np; ++i) vf.setParameter(pn[i].c_str(), pv[i]);
  IVector xc(n); for (int i = 0; i < n; ++i) xc[i] = interval(x[i].mid().leftBound());
  C0Rect2Set S(xc, C, r0, box + (x - xc));
  double done = 0;
  for (int c = 0; c < nch; ++c) {
    for (int i = 0; i < npc; ++i) vf.setParameter(pcn[i].c_str(), P[c][i]);
    IOdeSolver sb(vf, order); sb.setStep(Hh[c]); sb.turnOffStepControl();
    ITimeMap tm(sb); tm(interval(done + L[c]), S); done += L[c];
    IVector xx = S.get_x(), cr = S.get_C() * S.get_r0(), br = S.get_B() * S.get_r();
    double ax = 0, ac = 0, ab = 0;
    for (int i = 2; i < n; ++i) { ax = std::max(ax, (double)rightBound(abs(xx[i]))); ac = std::max(ac, (double)rightBound(abs(cr[i]))); ab = std::max(ab, (double)rightBound(abs(br[i]))); }
    if (c % 3 == 2) cout << "sigma=" << 32 - done << " |x|=" << ax << " |Cr0|=" << ac << " |Br|=" << ab << " brX=" << br[0] << " brE=" << br[1] << "\n";
  }
}
