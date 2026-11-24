// Generic validated inward integration for the step-3 computations (CAPD, C0Rect2Set Lohner sets).
// Stage B: sigma 32 -> 8 (time tauB), chunks with fixed step h = min(HMAX, 3.2/sigma); after each chunk the
//          box part is re-boxed (parameter columns kept).
// Stage C: sigma 8 -> tip.  Profile rows 0..3 are replaced by the tight step-1 set; all other rows are boxed.
//          Integrate to TLO, then over the tip window [TLO, THI]; print the enclosure of every component over
//          the window.  Optionally monitor components along the whole path (enclosures over subintervals).
// Input (whitespace separated):
//   n  np  {name lo hi}^np  VFB-string-line  VFC-string-line  tauB
//   x: n x (lo hi)   ncol   C: n x ncol doubles   r0: ncol x (lo hi)   box: n x (lo hi)
//   profile set: 4 x (lo hi), 4x4 doubles, 4 x (lo hi)   TLO THI   nmon idx...   dsub
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>
#include "capd/capdlib.h"
using namespace std;
using namespace capd;
#ifndef SETTYPE
#define SETTYPE C0Rect2Set
#endif
typedef capd::dynset::C0DoubletonSet<IMatrix, capd::dynset::FullQRWithPivoting<> > C0Rect2NoReorg;
typedef SETTYPE SetT;

static void monitor(const char* tag, ITimeMap::SolutionCurve& sol, double t0, double t1, double dsub,
                    const vector<int>& mon, double sig_at_0, double dir) {
  int N = std::max(1, (int)std::ceil((t1 - t0) / dsub));
  for (int i = 0; i < N; ++i) {
    double a = t0 + (t1 - t0) * i / N, b = t0 + (t1 - t0) * (i + 1) / N;
    IVector v = sol(interval(a, b));
    cout << "MON " << tag << " " << sig_at_0 + dir * a << " " << sig_at_0 + dir * b;
    for (int j : mon) cout << " " << v[j].leftBound() << " " << v[j].rightBound();
    cout << "\n";
  }
}

int main(int argc, char* argv[]) {
  cout.precision(17);
  if (argc < 2) { cerr << "usage: hsrun input [order]\n"; return 2; }
  int order = (argc > 2) ? atoi(argv[2]) : 20;
  ifstream in(argv[1]);
  int n, np; in >> n >> np;
  vector<string> pn(np); vector<interval> pv(np);
  for (int i = 0; i < np; ++i) { double lo, hi; in >> pn[i] >> lo >> hi; pv[i] = interval(lo, hi); }
  string line; getline(in, line);
  string VFB, VFC; getline(in, VFB); getline(in, VFC);
  double tauB; in >> tauB;
  int nch, npc; in >> nch >> npc;
  vector<string> pcn(npc); for (int i = 0; i < npc; ++i) in >> pcn[i];
  vector<double> chlen(nch), chh(nch); vector<vector<interval>> chp(nch, vector<interval>(npc));
  for (int c = 0; c < nch; ++c) { in >> chlen[c] >> chh[c];
    for (int i = 0; i < npc; ++i) { double lo, hi; in >> lo >> hi; chp[c][i] = interval(lo, hi); } }
  IVector x(n), box(n); for (int i = 0; i < n; ++i) { double lo, hi; in >> lo >> hi; x[i] = interval(lo, hi); }
  int ncol; in >> ncol;
  IMatrix C(n, n); IVector r0(n);
  for (int i = 0; i < n; ++i) for (int j = 0; j < ncol; ++j) { double v; in >> v; C[i][j] = interval(v); }
  for (int j = 0; j < ncol; ++j) { double lo, hi; in >> lo >> hi; r0[j] = interval(lo, hi); }
  for (int i = 0; i < n; ++i) { double lo, hi; in >> lo >> hi; box[i] = interval(lo, hi); }
  IVector px(4), pr0(4); IMatrix pC(4, 4);
  for (int i = 0; i < 4; ++i) { double lo, hi; in >> lo >> hi; px[i] = interval(lo, hi); }
  for (int i = 0; i < 4; ++i) for (int j = 0; j < 4; ++j) { double v; in >> v; pC[i][j] = interval(v); }
  for (int i = 0; i < 4; ++i) { double lo, hi; in >> lo >> hi; pr0[i] = interval(lo, hi); }
  double TLO, THI; in >> TLO >> THI;
  int nmon; in >> nmon; vector<int> mon(nmon); for (int i = 0; i < nmon; ++i) in >> mon[i];
  double dsub; in >> dsub;
  if (!in) { cerr << "bad input\n"; return 2; }
  const double HMAX = getenv("HMAX") ? atof(getenv("HMAX")) : 0.1;
  try {
    IMap vfb(VFB.c_str()), vfc(VFC.c_str());
    for (int i = 0; i < np; ++i) { vfb.setParameter(pn[i].c_str(), pv[i]); vfc.setParameter(pn[i].c_str(), pv[i]); }
    IOdeSolver solB(vfb, order), solC(vfc, order);
    // ---------------- stage B
    IVector xc(n), rr(n);
    for (int i = 0; i < n; ++i) { xc[i] = interval(x[i].mid().leftBound()); rr[i] = (x[i] - xc[i]) + box[i]; }
    const double CH = getenv("CHUNK") ? atof(getenv("CHUNK")) : 2.0;
    double done = 0.0;
    SetT S(xc, C, r0, rr);
    for (int c = 0; c < nch; ++c) {
      for (int i = 0; i < npc; ++i) vfb.setParameter(pcn[i].c_str(), chp[c][i]);
      IOdeSolver sb(vfb, order);
      sb.setStep(chh[c]); sb.turnOffStepControl();
      ITimeMap tm(sb);
      ITimeMap::SolutionCurve sol(done);
      tm(interval(done + chlen[c]), S, sol);
      if (nmon > 0) monitor("B", sol, done, done + chlen[c], dsub, mon, 32.0, -1.0);
      done += chlen[c];
      if (getenv("SETDEBUG")) {
        IVector xx = S.get_x(), cr = S.get_C() * S.get_r0(), br = S.get_B() * S.get_r();
        double ax = 0, ac = 0, ab = 0;
        for (int i = 6; i < n; ++i) { ax = std::max(ax, (double)rightBound(abs(xx[i]))); ac = std::max(ac, (double)rightBound(abs(cr[i]))); ab = std::max(ab, (double)rightBound(abs(br[i]))); }
        cout << "SET sigma=" << 32.0 - done << " |x|=" << ax << " |C r0|=" << ac << " |B r|=" << ab << "  X=" << xx[4] << " E=" << xx[5] << " Br_E=" << br[5] << "\n";
      }
    }
    if (std::fabs(done - tauB) > 1e-9) { cout << "RESULT FAIL chunk lengths\n"; return 1; }
    IVector xs = S.get_x(), r0s = S.get_r0(); IMatrix Cs = S.get_C(); IVector Br = S.get_B() * S.get_r();
    IVector hullB = xs + Cs * r0s + Br;
    cout << "STAGEB";
    for (int i = 0; i < n; ++i) cout << " " << hullB[i].leftBound() << " " << hullB[i].rightBound();
    cout << "\n";
    // ---------------- stage C set: profile rows from the step-1 set, all other rows boxed
    IVector xn(n), rn(n), r0n(n); IMatrix Cn(n, n);
    for (int i = 0; i < 4; ++i) { xn[i] = interval(px[i].mid().leftBound()); rn[i] = px[i] - xn[i];
                                  for (int j = 0; j < 4; ++j) Cn[i][j] = pC[i][j]; r0n[i] = pr0[i]; }
    for (int i = 4; i < n; ++i) { xn[i] = interval(hullB[i].mid().leftBound()); rn[i] = hullB[i] - xn[i]; }
    SetT SC(xn, Cn, r0n, rn);
    if (getenv("HC")) { solC.setStep(atof(getenv("HC"))); solC.turnOffStepControl(); }
    ITimeMap tmc(solC);
    ITimeMap::SolutionCurve solc(0.0);
    tmc(interval(TLO), SC, solc);
    if (nmon > 0) monitor("C", solc, 0.0, TLO, dsub, mon, 8.0, -1.0);
    SetT SD = SC;
    ITimeMap::SolutionCurve solw(TLO);
    tmc(interval(THI), SD, solw);
    IVector K = solw(interval(TLO, THI));
    cout << "TIP";
    for (int i = 0; i < n; ++i) cout << " " << K[i].leftBound() << " " << K[i].rightBound();
    cout << "\n";
    if (nmon > 0) {
      cout << "MON W " << 8.0 - TLO << " " << 8.0 - THI;
      for (int j : mon) cout << " " << K[j].leftBound() << " " << K[j].rightBound();
      cout << "\n";
    }
    cout << "RESULT OK\n";
  } catch (std::exception& e) {
    cout << "RESULT EXCEPTION: " << e.what() << "\n";
    return 1;
  }
  return 0;
}
