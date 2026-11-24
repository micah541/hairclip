// Validated inward integration, version 2 (time-dependent vector fields; CAPD Lohner sets).
// Stage B: sigma 32 -> 8 = CAPD time 0..24 (sigma = 32 - t), state = nB non-profile rows; profile terms enter as
//          per-chunk interval parameters (step-1 tail theorem).  One Lohner set for the whole stage.
// Stage C: sigma 8 -> tip = CAPD time 0..tau_tip (sigma = 8 - t); state = 4 profile rows (tight step-1 set) +
//          the nB rows (boxed from the stage-B hull).  Enclosure over the tip window [TLO, THI] is printed.
// Monitoring (optional): enclosures of selected stage-C components (indices in the stage-C state) over
//          subintervals, and of the corresponding stage-B components (index - 4).
// Input:  nB  np  {name lo hi}^np  VFB-line  VFC-line
//         nch npc {pname}^npc  {len h {lo hi}^npc}^nch
//         x: nB (lo hi)  ncol  C: nB x ncol  r0: ncol (lo hi)  box: nB (lo hi)
//         profile set: 4 (lo hi), 4x4, 4 (lo hi)   TLO THI   nmon idx...  dsub
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include "capd/capdlib.h"
using namespace std;
using namespace capd;

static void monitor(const char* tag, ITimeMap::SolutionCurve& sol, double t0, double t1, double dsub,
                    const vector<int>& mon, int off, double sig0) {
  int N = std::max(1, (int)std::ceil((t1 - t0) / dsub));
  for (int i = 0; i < N; ++i) {
    double a = t0 + (t1 - t0) * i / N, b = (i + 1 == N) ? t1 : t0 + (t1 - t0) * (i + 1) / N;
    if (b > t1) b = t1; if (a < t0) a = t0;
    IVector v = sol(interval(a, b));
    cout << "MON " << tag << " " << sig0 - a << " " << sig0 - b;
    for (int j : mon) cout << " " << v[j - off].leftBound() << " " << v[j - off].rightBound();
    cout << "\n";
  }
}

int main(int argc, char* argv[]) {
  cout.precision(17);
  if (argc < 2) { cerr << "usage: hsrun2 input [order]\n"; return 2; }
  int order = (argc > 2) ? atoi(argv[2]) : 20;
  ifstream in(argv[1]);
  int nB, np; in >> nB >> np;
  int n = nB + 4;
  vector<string> pn(np); vector<interval> pv(np);
  for (int i = 0; i < np; ++i) { double lo, hi; in >> pn[i] >> lo >> hi; pv[i] = interval(lo, hi); }
  string line; getline(in, line);
  string VFB, VFC; getline(in, VFB); getline(in, VFC);
  int nch, npc; in >> nch >> npc;
  vector<string> pcn(npc); for (int i = 0; i < npc; ++i) in >> pcn[i];
  vector<double> chlen(nch), chh(nch); vector<vector<interval>> chp(nch, vector<interval>(npc));
  for (int c = 0; c < nch; ++c) { in >> chlen[c] >> chh[c];
    for (int i = 0; i < npc; ++i) { double lo, hi; in >> lo >> hi; chp[c][i] = interval(lo, hi); } }
  IVector x(nB), box(nB); for (int i = 0; i < nB; ++i) { double lo, hi; in >> lo >> hi; x[i] = interval(lo, hi); }
  int ncol; in >> ncol;
  IMatrix C(nB, nB); IVector r0(nB);
  for (int i = 0; i < nB; ++i) for (int j = 0; j < ncol; ++j) { double v; in >> v; C[i][j] = interval(v); }
  for (int j = 0; j < ncol; ++j) { double lo, hi; in >> lo >> hi; r0[j] = interval(lo, hi); }
  for (int i = 0; i < nB; ++i) { double lo, hi; in >> lo >> hi; box[i] = interval(lo, hi); }
  IVector px(4), pr0(4); IMatrix pC(4, 4);
  for (int i = 0; i < 4; ++i) { double lo, hi; in >> lo >> hi; px[i] = interval(lo, hi); }
  for (int i = 0; i < 4; ++i) for (int j = 0; j < 4; ++j) { double v; in >> v; pC[i][j] = interval(v); }
  for (int i = 0; i < 4; ++i) { double lo, hi; in >> lo >> hi; pr0[i] = interval(lo, hi); }
  double TLO, THI; in >> TLO >> THI;
  int nmon; in >> nmon; vector<int> mon(nmon); for (int i = 0; i < nmon; ++i) in >> mon[i];
  double dsub; in >> dsub;
  if (!in) { cerr << "bad input\n"; return 2; }
  try {
    IMap vfb(VFB.c_str()), vfc(VFC.c_str());
    for (int i = 0; i < np; ++i) { vfb.setParameter(pn[i].c_str(), pv[i]); vfc.setParameter(pn[i].c_str(), pv[i]); }
    // ---------------- stage B
    IVector xc(nB);
    for (int i = 0; i < nB; ++i) xc[i] = interval(x[i].mid().leftBound());
    C0Rect2Set S(xc, C, r0, (x - xc) + box);
    double done = 0.0;
    bool monB = (nmon > 0) && !getenv("NOMONB");
    for (int c = 0; c < nch; ++c) {
      for (int i = 0; i < npc; ++i) vfb.setParameter(pcn[i].c_str(), chp[c][i]);
      IOdeSolver sb(vfb, order);
      sb.setStep(chh[c]); sb.turnOffStepControl();
      ITimeMap tm(sb);
      if (monB) {
        double t0c = S.getCurrentTime().rightBound();
        ITimeMap::SolutionCurve sol(t0c);
        tm(interval(done + chlen[c]), S, sol);
        double t1c = std::min(done + chlen[c], (double)sol.getRightDomain());
        monitor("B", sol, t0c, t1c, dsub, mon, 4, 32.0);
      } else {
        tm(interval(done + chlen[c]), S);
      }
      done += chlen[c];
      if (getenv("SETDEBUG")) {
        IVector xx = S.get_x(), cr = S.get_C() * S.get_r0(), br = S.get_B() * S.get_r();
        double ax = 0, ac = 0, ab = 0;
        for (int i = 0; i < nB; ++i) { ax = std::max(ax, (double)rightBound(abs(xx[i]))); ac = std::max(ac, (double)rightBound(abs(cr[i]))); ab = std::max(ab, (double)rightBound(abs(br[i]))); }
        cout << "SET sigma=" << 32.0 - done << " |x|=" << ax << " |C r0|=" << ac << " |B r|=" << ab << "\n";
      }
    }
    IVector hullB = S.get_x() + S.get_C() * S.get_r0() + S.get_B() * S.get_r();
    cout << "STAGEB";
    for (int i = 0; i < nB; ++i) cout << " " << hullB[i].leftBound() << " " << hullB[i].rightBound();
    cout << "\n";
    // ---------------- stage C
    IVector xn(n), rn(n), r0n(n); IMatrix Cn(n, n);
    for (int i = 0; i < 4; ++i) { xn[i] = interval(px[i].mid().leftBound()); rn[i] = px[i] - xn[i];
                                  for (int j = 0; j < 4; ++j) Cn[i][j] = pC[i][j]; r0n[i] = pr0[i]; }
    for (int i = 0; i < nB; ++i) { xn[4 + i] = interval(hullB[i].mid().leftBound()); rn[4 + i] = hullB[i] - xn[4 + i]; }
    C0Rect2Set SC(xn, Cn, r0n, rn);
    IOdeSolver solC(vfc, order);
    if (getenv("HC")) { solC.setStep(atof(getenv("HC"))); solC.turnOffStepControl(); }
    ITimeMap tmc(solC);
    ITimeMap::SolutionCurve solc(0.0);
    tmc(interval(TLO), SC, solc);
    if (nmon > 0) monitor("C", solc, 0.0, TLO, dsub, mon, 0, 8.0);
    C0Rect2Set SD = SC;
    ITimeMap::SolutionCurve solw(TLO);
    tmc(interval(THI), SD, solw);
    IVector K = solw(interval(TLO, THI));
    cout << "TIP";
    for (int i = 0; i < n; ++i) cout << " " << K[i].leftBound() << " " << K[i].rightBound();
    cout << "\n";
    cout << "RESULT OK\n";
  } catch (std::exception& e) {
    cout << "RESULT EXCEPTION: " << e.what() << "\n";
    return 1;
  }
  return 0;
}
