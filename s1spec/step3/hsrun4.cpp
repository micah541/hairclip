// Validated inward integration, version 3.
//  * the tip parameter st is a STATE (d st/dt = 0), so its uncertainty is carried by the Lohner set (no
//    per-step re-wrapping of an interval parameter);
//  * stage B -> C switch keeps the full stage-B doubleton representation x + C r0 + B r (no boxing);
//  * stage C runs in rescaled time t in [0, T0] with sigma = 8 - c t, c = (8 - st)/T0 (built into VFC by the
//    generator), so that t = T0 is exactly the tip sigma = st for every st in the set;
//  * output at t = T0: the full doubleton representation (x, C, r0, B, r) for mean-value post-processing.
// Stage B: sigma 32 -> 8 = CAPD time 0..24 (sigma = 32 - t), state = nB rows; profile terms enter as per-chunk
//          interval parameters (step-1 tail theorem).
// Input:  nB  np  {name lo hi}^np  VFB-line  VFC-line
//         nch npc {pname}^npc  {len h {lo hi}^npc}^nch
//         x: nB (lo hi)  ncol  C: nB x ncol  r0: ncol (lo hi)  box: nB (lo hi)
//         profile set: 4 (lo hi), 4x4, 4 (lo hi)   T0   nmon idx...  dsub
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include "capd/capdlib.h"
using namespace std;
using namespace capd;

static void monitor(const char* tag, ITimeMap::SolutionCurve& sol, double t0, double t1, double dsub,
                    const vector<int>& mon, int off) {
  int N = std::max(1, (int)std::ceil((t1 - t0) / dsub));
  for (int i = 0; i < N; ++i) {
    double a = t0 + (t1 - t0) * i / N, b = (i + 1 == N) ? t1 : t0 + (t1 - t0) * (i + 1) / N;
    if (b > t1) b = t1; if (a < t0) a = t0;
    IVector v = sol(interval(a, b));
    cout << "MON " << tag << " " << a << " " << b;
    for (int j : mon) cout << " " << v[j - off].leftBound() << " " << v[j - off].rightBound();
    cout << "\n";
  }
}

static void monitorAll(const char* tag, ITimeMap::SolutionCurve& sol, double t0, double t1, double dsub) {
  int N = std::max(1, (int)std::ceil((t1 - t0) / dsub));
  for (int i = 0; i < N; ++i) {
    double a = t0 + (t1 - t0) * i / N, b = (i + 1 == N) ? t1 : t0 + (t1 - t0) * (i + 1) / N;
    IVector v = sol(interval(a, b));
    cout << tag << " " << a << " " << b;
    for (int j = 0; j < v.dimension(); ++j) cout << " " << v[j].leftBound() << " " << v[j].rightBound();
    cout << "\n";
  }
}
static void pv(const char* tag, const IVector& v) {
  cout << tag;
  for (int i = 0; i < v.dimension(); ++i) cout << " " << v[i].leftBound() << " " << v[i].rightBound();
  cout << "\n";
}
static void pm(const char* tag, const IMatrix& M) {
  cout << tag << " " << M.numberOfRows() << " " << M.numberOfColumns();
  for (int i = 0; i < (int)M.numberOfRows(); ++i)
    for (int j = 0; j < (int)M.numberOfColumns(); ++j) cout << " " << M[i][j].leftBound() << " " << M[i][j].rightBound();
  cout << "\n";
}


template <class S> static auto setRF_impl(S& s, double f, int) -> decltype(s.setFactor(f), void()) { s.setFactor(f); }
template <class S> static void setRF_impl(S&, double, long) {}
template <class S> static void setRF(S& s, double f) { setRF_impl(s, f, 0); }

template <class SetT>
static void stageC(IMap& vfc, IMap& vfd, int order, int n, double T0, double TS, double RF, const IVector& xn, const IMatrix& Cn,
                   const IVector& r0n, const IMatrix& Bn, const IVector& rn) {
  SetT SC(xn, Cn, r0n, Bn, rn, interval(0.0));
  setRF(SC, RF);
  double dm = getenv("MONALL") ? atof(getenv("MONALL")) : 0.0;
  {
    IOdeSolver solC(vfc, order);
    if (getenv("HC")) { solC.setStep(atof(getenv("HC"))); solC.turnOffStepControl(); }
    ITimeMap tmc(solC);
    if (dm > 0) { ITimeMap::SolutionCurve c(0.0); tmc(interval(TS), SC, c); monitorAll("MC", c, 0.0, TS, dm); }
    else tmc(interval(TS), SC);
  }
  {
    IOdeSolver solD(vfd, order);
    if (getenv("HC")) { solD.setStep(atof(getenv("HC"))); solD.turnOffStepControl(); }
    ITimeMap tmd(solD);
    if (dm > 0) { ITimeMap::SolutionCurve c(TS); tmd(interval(T0), SC, c); monitorAll("MC", c, TS, T0, dm); }
    else tmd(interval(T0), SC);
  }
  IVector hull = SC.get_x() + SC.get_C() * SC.get_r0() + SC.get_B() * SC.get_r();
  pv("TIP", hull);
  pv("TX", SC.get_x()); pm("TC", SC.get_C()); pv("TR0", SC.get_r0()); pm("TB", SC.get_B()); pv("TR", SC.get_r());
}

int main(int argc, char* argv[]) {
  cout.precision(17);
  if (argc < 2) { cerr << "usage: hsrun3 input [order]\n"; return 2; }
  int order = (argc > 2) ? atoi(argv[2]) : 20;
  ifstream in(argv[1]);
  int nB, np; in >> nB >> np;
  int n = nB + 4;
  vector<string> pn(np); vector<interval> pvv(np);
  for (int i = 0; i < np; ++i) { double lo, hi; in >> pn[i] >> lo >> hi; pvv[i] = interval(lo, hi); }
  string line; getline(in, line);
  string VFB, VFC, VFD; getline(in, VFB); getline(in, VFC); getline(in, VFD);
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
  double T0, TS; in >> T0 >> TS;
  int nmon; in >> nmon; vector<int> mon(nmon); for (int i = 0; i < nmon; ++i) in >> mon[i];
  double dsub; in >> dsub;
  if (!in) { cerr << "bad input\n"; return 2; }
  try {
    IMap vfb(VFB.c_str()), vfc(VFC.c_str()), vfd(VFD.c_str());
    for (int i = 0; i < np; ++i) { vfb.setParameter(pn[i].c_str(), pvv[i]); vfc.setParameter(pn[i].c_str(), pvv[i]); vfd.setParameter(pn[i].c_str(), pvv[i]); }
    // ---------------- stage B
    IVector xc(nB);
    for (int i = 0; i < nB; ++i) xc[i] = interval(x[i].mid().leftBound());
    C0Rect2Set S(xc, C, r0, (x - xc) + box);
    double RF = getenv("REORG") ? atof(getenv("REORG")) : 1e300;
    S.setFactor(RF);
    double done = 0.0;
    bool monB = (nmon > 0) && !getenv("NOMONB");
    for (int c = 0; c < nch; ++c) {
      for (int i = 0; i < npc; ++i) vfb.setParameter(pcn[i].c_str(), chp[c][i]);
      IOdeSolver sb(vfb, order);
      { double bs = getenv("BSUB") ? atof(getenv("BSUB")) : 1.0; sb.setStep(chh[c] / bs); sb.turnOffStepControl(); }
      ITimeMap tm(sb);
      if (getenv("MONALL")) {
        double t0c = S.getCurrentTime().rightBound();
        ITimeMap::SolutionCurve sol(t0c);
        tm(interval(done + chlen[c]), S, sol);
        monitorAll("MB", sol, t0c, std::min(done + chlen[c], (double)sol.getRightDomain()), atof(getenv("MONALL")));
      } else if (monB) {
        double t0c = S.getCurrentTime().rightBound();
        ITimeMap::SolutionCurve sol(t0c);
        tm(interval(done + chlen[c]), S, sol);
        double t1c = std::min(done + chlen[c], (double)sol.getRightDomain());
        monitor("B", sol, t0c, t1c, dsub, mon, 4);
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
    pv("STAGEB", hullB);
    // ---------------- stage C: keep the doubleton representation
    IVector xB = S.get_x(), r0B = S.get_r0(), rB = S.get_r();
    IMatrix CB = S.get_C(), BB = S.get_B();
    IVector xn(n), rn(n), r0n(n); IMatrix Cn(n, n), Bn(n, n);
    for (int i = 0; i < n; ++i) for (int j = 0; j < n; ++j) { Cn[i][j] = 0.0; Bn[i][j] = 0.0; }
    for (int i = 0; i < 4; ++i) {
      xn[i] = interval(px[i].mid().leftBound()); rn[i] = px[i] - xn[i]; r0n[i] = pr0[i]; Bn[i][i] = 1.0;
      for (int j = 0; j < 4; ++j) Cn[i][j] = pC[i][j];
    }
    for (int i = 0; i < nB; ++i) {
      xn[4 + i] = xB[i]; rn[4 + i] = rB[i]; r0n[4 + i] = r0B[i];
      for (int j = 0; j < nB; ++j) { Cn[4 + i][4 + j] = CB[i][j]; Bn[4 + i][4 + j] = BB[i][j]; }
    }
    const char* cs = getenv("CSET");
    string cset = cs ? cs : "rect2";
    if (cset == "rect2") stageC<C0Rect2Set>(vfc, vfd, order, n, T0, TS, RF, xn, Cn, r0n, Bn, rn);
    else if (cset == "trip") stageC<C0TripletonSet>(vfc, vfd, order, n, T0, TS, RF, xn, Cn, r0n, Bn, rn);
    else { cerr << "bad CSET\n"; return 2; }
    cout << "RESULT OK\n";
  } catch (std::exception& e) {
    cout << "RESULT EXCEPTION: " << e.what() << "\n";
    return 1;
  }
  return 0;
}
