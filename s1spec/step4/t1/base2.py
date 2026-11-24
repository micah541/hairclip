"""Access to a stored float base (run2 format): slabs_<tag>.pkl -> per-slab Rep, Y(tau), k(tau), tau(k)."""
import pickle, numpy as np
from rep import Rep, Mesh, K0, DS1
from ops import Dyn
from radau import Radau

class Base:
    def __init__(self, tag):
        z = pickle.load(open(f'/home/claude/s1spec/step4/t1/slabs_{tag}.pkl', 'rb'))
        self.cfg = cfg = z['cfg']; self.slabs = z['slabs']; self.RD = Radau(cfg['s'])
        m0 = Mesh(JT=cfg['JT'], h0=cfg['h0'], cT=cfg['cT'], C=np.zeros(1))
        self.tau0 = np.array([q['tau0'] for q in self.slabs]); self.dt = np.array([q['dt'] for q in self.slabs])
        self._R = {}
        self.lnlam_end = np.array([q['nodes'][-1, -3] for q in self.slabs])
    def rep(self, i):
        q = self.slabs[i]; key = len(q['C'])
        if i not in self._R:
            cfg = self.cfg; m = Mesh(JT=cfg['JT'], h0=cfg['h0'], cT=cfg['cT'], C=q['C'])
            R = Rep(cfg['P'], m); self._R[i] = (R, Dyn(R))
        return self._R[i]
    def slab_of(self, tau):
        i = int(np.searchsorted(self.tau0, tau, side='right') - 1); return max(0, min(i, len(self.tau0) - 1))
    def Y(self, tau, i=None):
        i = self.slab_of(tau) if i is None else i; th = (tau - self.tau0[i]) / self.dt[i]
        return (self.RD.interp_matrix(np.array([th])) @ self.slabs[i]['nodes'])[0]
    def Ydot(self, tau, i=None):
        i = self.slab_of(tau) if i is None else i; th = (tau - self.tau0[i]) / self.dt[i]
        return (self.RD.deriv_matrix(np.array([th])) @ self.slabs[i]['nodes'])[0] / self.dt[i]
    def k(self, tau): return K0 * np.exp(-self.Y(tau)[-3])
    def tau_of_k(self, kk):
        lk = np.log(K0 / kk); i = int(np.argmax(self.lnlam_end <= lk))
        if self.lnlam_end[i] > lk: raise ValueError('k beyond the stored base')
        lo, hi = self.tau0[i], self.tau0[i] + self.dt[i]
        for _ in range(60):
            mid = 0.5 * (lo + hi)
            if self.Y(mid, i)[-3] > lk: lo = mid
            else: hi = mid
        return 0.5 * (lo + hi), i
    def kappa(self, tau, xs, m=0, i=None):
        i = self.slab_of(tau) if i is None else i; R, D = self.rep(i); Y = self.Y(tau, i)
        return R.eval_at(Y[:R.n], Y[R.n], np.asarray(xs, float), m)
