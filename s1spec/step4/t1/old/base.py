"""Access to a stored float base: slabs_<tag>.npz -> Y(tau), k(tau), tau(k), geometry."""
import json, numpy as np
from rep import Rep, Mesh, K0, DS1
from ops import Dyn
from radau import Radau


class Base:
    def __init__(self, tag, path=None):
        z = np.load(path or f'/home/claude/s1spec/step4/t1/slabs_{tag}.npz', allow_pickle=True)
        self.cfg = cfg = json.loads(str(z['cfg']))
        self.R = Rep(cfg['P'], Mesh(JT=cfg['JT'], h0=cfg['h0'], cT=cfg['cT'], NO=cfg['NO'], c=cfg['c'], beta=cfg['beta']))
        self.D = Dyn(self.R); self.RD = Radau(cfg['s'])
        self.tau0 = z['tau0']; self.dt = z['dt']; self.nodes = z['nodes']; self.n = self.R.n
        self.tau_end = self.tau0[-1] + self.dt[-1]
        self.lnlam_end = self.nodes[:, -1, self.n + 1]

    def slab_of(self, tau):
        i = int(np.searchsorted(self.tau0, tau, side='right') - 1); return max(0, min(i, len(self.tau0) - 1))

    def Y(self, tau):
        i = self.slab_of(tau); th = (tau - self.tau0[i]) / self.dt[i]
        W = self.RD.interp_matrix(np.array([th]))
        return (W @ self.nodes[i])[0]

    def Ydot(self, tau):
        i = self.slab_of(tau); th = (tau - self.tau0[i]) / self.dt[i]
        Wd = self.RD.deriv_matrix(np.array([th]))
        return (Wd @ self.nodes[i])[0] / self.dt[i]

    def k(self, tau): return K0 * np.exp(-self.Y(tau)[self.n + 1])

    def tau_of_k(self, kk):
        lk = np.log(K0 / kk)            # target ln lam
        i = int(np.argmax(self.lnlam_end <= lk))
        if self.lnlam_end[i] > lk: raise ValueError('k beyond the stored base')
        lo, hi = self.tau0[i], self.tau0[i] + self.dt[i]
        for _ in range(60):
            mid = 0.5 * (lo + hi)
            if self.Y(mid)[self.n + 1] > lk: lo = mid
            else: hi = mid
        return 0.5 * (lo + hi)

    def kappa(self, tau, xs, m=0):
        Y = self.Y(tau); return self.R.eval_at(Y[:self.n], Y[self.n], np.asarray(xs, float), m)
