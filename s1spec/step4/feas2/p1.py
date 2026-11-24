"""P1: float prototype of a Galerkin + tail coercivity certificate for  -sym_rho(L_Gamma) - omega  on the rho-complement of 3 constraint functions.
Reference operator: A_ref = -sym_rho(L_S1) (self-adjoint), basis e_1..e_N = its first N rho-orthonormal eigenvectors (nu_1<=nu_2<=...).
u = sum_{j<=N} c_j e_j + v,  v ⟂ e_1..e_N,  c_{1..3} = -<u, e_i - phi_i> (constraints phi_i = chi e_i, i<=3).
 -S_G(u,u) = sum nu_j c_j^2 + E(v) + Delta(Pu,Pu) + 2 Delta_s(Pu,v) + Delta(v,v),   E(v) = -S_ref(v,v) >= nu_{N+1}|v|^2,  Delta(a,b) = -<dA a,b>_rho (sym),  dA = L_G - L_ref (2nd order)
 tail:  Delta(v,v) >= -(d2 |v''||v| + d1 |v'||v| + d0 |v|^2)  with Garding  |v''|^2 <= 2E + 2 C_G |v|^2,  |v'|^2 <= 2|v||v''| + m1^2|v|^2
 cross: 2|Delta_s(Pu,v)| <= |v| (|dA Pu| + |dA^* Pu|)
 =>  Q(u) >= c^T A c + t x^2 - x sqrt(2 c^T (G+G*) c),  A = diag(nu-omega)+D,  t = min_{E>=nu_{N+1}}[E - d2 sqrt(2E+2C_G) - d1 sqrt(2 sqrt(2E+2C_G)+m1^2) - d0] - omega
 => closes iff  A_bb - (G+G*)_bb/(2t) >= 0 (b = indices 4..N).   Float only: FEM eigenpairs treated as exact; sup-norms over a fine grid."""
import sys, json, numpy as np, scipy.linalg as sla
sys.dont_write_bytecode = True
sys.path.insert(0, '/home/claude/s1spec/step4/feas2')
from f3geom import *; from f3fem import *
from scipy.interpolate import CubicSpline
import sympy as sp_

# ---- analytic rho and derivatives up to order 4
_x = sp_.symbols('x'); _rho = sp_.exp(-(1 + _x**2) ** sp_.Rational(2, 3) / 4)
_RD = [sp_.lambdify(_x, sp_.diff(_rho, _x, n), 'numpy') for n in range(5)]
def rhod(x, n): return _RD[n](np.asarray(x, float)) * np.ones_like(x, dtype=float)

def smoothstep_cut(x, a, b):
    """1 for x<=a, 0 for x>=b, quintic C^2 in between"""
    t = np.clip((x - a) / (b - a), 0, 1); return 1 - t**3 * (10 - 15 * t + 6 * t * t)

class Setup:
    def __init__(self, Sigma=35.0, h=0.1, nmax=70):
        self.mesh = Mesh(Sigma, h); m = self.mesh; self.M = m.mass(); self.Sigma = Sigma
        self.xg = np.arange(0, Sigma + 1e-9, 0.01)
        self.sp1 = s1_splines(Smax=Sigma + 5, dx=0.01)
        K1 = m.form(coeffs(m, self.sp1, use_tau=False)); self.K1 = K1
        w, V = eig_sym(K1, self.M); self.nu = -w[:nmax]; self.V = V[:, :nmax]
        ne = m.ne; self.idx = (2 * np.arange(ne)[:, None] + np.arange(4)[None, :])
        self.W = m.wq * m.rho
        self.x = m.xq
    def funcs(self, c):
        """values, first, second derivative of FEM function with free-dof vector c at the quadrature points"""
        m = self.mesh; val = np.zeros(m.ndof); val[m.free] = c; loc = val[self.idx]
        return [np.einsum('ei,iq->eq', loc, B) for B in m.B]
    def coeff_grid(self, sp, tau, x=None):
        x = self.xg if x is None else x; k, k1, k2, q = sp['kap'](x), sp['kap1'](x), sp['kap2'](x), sp['q'](x)
        return np.array([-k**2, -5 * k * k1 - q / 4 + (sp['tau'](x) if tau else 0.0), 0.25 - 3 * k1**2 - 4 * k * k2])
    def ref_constants(self):
        """S1: m1 = sup|rho'/rho|, Bbar = sup (-B2/rho)_+, cbar = sup (c0/rho)_+ over the whole domain;  C_G = 2 Bbar^2 + Bbar m1^2 + cbar"""
        x = self.xg; b = self.coeff_grid(self.sp1, False); r = [rhod(x, n) for n in range(5)]
        bs = [CubicSpline(x, b[i]) for i in range(3)]
        b2, b2d, b2dd = b[0], bs[0](x, 1), bs[0](x, 2); b1, b1d = b[1], bs[1](x, 1); b0 = b[2]
        m1 = np.max(np.abs(r[1] / r[0]))
        B2 = b2 - 2 * r[2] / r[0]
        phi_dd = (r[2] * b2 + 2 * r[1] * b2d + r[0] * b2dd) - r[4]               # (rho b2 - rho'')''
        c0 = 0.5 * phi_dd - 0.5 * (r[1] * b1 + r[0] * b1d) + r[0] * b0
        Bbar = max(np.max(-B2), 0.0); cbar = max(np.max(c0 / r[0]), 0.0)
        return dict(m1=m1, Bbar=Bbar, cbar=cbar, CG=2 * Bbar**2 + Bbar * m1**2 + cbar)
    def diff_coeffs(self, spG, tauG=True, Rc=None, taper=4.0, xmaxG=40.0):
        """Delta b_j = (b_j^G - b_j^S1) * chi_Rc (chi = 1 on [0,Rc-taper], 0 beyond Rc), as splines (values and 2 derivatives)"""
        x = self.xg; d = np.zeros((3, len(x))); mk = x <= xmaxG
        d[:, mk] = self.coeff_grid(spG, tauG, x[mk]) - self.coeff_grid(self.sp1, False, x[mk])      # Gamma data only on [0, xmaxG]; difference set to 0 beyond (taper must end before)
        chi = np.ones_like(x) if Rc is None else smoothstep_cut(x, Rc - taper, Rc)
        d = d * chi[None, :]; return [CubicSpline(x, d[i]) for i in range(3)], d

def certificate(st, spG, tauG=True, Rc=None, omega=0.05, Ns=tuple(range(3, 81)), Rphi=(20.0, 28.0), const=None, verbose=True, ncon=3):
    m = st.mesh; x = st.x; const = const or st.ref_constants()
    dsp, dgrid = st.diff_coeffs(spG, tauG, Rc)
    d2, d1, d0 = [float(np.max(np.abs(dgrid[i]))) for i in range(3)]
    # coefficient differences at quadrature points
    D2, D1, D0 = [dsp[i](x) for i in range(3)]; D2p, D2pp, D1p = dsp[0](x, 1), dsp[0](x, 2), dsp[1](x, 1)
    r1, r2 = m.rho1 / m.rho, m.rho2 / m.rho
    nmax = max(Ns) + 1; W = st.W
    Fs = [st.funcs(st.V[:, j]) for j in range(nmax)]
    X = np.array([D2 * f[2] + D1 * f[1] + D0 * f[0] for f in Fs])                      # dA e_j
    Y = np.array([D2 * f[2] + (2 * D2p + 2 * r1 * D2 - D1) * f[1] + (D2pp + 2 * r1 * D2p + r2 * D2 - D1p - r1 * D1 + D0) * f[0] for f in Fs])   # dA^* e_j
    E0 = np.array([f[0] for f in Fs])
    Gm = np.einsum('jeq,ieq,eq->ij', X, E0, W)          # <dA e_j, e_i>
    Gs = np.einsum('jeq,ieq,eq->ij', Y, E0, W)          # <dA^* e_j, e_i> (should be Gm^T)
    adj_err = float(np.max(np.abs(Gs - Gm.T)) / max(np.max(np.abs(Gm)), 1e-300))
    D = -0.5 * (Gm + Gm.T)
    GG = np.einsum('ieq,jeq,eq->ij', X, X, W); GS = np.einsum('ieq,jeq,eq->ij', Y, Y, W)
    # constraint leak eps_i = |e_i - phi_i|_rho
    chi_phi = smoothstep_cut(m.xq, Rphi[0], Rphi[1]); eps = [float(np.sqrt(np.sum(W * ((1 - chi_phi) * Fs[i][0])**2))) for i in range(3)]
    out = []
    for N in Ns:
        nuN1 = st.nu[N]; nu = st.nu[:N]
        def f(E): y = np.sqrt(2 * E + 2 * const['CG']); return E - d2 * y - d1 * np.sqrt(2 * y + const['m1']**2) - d0
        Eg = nuN1 * np.exp(np.linspace(0, np.log(1e4), 4000)) if nuN1 > 0 else np.linspace(1e-3, 1e4, 4000)
        fv = f(Eg); iE = int(np.argmin(fv)); t = float(fv[iE] - omega)
        Em = Eg[iE]; y = np.sqrt(2 * Em + 2 * const['CG'])
        pen = dict(d2=d2 * y, d1=d1 * np.sqrt(2 * y + const['m1']**2), d0=d0)
        b = np.arange(ncon, N); A = np.diag(nu - omega) + D[:N, :N]
        Abb = A[np.ix_(b, b)]; Gb = (GG + GS)[:N, :N][np.ix_(b, b)]
        row = dict(N=N, nuN1=float(nuN1), t=t, pen_d2=float(pen['d2']), pen_d1=float(pen['d1']), pen_d0=float(pen['d0']), Estar=float(Em))
        if len(b):
            lam_A = float(np.linalg.eigvalsh(Abb)[0]); row['lamA'] = lam_A
            if t > 0:
                # margin m*: largest m in [0,t) with lam_min(A - m - G/(2(t-m))) >= 0
                red = lambda mm: np.linalg.eigvalsh(Abb - mm * np.eye(len(b)) - Gb / (2 * (t - mm)))[0]
                r0 = float(red(0.0)); row['lamR'] = r0
                if r0 >= 0:
                    lo, hi = 0.0, t * (1 - 1e-9)
                    for _ in range(60):
                        mid = 0.5 * (lo + hi)
                        if red(mid) >= 0: lo = mid
                        else: hi = mid
                    row['margin'] = float(min(lo, t))
                else: row['margin'] = r0
            else: row['lamR'] = float('nan'); row['margin'] = t
        else:
            row['lamA'] = float('nan'); row['lamR'] = float('nan'); row['margin'] = t
        # constraint leak (3 constrained coefficients are not exactly 0): subtract eps*(2|A_ab| + eps|A_aa|) crude
        if len(b) and t > 0:
            Aab = np.linalg.norm(A[:ncon, ncon:N], 2) if N > ncon else 0.0; Aaa = np.linalg.norm(A[:ncon, :ncon], 2); ee = max(eps)
            row['leak'] = float(ee * (2 * Aab + ee * Aaa) + 0)
        else: row['leak'] = float(max(eps) ** 2 * np.linalg.norm(np.diag(st.nu[:ncon] - omega) + D[:ncon, :ncon], 2))
        row['margin_net'] = row['margin'] - row['leak']
        out.append(row)
    info = dict(d2=d2, d1=d1, d0=d0, adj_err=adj_err, eps=eps, const={k: float(v) for k, v in const.items()}, Rc=Rc)
    return info, out

def show(name, info, rows):
    print(f"== {name}: Rc={info['Rc']} d2={info['d2']:.3e} d1={info['d1']:.3e} d0={info['d0']:.3e}  adj_err={info['adj_err']:.1e} eps_phi={max(info['eps']):.1e}  consts {info['const']}", flush=True)
    print("   N  nu_{N+1}     t     pen(d2,d1,d0)           lam(A_bb)   lam(R)    margin   net", flush=True)
    for r in rows:
        if r['N'] not in (3,4,6,8,10,12,16,20,24,30,40,50,60,70,80): continue
        print(f"  {r['N']:3d} {r['nuN1']:8.3f} {r['t']:8.3f}  ({r['pen_d2']:.3f},{r['pen_d1']:.3f},{r['pen_d0']:.3f})  {r['lamA']:9.3f} {r['lamR']:9.3f} {r['margin']:9.4f} {r['margin_net']:9.4f}", flush=True)
