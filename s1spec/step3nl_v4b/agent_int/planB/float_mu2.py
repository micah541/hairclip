from float_mu import *
c1t = 2*D(c2) + c1; c0t = D(D(c2)) + 2*D(c1) + c0
e1 = sp.expand(D(D(c1)) + 2*D(c0)); e0 = sp.expand(D(D(c0)))
def forms(a2_, a1_, a0_):
    A = 2*(P1s**2 - P2s) - a2_
    V = a0_ - (D(a1_) - P1s*a1_)/2 + rr2(a2_)/2 - rr2(P1s**2 - P2s)/2
    return sp.expand(-A**2/2 - rr2(A)/2 - V)
mu2 = forms(c2, c1t, c0t); f2 = sp.lambdify((P1s, P2s, P3s, P4s, k, m, p, q), mu2, 'numpy')
g = sp.log(1 + (x/60)**4); gd = sp.lambdify(x, [sp.diff(g, x, j) for j in range(1, 5)], 'numpy')
def mu2_at(xv, weighted=True):
    r = 'lo' if xv <= 40 else ('mid' if xv <= 60 else 'hi')
    P = np.array(dphi[r](xv), float)
    if weighted: P = P - np.array(gd(xv), float)
    return f2(*P, 0, 0, 0, xv - 2.2490776511)
if __name__ == "__main__":
    xs = np.concatenate([np.arange(10.25, 60, 0.05), np.arange(60, 262, 0.5)])
    for w in (False, True):
        v = np.array([mu2_at(a, w) for a in xs]); i = v.argmin(); print("weighted" if w else "plain W", "min mu2", v[i], "at", xs[i], "| at 100,200,260:", [round(mu2_at(a, w), 4) for a in (100, 200, 260)])
