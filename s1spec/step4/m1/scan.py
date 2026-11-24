import numpy as np, json, sys, en
THR_W = [2e-7, 3e-6, 3e-6, 3e-6, 3e-6]
def first_violation(tag):
    d = dict(np.load(f'p_{tag}.npz')); B = en.build(d); out, slide, z = en.shell(B)
    yg, sg = B['yg'], B['sg']; om = 28.0 / yg
    ystar = []
    yas = np.exp(np.arange(np.log(14), np.log(253) + 1e-9, 0.01))
    E = np.full((5, len(yas)), np.nan)
    for i, ya in enumerate(yas):
        yb = 1.5 * ya
        sel = (yg >= ya) & (yg <= yb) & np.isfinite(z[4])
        if sel.sum() < 10 or yg[sel][-1] < 0.995 * yb or yg[sel][0] > 1.005 * ya: continue
        for j in range(5): E[j, i] = en.wint(om * np.nan_to_num(z[j])**2, sg, np.log(ya), np.log(yb))
    for j in range(5):
        bad = np.where(E[j] > THR_W[j])[0]
        # first start y >= 50 with violation
        bad = bad[yas[bad] >= 50 - 1e-9]
        ystar.append(float(yas[bad[0]]) if len(bad) else None)
    return ystar, d
if __name__ == '__main__':
    for tag in sys.argv[1:]:
        ys, d = first_violation(tag); print(tag, "first window start y violating E0..E4:", ys)
