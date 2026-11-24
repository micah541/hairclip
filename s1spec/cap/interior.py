"""
Analyticity of the Evans functions inside the contour: the Lyapunov-Perron bound at sigma0 = 32 (lintail2.eta_bound)
must hold for every lam in a neighbourhood of the closed region bounded by the contour.  By the conjugation symmetry
(all bounds are invariant under lam -> conj(lam)) it suffices to cover the upper half
    { -0.125 <= Re lam <= 1.85,  0 <= Im lam <= h(Re lam) },  h linear from 4.1 (Re = -0.125) to 0.3 (Re = 1.85),
which we cover by closed boxes (adaptive bisection), each slightly enlarged.
"""
import sys, json, time
from multiprocessing import Pool
sys.path.insert(0, '/home/claude/s1spec/lin')
import lintail2 as LT

X0, X1 = -0.125, 1.85
H0, H1 = 4.1, 0.3
PAD = 0.01


def h(x):
    return H0 + (H1 - H0) * (x - X0) / (X1 - X0)


def check(box):
    rl, rh, il, ih = box
    try:
        lam = LT.cbox(rl - PAD, rh + PAD, il - PAD, ih + PAD)
        lc = complex((rl + rh) / 2, (il + ih) / 2)
        eta, Lip, g = LT.eta_bound(0.25 - lam, 0.25 - lc)
        return box, True, float(eta.upper()), float(Lip.upper())
    except Exception as ex:
        return box, False, repr(ex), None


def main(size=0.2, nproc=2):
    t0 = time.time()
    boxes = []
    x = X0
    while x < X1 - 1e-12:
        xr = min(x + size, X1)
        top = max(h(x), h(xr))
        y = 0.0
        while y < top - 1e-12:
            boxes.append((x, xr, y, min(y + size, top)))
            y += size
        x = xr
    good, depth = [], 0
    while boxes:
        with Pool(nproc) as P:
            res = P.map(check, boxes)
        boxes = []
        for box, ok, eta, Lip in res:
            if ok:
                good.append(dict(box=box, eta=eta, Lip=Lip))
            else:
                rl, rh, il, ih = box
                if rh - rl < 0.004:
                    raise RuntimeError(f"cannot certify box {box}: {eta}")
                mr, mi = (rl + rh) / 2, (il + ih) / 2
                boxes += [(rl, mr, il, mi), (mr, rh, il, mi), (rl, mr, mi, ih), (mr, rh, mi, ih)]
        depth += 1
        print(f"level {depth}: {len(good)} boxes certified, {len(boxes)} to split ({time.time()-t0:.0f}s)", flush=True)
    json.dump(good, open("interior_result.json", "w"))
    print(f"INTERIOR OK: {len(good)} boxes, max eta = {max(g['eta'] for g in good):.4g}, max Lip = {max(g['Lip'] for g in good):.4g}")


if __name__ == "__main__":
    main()
