"""Readers for T6 outputs (for T7, T8 and the P2 driver; numpy only).

  R = residuals()                       # dict from residuals.json
  norm('R2', 37)                        # certified sup over slab 37 (float, rounded up)
  norm_int('RF3_T5', 37)                # certified (int_slab N^2 dtau)^(1/2)
  interval_max('R2', 30, 34)            # max over slabs 30..34 (e.g. a T7 certificate interval)
  J = jets(37); J['k2'][e, c]           # [lo, hi] of kappa'' on sub-cell c of element e, over the whole slab
  hull(37, 'V', x_lo, x_hi)             # [lo, hi] of V~ over all sub-cells meeting [x_lo, x_hi]
Exact base polynomials (same definition as every T6 bound): base.SlabBase(i) and engine.element_a(...)."""
import os, json
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
_R = None


def residuals(path=None):
    global _R
    if _R is None or path:
        _R = json.load(open(path or os.path.join(HERE, 'residuals.json')))
        _R['_by'] = {r['slab']: r for r in _R['per_slab']}
    return _R


def norm(name, slab):
    return residuals()['_by'][slab]['norms'][name]


def norm_int(name, slab):
    return residuals()['_by'][slab]['norms_int'][name]


def interval_max(name, s0, s1, key='norms'):
    by = residuals()['_by']
    return max(by[s][key][name] for s in range(s0, s1 + 1))


def jets(slab, d=None):
    z = np.load(os.path.join(d or os.path.join(HERE, 'jets'), f'jets_{slab:03d}.npz'))
    out = {}
    for k in z.files:                       # hulls are stored round-to-nearest: widen [lo, hi] pairs outward by 2 ulps
        v = z[k]
        if not k.startswith('_') and v.dtype.kind == 'f' and v.ndim >= 1 and v.shape[-1] == 2:
            v = v.copy(); v[..., 0] = np.nextafter(np.nextafter(v[..., 0], -np.inf), -np.inf); v[..., 1] = np.nextafter(np.nextafter(v[..., 1], np.inf), np.inf)
        out[k] = v
    return out


def hull(slab, name, x_lo, x_hi, d=None):
    J = jets(slab, d)
    X = J['x']
    m = (X[..., 1] >= x_lo) & (X[..., 0] <= x_hi)
    H = J[name][m]
    return float(H[:, 0].min()), float(H[:, 1].max())
