"""usage: python3 drive.py SNAP [up] [tag]   -> res_<tag>.json (all hand-off quantities)"""
import sys, json, os, numpy as np
import prof, en, jets
from geo import *
snap = sys.argv[1]; up = int(sys.argv[2]) if len(sys.argv) > 2 else 8
tag = sys.argv[3] if len(sys.argv) > 3 else os.path.basename(snap).replace('.npz', '')
pf = f"p_{tag}.npz"
if not os.path.exists(pf):
    o, P = prof.analyse(snap, up=up, save=pf); del P
d = dict(np.load(pf))
I = en.inner(d); B = en.build(d); sh, slide, z = en.shell(B)
J = jets.run(snap, float(d['lam']), complex(d['xc']))
cr = J['cone_rows']; fr = J['far_rows']
def cone_sup(r0, r1):
    m = (cr[:, 0] >= r0 - 1e-6) & (cr[:, 0] <= r1)
    if not m.any(): return None
    return dict(zip(('slope', 'K2', 'K3', 'K4', 'N3'), [float(cr[m][:, i].max()) for i in range(2, 7)]))
res = dict(tag=tag, kmax=float(d['kmax']), kt=float(d['kt']), N=int(d['N']), lam=float(d['lam']), alpha=float(d['alpha']), alpha_k=float(d['alpha'] * d['kt']),
           xc=[float(complex(d['xc']).real), float(complex(d['xc']).imag)], fit_r2=float(d['r2']), odd_norm=float(d['oddn']), nfev=int(d['nfev']),
           chart_ext=[[float(a), float(b)] for a, b in d['chart_ext']], consist=B['cons'], ymax_data=float(B['ymax']),
           inner={k: float(v) for k, v in I.items() if not isinstance(v, np.ndarray)},
           shell={f"{a}-{b}": v for (a, b), v in sh.items()}, slide={j: [float(e), (float(ya) if ya is not None else None)] for j, (e, ya) in slide.items()},
           v_samples={str(yy): float(B['v'][np.argmin(abs(B['yg'] - yy))]) for yy in (14, 30, 50, 75, 112, 169, 253, 350)},
           U_samples={str(yy): float(B['U'][np.argmin(abs(B['yg'] - yy))]) for yy in (14, 50, 100, 253, 350)},
           f_odd_max=float(np.max(np.abs(d['f_odd']))),
           cone={'rows': cr.tolist()}, cone_sup={'200-226': cone_sup(200, 226), '200-400': cone_sup(200, 400), '200-700': cone_sup(200, 700), '200': cone_sup(200, 201)},
           far=dict(r0=float(fr[0, 0]), slope=float(fr[:, 1].max()), h2=float(fr[:, 2].max()), h3=float(fr[:, 3].max()), h4=float(fr[:, 4].max()), N3=float(fr[:, 5].max())),
           pt=dict(kappa=float(J['pt'][0]), kappa_s=float(J['pt'][1]), kappa_ss=float(J['pt'][2]), at=J['pt_r']))
json.dump(res, open(f"res_{tag}.json", 'w'), indent=1)
print(json.dumps({k: res[k] for k in ('tag', 'kmax', 'lam', 'alpha', 'alpha_k', 'inner', 'chart_ext')}, indent=0)[:1500])
