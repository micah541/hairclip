import json, numpy as np
Ks = [10, 30, 100, 400, 1600]; R = {K: json.load(open(f'logs/f1_{K}.json')) for K in Ks}
out = []
def line(name, getter, rs=('0.25', '0.5', '1.0')):
    out.append(f"\n{name}")
    out.append("| r | " + " | ".join(f"k={K}" for K in Ks) + " | slope per ln k | DOF(k=2900) |")
    out.append("|---|" + "---|" * (len(Ks) + 2))
    for r in rs:
        v = np.array([getter(R[K]['runs'][r]) for K in Ks], float)
        b, a = np.polyfit(np.log(Ks[2:]), v[2:], 1)          # fit on k >= 100
        out.append(f"| {r} | " + " | ".join(f"{x:.0f}" for x in v) + f" | {b:.1f} | {a + b*np.log(2900):.0f} |")
out.append("sigma_end (half-curve length in rescaled units): " + ", ".join(f"k={K}: {R[K]['sigma_end']:.0f}" for K in Ks))
out.append("elements (r=0.25/0.5/1): " + "; ".join(f"k={K}: " + "/".join(str(R[K]['runs'][r]['nel']) for r in ('0.25', '0.5', '1.0')) for K in Ks))
for tol in ('1e-10', '1e-11', '1e-12', '1e-13'):
    line(f"DOF with all trailing Chebyshev coefficients of theta < {tol}", lambda e, t=tol: e['tol'][t]['dof'])
out.append("\nmax degree p (tol 1e-13): " + "; ".join(f"r={r}: " + ",".join(str(R[K]['runs'][r]['tol']['1e-13']['pmax']) for K in Ks) for r in ('0.25', '0.5', '1.0')))
out.append("max over elements of sup-error of theta'''' (truncation at tol 1e-13, bound, noise-cut 3e-15): " + "; ".join(f"r={r}: " + ",".join(f"{R[K]['runs'][r]['tol']['1e-13']['d4err_sup']:.1e}" for K in Ks) for r in ('0.25', '0.5', '1.0')))
out.append("same for tol 1e-12: " + "; ".join(f"r={r}: " + ",".join(f"{R[K]['runs'][r]['tol']['1e-12']['d4err_sup']:.1e}" for K in Ks) for r in ('0.25', '0.5', '1.0')))
out.append("sup|theta''''| itself (max over elements): " + ",".join(f"{R[K]['runs']['0.5']['d4sup_ref']:.2f}" for K in Ks))
for t4 in ('1e-08', '1e-10', '1e-12'):
    line(f"DOF so that sup error of theta'''' <= {t4} (geometric extrapolation of the coefficient tail below the 3e-15 noise floor)", lambda e, t=t4: e[f'dof_d4x_{t}'])
open('logs/f1_summary.txt', 'w').write("\n".join(out)); print("\n".join(out))
