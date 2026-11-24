"""For each tag: alpha*k, r2/|alpha|^3 at its snapshots (axis through the tip), predicted k* (|alpha|<=2e-3 and r2<=3e-7)."""
import sys, glob, re, numpy as np
sys.path.insert(0, '/home/claude/s1spec/step4/agent_handoff2')
from handoff2 import analyse
out = open('design_table.txt', 'a')
def P(s): print(s, flush=True); out.write(s + '\n')
for tag in sys.argv[1:]:
    files = sorted(glob.glob(f'snap_{tag}_*.npz'), key=lambda f: int(re.findall(r'_(\d+)\.npz', f)[0]))
    rows = []
    for f in files:
        try: d = analyse(f)
        except Exception as e: P(f"{tag} {f} failed {e!r}"); continue
        z = np.load(f); rows.append((d['kmax'], d['alpha'], d['alpha'] * d['kmax'], d['ratio'], d['r2'], d['r2'] / abs(d['alpha'])**3, float(z['L']), int(z['N']), d['xend']))
        P(f"{tag:24s} k={d['kmax']:7.1f} alpha={d['alpha']:+.3e} a*k={d['alpha']*d['kmax']:+7.3f} b2/a2={d['ratio']:+6.2f} r2={d['r2']:.3e} r2/|a|^3={d['r2']/abs(d['alpha'])**3:8.1f} L={float(z['L']):.3f} N={int(z['N'])} chart_end={d['xend']:.1f}")
    if rows:
        k, a, ak, b, r2, c3, L, N, xe = rows[-1]
        astar = min(2e-3, (3e-7 / c3)**(1 / 3)); kstar = abs(ak) / astar
        P(f"  -> {tag}: from last snapshot: |a*k|={abs(ak):.2f}, r2/|a|^3={c3:.1f} -> |alpha*|={astar:.2e}, k*={kstar:.0f}; L={L:.2f}; N(k*)~{N*kstar/k:.0f} (N/k={N/k:.1f})")
