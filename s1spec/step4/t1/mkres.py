"""residual table for the report from polish_P11.log (full polish from the seed) and resfeas outputs (chains from the float base)"""
import json, re
rows = []
for l in open('polish_P11.log'):
    if l.startswith('NORMS'): d = json.loads(l[6:]); d['src'] = 'seed'; rows.append(d)
seed0 = {}
for f in ('resfeas_P11_a.out',):
    for l in open(f):
        m = re.match(r'slab\s+(\d+) k=\s*([\d.]+) n=(\d+): (.*?)\s+\(polish', l)
        if m and int(m.group(1)) in (0, 1): d = dict(slab=int(m.group(1)), k=float(m.group(2)), n=int(m.group(3)), src='seed')
        elif m and int(m.group(1)) == 59: d = dict(slab=59, k=float(m.group(2)), n=int(m.group(3)), src='float+5')
        else: continue
        for kv in m.group(4).split(): k, v = kv.split('='); d[k] = float(v)
        rows.append(d)
for l in open('resfeas_P11_b.out'):
    m = re.match(r'slab\s+(\d+) k=\s*([\d.]+) n=(\d+): (.*?)\s+\(polish', l)
    if m and int(m.group(1)) in (122, 184):
        d = dict(slab=int(m.group(1)), k=float(m.group(2)), n=int(m.group(3)), src='float+4')
        for kv in m.group(4).split(): k, v = kv.split('='); d[k] = float(v)
        rows.append(d)
rows = [d for d in rows if not (d['src'] == 'seed' and d['slab'] in (90, 150))]
rows.sort(key=lambda d: (d['src'] != 'seed', d['k']))
K = ['R2', 'R1', 'Rs', 'RIs', 'RF0', 'RF1', 'RF2', 'RF3', 'RF4']
T = dict(R2=2e-14, R1=5e-13, Rs=5e-11, RIs=1e-10, RF0=1e-7, RF1=1e-7, RF2=1e-7, RF3=1e-7, RF4=1e-5)
out = ['| slab (k) | start | n | R2 | R1 | Rs | Rᴵs | RF0 | RF1 | RF2 | RF3 | RF4 |', '|---|---|---|' + '---|' * 9,
       '| **target** | | | 2e-14 | 5e-13 | 5e-11 | 1e-10 | 1e-7 | 1e-7 | 1e-7 | 1e-7 | 1e-5 |']
for d in rows:
    cells = [(f'**{d[k]:.1e}**' if d[k] > T[k] else f'{d[k]:.1e}') for k in K]
    out.append(f"| {d['slab']} ({d['k']:.0f}) | {d['src']} | {d['n']} | " + ' | '.join(cells) + ' |')
print('\n'.join(out))
