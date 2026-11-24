import json, numpy as np
R = json.load(open('logs/f2_results_partial.json'))
L = ["| kj | shape | amp | size_rho/d | k | A_rho=|dw|/size | A_mod=|dw|/|dw(kj)| | A_W=|dw|_W/size_W | dalpha/d | lin ratio (2d/d) |", "|---|---|---|---|---|---|---|---|---|---|"]
for kj in (10, 30, 100, 400):
    for sh in ('tip', 'chamber'):
        for a in (1, 2):
            t = f'p{kj}_{sh}_{a}'
            if t not in R: continue
            r = R[t]; d = r['delta']; info = r['info']; rows = {x['lev']: x for x in r['rows']}
            r1 = {x['lev']: x for x in R[f'p{kj}_{sh}_1']['rows']}
            for lev in (100, 400, 1600):
                if lev < kj or lev not in rows: continue
                x = rows[lev]; base = rows.get(kj)
                am = x['dw_rho'] / base['dw_rho'] if base else float('nan')
                lin = (x['dw_rho'] / d) / (r1[lev]['dw_rho'] / 1e-6) if a == 2 else float('nan')
                L.append(f"| {kj} | {sh} | {d:.0e} | {info['size_rho']/d:.2f} | {lev} | {x['dw_rho']/info['size_rho']:.2e} | {am:.2e} | {x['dw_W']/info['size_W']:.2e} | {x['dalpha']/d:+.2e} | {lin:.2f} |")
open('logs/f2_table.md', 'w').write("\n".join(L)); print("\n".join(L))
