import json, numpy as np, sys
TAGS = sys.argv[1].split(',')           # result tags in k order
R = [json.load(open(f'res_{t}.json')) for t in TAGS]
def rows(r):
    S = r['shell']; I = r['inner']; cs = r['cone_sup']; fa = r['far']; pt = r['pt']
    out = []   # (group, name, value, threshold, kind)  pass if value <= thr
    out.append(('H1', '|alpha|  (thr 5e-4)', abs(r['alpha']), 5e-4))
    out.append(('H1', 'sqrt(Lambda)  (thr 5e-8)', I['sqrtLam'], 5e-8))
    out.append(('', '  of which sqrt(E) = (A^2+eps B^2)^.5', I['sqrtE'], None))
    out.append(('', '  of which sqrt(mu E_I)', I['sqrtmuEI'], None))
    for j, thr in ((0, 2e-7), (1, 1e-4), (2, 0.2)):
        out.append(('H2', f'E{j} on [14,60]  (thr {thr:g})', S['14-60'][j], thr))
    for w in ('50-75', '75-112', '112-169', '169-253'):
        for j, thr in ((0, 2e-7), (1, 3e-6), (2, 3e-6), (3, 3e-6), (4, 3e-6)):
            out.append(('H3', f'E{j} on [{w.replace("-", ", ")}]  (thr {thr:g})', S[w][j], thr))
    for j, thr in ((0, 2e-7), (1, 3e-6), (2, 3e-6), (3, 3e-6), (4, 3e-6)):
        out.append(('H3s', f'E{j} sliding max, windows [y,1.5y], y in [50,253]  (thr {thr:g})', r['slide'][str(j)][0], thr))
    out.append(('H4', 'chart: first y with |f|>0.5  (need >= 380)', 380.0, r['chart_ext'][0][0]))   # inverted: value=need, thr=have
    for q, thr in (('slope', 0.12), ('K2', 0.12), ('K3', 0.12), ('K4', 0.12), ('N3', 0.14)):
        out.append(('H5c', f'cone r_c=200, {q}  (thr {thr:g})', cs['200'][q], thr))
    for q, thr in (('slope', 0.12), ('K2', 0.12), ('K3', 0.12), ('K4', 0.12), ('N3', 0.14)):
        out.append(('H5c', f'cone r_c in [200,226], {q}  (thr {thr:g})', cs['200-226'][q], thr))
    for q, thr in (('slope', 0.21), ('h2', 8.9e-3), ('h3', 7.0e-5), ('h4', 6.2e-7), ('N3', 1.0)):
        out.append(('H5f', f'far cyl. (X=2,L=113, r_c>=376), {q}  (thr {thr:g})', fa[q], thr))
    for q, thr in (('kappa', 8.9e-3), ('kappa_s', 7.0e-5), ('kappa_ss', 6.2e-7)):
        out.append(('H5f', f'pointwise {q} at r>=226  (thr {thr:g})', pt[q], thr))
    return out
allr = [rows(r) for r in R]
hdr = '| quantity | ' + ' | '.join(f"k={r['kmax']:.0f}" for r in R) + ' |\n|---|' + '---|' * len(R)
print(hdr)
for i in range(len(allr[0])):
    g, name, _, thr = allr[0][i]
    cells = []
    for rr in allr:
        _, _, v, t = rr[i]
        if t is None: cells.append(f'{v:.2e}'); continue
        if g == 'H4':
            cells.append(f'{t:.0f} ' + ('OK' if t >= v else 'FAIL')); continue
        m = t / v if v > 0 else float('inf')
        cells.append(f'{v:.2e} ' + (f'OK x{m:.1f}' if v <= t else f'FAIL /{v/t:.2f}'))
    print(f'| {name} | ' + ' | '.join(cells) + ' |')
