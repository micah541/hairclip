"""[kit file] The two aggregation steps of the slab-180 hand-off chain that the proof tree did inline (no script): writes out/h1slab/summary.json
(16 exact H1 points + 32 derivative balls) and out/h1slab/win_summary.json (32 window/shell balls).
usage: python h1slab_summaries.py h1|win [OUT_DIR]      (OUT_DIR = the p2/out directory; default: ../out relative to this file)"""
import os, sys, json, math
HERE = os.path.dirname(os.path.abspath(__file__))
O = sys.argv[2] if len(sys.argv) > 2 else os.path.join(HERE, '..', 'out')
up = lambda x: math.nextafter(x, math.inf)


def h1():
    pts = {k: json.load(open(os.path.join(O, f't9_h1_p{-1 + (2 * k + 1) / 16}_coh_s180.json'))) for k in range(16)}
    dt = {}
    for k in range(16):
        tk = -1 + (2 * k + 1) / 16
        for sg in (-1, 1): dt[(k, sg)] = json.load(open(os.path.join(O, 'h1slab', f'dt_{tk + sg * 0.03125}_0.03125.json')))
    rows = []; amax = 0.0
    for k in range(16):
        tk = -1 + (2 * k + 1) / 16
        lam = pts[k]['Lambda']['sqrtLam_total']; lip = max(dt[(k, s)]['lip_term'] for s in (-1, 1)); tot = up(lam + lip)
        rows.append([tk, lam, lip, tot, tot <= 5e-8])
        a = pts[k]['alpha']; ab = max(abs(a[0]), abs(a[1])); pd = max(abs(v) for s in (-1, 1) for v in dt[(k, s)]['pdot'][2])
        amax = max(amax, up(ab + pd / 16))
    ph = {key: [min(d['P_B'][key][0] for d in dt.values()), max(d['P_B'][key][1] for d in dt.values())] for key in ('xc', 'lam', 'alpha')}
    out = dict(rows=rows, sup_sqrtLam_base=max(r[3] for r in rows), P_hull=ph, alpha_abs_max=amax, args={'points': 16, 'balls': 32, 'slab': 180, 'u4_cell': 'exact (review B5)'})
    json.dump(out, open(os.path.join(O, 'h1slab', 'summary.json'), 'w'), indent=1)
    print('summary: sup sqrtLam base', out['sup_sqrtLam_base'], 'alpha_abs_max', amax)


def win():
    w = [json.load(open(os.path.join(O, 'h1slab', f'win_{m}.json'))) for m in range(32)]
    out = dict(windows_max=[max(x['windows_max'][j] for x in w) for j in range(4)], shells=[max(x['shells'][j] for x in w) for j in range(3)],
               balls=[f'h1slab/win_{m}.json' for m in range(32)],
               args={'balls': 32, 'radius': 0.03125, 'slab': 180, 'env_note': 'win.py per ball with T9TC, T9IFT, WINOUT (h1slab/win_cmds.txt); envelope S1 chart (no S1TT)'})
    json.dump(out, open(os.path.join(O, 'h1slab', 'win_summary.json'), 'w'), indent=1)
    print('win_summary: windows', out['windows_max'], 'shells', out['shells'])


if __name__ == '__main__':
    {'h1': h1, 'win': win}[sys.argv[1]]()
