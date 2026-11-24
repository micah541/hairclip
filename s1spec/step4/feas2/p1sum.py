import json, glob, sys
fs = sorted(glob.glob('logs/p1_*_S100.json'))
print('files', fs)
for f in fs:
    R = json.load(open(f))
    for k, v in R.items():
        rows = v['rows']; info = v['info']; ok = [r['N'] for r in rows if r['margin_net'] > 0]
        c = info['const']
        s = f"{k:28s} d2={info['d2']:.3g} d1={info['d1']:.3g} d0={info['d0']:.3g} CG={c['CG']:.1f} | closes N={ok[0] if ok else 'none<=80'}"
        if ok:
            r = [x for x in rows if x['N'] == ok[0]][0]; s += f" margin={r['margin_net']:.4f} t={r['t']:.3f} nu={r['nuN1']:.2f} lamA={r['lamA']:.3f} lamR={r['lamR']:.3f} pen(d2,d1,d0)=({r['pen_d2']:.2f},{r['pen_d1']:.2f},{r['pen_d0']:.2f})"
        else:
            r = rows[-1]; s += f" (N=80: t={r['t']:.3f} lamA={r['lamA']:.3f} lamR={r['lamR']:.3f})"
        print(s)
