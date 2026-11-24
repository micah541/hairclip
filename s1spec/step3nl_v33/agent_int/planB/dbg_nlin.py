import sys; sys.argv=['x','1e-3','1e-2','1e-2','1']
src=open('pb_nlin.py').read().split("if __name__")[0]
exec(src)
bl = make_segments(BPS, [w * 1 for w in (0.02, 0.02, 0.025, 0.05, 0.1, 0.2, 0.2, 0.2)])
for n,(sa,sb) in enumerate(bl):
    b = make_box(sa, sb)
    y = PSI.jets_any(b.X, b.xc, 5); h = H2.jets_any(b.X, b.xc, 5)
    try:
        d = box_vals(y, h, b.p, b.q, b.k, b.m)
    except AssertionError as e:
        print('assert', n, b.xa, e); break
    bad=[k_ for k_,v in d.items() if v.is_nan()]
    if bad: print(n,b.xa,bad,b.k,b.m,b.p,b.q,[str(t) for t in y[:3]]); break
