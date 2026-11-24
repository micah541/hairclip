import sys; sys.argv=['x']
exec(open('lint_cert_p.py').read().split("grid = [XA")[0])
for a in [30.2,30.25,34,38.0,38.2,38.24,38.245,38.25,38.3]:
    r,lo,rg=piece(a,a+0.005); print(a,rg,r.str(5), )
X=arb(38.2475,0.0025)
print([x.str(4) for x in Pj(X,'blend')]); print([x.str(4) for x in Pj(X,'rho')])
print([x.str(4) for x in [SD[j](X) for j in range(4)]])
