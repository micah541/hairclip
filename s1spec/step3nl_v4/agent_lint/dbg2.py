import sys; sys.argv=['x']
exec(open('lint_cert_p.py').read().split("# x >= grid[-1]")[0].split("reg = {}")[0])
u0 = arb(68.25) ** (arb(-1) / 3); N=400
for j in list(range(0,4))+[200,399]:
    ulo, uhi = float(u0.upper()) * j / N, float(u0.upper()) * (j + 1) / N
    if j == 0: ulo = 1e-12
    xl, xh = uhi ** -3, ulo ** -3
    if xh > 1e12: xh = 1e12
    X = arb((xl + xh) / 2, (xh - xl) / 2)
    kb, mb, pb, qs = H.profile_bounds(arb(xl) + STIP.lower())
    P = Pj(X, "rho")
    print(j, xl, xh, kb.str(3), qs.str(3), [p.str(3) for p in P])
    v = MU(*P, arb(0, float(kb.upper())), arb(0, float(mb.upper())), arb(0, float(pb.upper())), X + STIP + arb(0, float(qs.upper())))
    print('   mu', v.str(5))
