import gen_hs4  # noqa
from flint import arb
import gram_post7 as P
from gram_post7 import Aff, DSet, na, SHARED_ALL, SHARED_Z
from gen_hs3 import QSHIFT
for g in ["wa", "wb", "w"]:
    Z = DSet(f"out9Z{g}.txt", SHARED_ALL + SHARED_Z); o = 5 + na([g])
    r1 = Z.ratio(o + 7, o + 4); r3 = Z.ratio(o + 9, o + 4)
    S = DSet(f"out9Q{g}{g}.txt", SHARED_ALL); oq = 5 + na([g, g])
    main = S.form({oq + 1: r3.c, oq + 3: -r1.c, oq + 4: arb(1)}).scale(arb(2) ** QSHIFT)
    y1 = S.form({oq + 1: arb(1)}); y3 = S.form({oq + 3: arb(1)})
    d1 = Aff(arb(0), r1.d, r1.e); d3 = Aff(arb(0), r3.d, r3.e)
    t3 = (d3 * y1).scale(arb(2) ** QSHIFT); t1 = (d1 * y3).scale(arb(2) ** QSHIFT)
    f = lambda a: f"{float(a.interval().mid()):+.4e}(+-{float(a.interval().rad()):.1e})"
    print(g, "rho1.c", r1.c, "rho3.c", r3.c, "| main", f(main), "| d3*y1", f(t3), "| d1*y3", f(t1))
