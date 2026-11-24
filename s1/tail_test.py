from tail import *
F=lambda x: float(x.mid())
for S,em in [(6,"2.0e-3"),(8,"1.7e-4"),(10,"1.1e-5")]:
    r=radii(S, arb(em))
    print(f"S={S}: kappa={F(r['kappa']):.5g} r0={F(r['r0']):.4g} r1={F(r['r1']):.4g} delta1={F(r['delta1']):.4g} (rel {F(r['delta1'])/float(em):.3g}) L={F(r['L']):.4f} Jbar={F(r['Jbar']):.3g} verified={r['verified']}")
    g=gamma1(S)
    print("   gamma1 =", F(g.real), F(g.imag), " rad", float(g.real.rad()), float(g.imag.rad()))
