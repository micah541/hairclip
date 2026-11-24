"""Linear triangular Lyapunov lemma with the slow direction phi3 = (I - Pi) psi3^ / N', N' = <(I - Pi) psi3^, f3>,
c = <u, f3>, v = u - c phi3 in {f1, f2, f3}^perp  (original certified constraints, omega = 3/64).
  c_tau = lamt c + l(u),  |l(u)| <= eps1 ||u||;   v_tau = Lcal v + c g - l(u) phi3,  g = (I - Pi) r^ / N'.
Condition for E = kappa c^2 + ||v||^2 to decay at rate 2 omega':
  (kappa eps1 + G0)^2 <= 4 kappa A0 B0,  G0 = ||g|| + eps1 ||phi3||^2,
  A0 = |lamt| - eps1 ||phi3|| - omega',  B0 = omega - eps2 - eps1 ||phi3|| - omega'.
All norms on the half line (as in the linear certificate)."""
import json
from flint import arb
L2 = json.load(open('/home/claude/s1spec/step3/lyap_consts2.json'))
IP = json.load(open('psi3ip.json')); RS = json.load(open('psi3res3.json'))
iv = lambda p: arb((p[0] + p[1]) / 2, (p[1] - p[0]) / 2)
eps1, eps2, Pin, Minv = arb(L2["eps1"]), arb(L2["eps2"]), arb(L2["Pi"]), arb(L2["Minv"])
psif3n = (arb(L2["psif3"][0]) ** 2 + arb(L2["psif3"][1]) ** 2).sqrt()
nf3 = arb(L2["normF"]) * (1 + arb("1e-4"))            # ||f3|| <= 0.0194760 + 4e-7 (table) -- generous
b1, b2, b3 = iv(IP["psi3_f1"]), iv(IP["psi3_f2"]), iv(IP["psi3_f3"])
npsi = iv(IP["psi3_norm2"]).sqrt(); rhat = arb(RS["r_half"])
bn = (abs(b1) ** 2 + abs(b2) ** 2).sqrt()
Np = b3.lower() - psif3n * Minv * bn.upper()           # N' >= <psi3^, f3> - |<Pi psi3^, f3>|
phi3 = Pin * npsi.upper() / Np
g = Pin * rhat / Np
lamabs = arb("0.2561532224577463"); OM = arb(3) / 64
G0 = g + eps1 * phi3 * phi3
def ok(wp):
    A0 = lamabs - eps1 * phi3 - wp; B0 = OM - eps2 - eps1 * phi3 - wp
    if not (A0 > 0 and B0 > 0): return False
    kap = G0 / eps1                                      # kappa = G0/eps1: condition (2 G0)^2 <= 4 kappa A0 B0
    return bool(G0 * G0 * 4 < 4 * kap * A0 * B0)
lo, hi = 0.0, float(OM)
for _ in range(50):
    mid = (lo + hi) / 2
    if ok(arb(mid)): lo = mid
    else: hi = mid
out = {"N_prime_lo": float(Np), "norm_psi3": float(npsi.upper()), "phi3": float(phi3.upper()), "g": float(g.upper()),
       "G0": float(G0.upper()), "eps1_phi3": float((eps1 * phi3).upper()), "omega_prime_max": lo}
for frac in (0.9, 0.95):
    wp = arb(lo * frac); A0 = lamabs - eps1 * phi3 - wp; B0 = OM - eps2 - eps1 * phi3 - wp; best = None
    for e in range(0, 400):
        kap = arb(10) ** (arb(e) / 40)
        if not bool((kap * eps1 + G0) ** 2 < 4 * kap * A0 * B0): continue
        Cc = ((phi3 * phi3 / kap + 1) * (kap * nf3 * nf3 + (1 + nf3 * phi3) ** 2)).sqrt()
        if best is None or Cc.upper() < best[1].upper(): best = (kap, Cc)
    out[f"omega_prime_{frac}"] = lo * frac; out[f"kappa_{frac}"] = float(best[0].mid()); out[f"C_{frac}"] = float(best[1].upper())
for k, v in out.items(): print(f"  {k:18s} {v}")
json.dump(out, open('lyap_psi3b.json', 'w'), indent=1)
