# FEASIBILITY PROBE for FLUX2_DESIGN.md (design stage, Nov 20 2026). arb ball arithmetic for the cutoff sups [cert, not sharp]; everything downstream is [est]. NOT part of any proof. Reads v35 outputs read-only.
import json, math
from fractions import Fraction
from flint import arb, ctx
ctx.prec = 160
def up(a): return math.nextafter(float(a.upper()), math.inf)
def amax(a, b): return a if a.upper() >= b.upper() else b
NC = 6000
def cells(lo, hi, n):
    for i in range(n):
        a = arb(lo) + (arb(hi) - arb(lo)) * i / n; b = arb(lo) + (arb(hi) - arb(lo)) * (i + 1) / n
        yield a.union(b)
def S7j(u):
    # derivatives of S7(u) = 35u^4 - 84u^5 + 70u^6 - 20u^7 w.r.t. u, orders 0..4
    return [35*u**4 - 84*u**5 + 70*u**6 - 20*u**7,
            140*u**3 - 420*u**4 + 420*u**5 - 140*u**6,
            420*u**2 - 1680*u**3 + 2100*u**4 - 840*u**5,
            840*u - 5040*u**2 + 8400*u**3 - 4200*u**4,
            840 - 10080*u + 25200*u**2 - 16800*u**3]
sup = {k: arb(0) for k in ('z1','z2','z3','z4','a3','a2','a1','a0','psi1','psi2','kapOR','b0lin')}
psi1min = None
for x in cells(57, 60, NC):
    u = (60 - x) / 3
    s = S7j(u); z = [s[i] * (arb(-1) / 3) ** i for i in range(5)]
    q = 1 + x * x; c = q.root(3)            # (1+x^2)^(1/3)
    p1 = x / (3 * c)                         # phi' = x/(3 (1+x^2)^(1/3))
    p2 = 1 / (3 * c) - 2 * x * x / (9 * c * q)   # phi''
    psi2 = p1 * p1 - p2
    A3 = 4 * z[1] - 2 * p1 * z[0]
    A2 = 6 * z[2] - 6 * p1 * z[1] + psi2 * z[0]
    A1 = 4 * z[3] - 6 * p1 * z[2] + 2 * psi2 * z[1]
    A0 = z[4] - 2 * p1 * z[3] + psi2 * z[2]
    Cq = (x - arb('2.2490')) / 4 + arb('1e-30')      # C = 7km + q/4 <= (sigma + tiny)/4, sigma <= x - 2.2490
    vals = dict(z1=z[1], z2=z[2], z3=z[3], z4=z[4], a3=A3, a2=A2, a1=A1, a0=A0, psi1=p1, psi2=psi2, kapOR=abs(z[1]) + p1 * z[0], b0lin=abs(z[4]) + Cq * abs(z[1]))
    for k, v in vals.items(): sup[k] = amax(sup[k], abs(v) if k != 'kapOR' and k != 'b0lin' else v)
    psi1min = p1 if psi1min is None or p1.lower() < psi1min.lower() else psi1min
S = {k: up(v) for k, v in sup.items()}
print('sup over OR [cert]:', json.dumps(S, indent=0))
print('exact: 35/48 =', 35/48, ' 280/27 =', 280/27, ' inf psi\' >=', float(psi1min.lower()))
# W(60)/W(57) = exp(-(phi(60)-phi(57))), phi = (1+x^2)^(2/3)/4
ph = lambda x: (1 + arb(x) ** 2) ** (arb(2) / 3) / 4
print('W(57)/W(60) <=', up((ph(60) - ph(57)).exp()))
# ---- v35 data (cert inputs) ----
B = '/home/claude/s1spec/writeup/lemmaU/bulk/'
a = json.load(open(B + 'fix2/n12/t2_A_v35/out/asm_design_N2ref4.json'))['res']['tight']
V = lambda o: float(Fraction(int(o['arb'][0])) * Fraction(2) ** o['arb'][1])
E = [V(e) for e in a['E']]; ex = a['export']
W57 = json.load(open('/home/claude/s1spec/writeup/lemmaS/lemmaS_cflux.json'))['W57_upper']
beta = 2.0; conv = 60 / 28; Eu = [beta * conv * e for e in E]
dfac = 1.0101 ** 2     # [est] d + chart Minkowski factor on N^2 (K0 with d / K0 without = 1.0201 at beta = 2, v35)
N0s = W57 * 57.76 * Eu[0] * dfac; N1s = W57 * Eu[1] / 54.75 * dfac; N2s = W57 * 2 * (Eu[2] + Eu[1]) / 54.75 ** 3 * dfac
N0, N1, N2 = [math.sqrt(v) for v in (N0s, N1s, N2s)]
print('E^A (v35, 4 N2ref) =', E, ' W57 =', W57)
print('N0,N1,N2 at beta_A=2 [est from cert]:', N0, N1, N2, '  M0,M1,M2:', N0 / math.sqrt(W57), N1 / math.sqrt(W57), N2 / math.sqrt(W57))
# level-1 export (cert v35), level-2 export GUESS
C0_1, C1_1 = V(ex['1']['C0']), V(ex['1']['C1']); r1, ac = V(ex['1']['r']), V(ex['1']['amin_c'])
ymax, ymin = 57.751, 54.7508
shell = ymax ** 5 / 28
def M3(C0, C1, boot=2.0):   # int_window M3^2 <= c + s T
    c = 3 * boot * C0 * shell / ymin ** 5 * dfac; s = 3 * (boot * C1 * shell + (ymax / 28) * beta * (9 * E[2] + 4 * E[1])) / ymin ** 5 * dfac
    return c, s
m3c, m3s = M3(C0_1, C1_1)
for C0_2, C1_2, tag in ((119.0, 0.055, 'guess r2=1/1000, D2c=E2'), (1190.0, 0.55, 'x10 pessimistic')):
    m4c = 4 * (2 * C0_2 * shell + 36 * 2 * C0_1 * shell) / ymin ** 7 * dfac
    m4s = 4 * (2 * C1_2 * shell + 36 * 2 * C1_1 * shell + (ymax / 28) * beta * (121 * E[2] + 36 * E[1])) / ymin ** 7 * dfac
    n3c, n3s, n4c, n4s = W57 * m3c, W57 * m3s, W57 * m4c, W57 * m4s
    # coefficients (alpha-level 2.5e-4, zone ZB cert values of constants_v4.json; s1 = S1V/54.75)
    cb, ca, cB = 3.3654761100478945e-3, 6.164049716517186e-05, 0.1651209793053568
    s1 = 0.6 / 54.75; C1c = ca + cB * s1
    d = dict(d010=1.5358652467512007e-4, d010a=1.5076683759218672e-06, d010b=3.7899841663602283e-3, d100=1.366435656180443e-06, d100a=6.571239316577328e-09, d100b=4.203762372867461e-06,
             d110=0.41142363427206885, d110a=3.7899841627222495e-3, d110b=10.152512157801542, d200=3.4108144063793597e-3, d200a=3.833563575028622e-06, d200b=2.5148201748379515e-3,
             d300=3.045600652694703, d300a=8.382733922189802e-4, d300b=2.2455452159047136, A1=2.4526616085318595e-09, B1=1.1424148909852663e-10)
    abar = 2.14041e-7
    # level 2: D = ||D2v||, X = ||X2||
    Dpt = S['a2'] * N2 + S['a1'] * N1 + S['a0'] * N0
    Xpt = (6 * S['z2'] * N2 + 4 * S['z3'] * N1 + S['b0lin'] * N0) + C1c * (3 * S['z1'] * N2 + 3 * S['z2'] * N1 + S['z3'] * N0) + abar * S['z1'] * 58.8 * N0 \
        + (d['d010a'] + d['d010b'] * s1 + d['d100']) * N2 + (d['d100a'] + d['d100b'] * s1) * N1 + 2 * d['d200'] * s1 * N2 + (d['d200a'] + d['d200b'] * s1) * s1 * N1 \
        + (d['d110a'] + d['d110b'] * s1) * s1 * N2 + 3 * d['d300'] * s1 ** 2 * N2 + (d['d300a'] + d['d300b'] * s1) * s1 ** 2 * N1 + d['A1'] * N1 + d['B1'] * N0
    lam1 = lam2 = 1.0
    K2pt = Dpt * Xpt + lam1 / 2 * Dpt ** 2 + lam2 / 2 * Xpt ** 2
    # integrated: Xint^2 <= 3 [ (4 z1 + d010 + d110 s1)^2 N3^2 + (d110 N2)^2 s2^2 ] (two terms -> factor 2), Dint^2 <= 2 (N4^2 + a3^2 N3^2)
    x3 = 4 * S['z1'] + d['d010'] + d['d110'] * s1
    M2 = N2 / math.sqrt(W57)
    lam = 1.0
    s2c = M2 * m3c / lam; s2s = M2 ** 2 / 3 + M2 * (lam + m3s / lam)          # int s2^2 <= s2c + s2s T
    Xi_c = 2 * (x3 ** 2 * n3c + (d['d110'] * N2) ** 2 * s2c); Xi_s = 2 * (x3 ** 2 * n3s + (d['d110'] * N2) ** 2 * s2s)
    Di_c = 2 * (n4c + S['a3'] ** 2 * n3c); Di_s = 2 * (n4s + S['a3'] ** 2 * n3s)
    phc = (1 / (2 * lam1) + 0.5) * Xi_c + (1 / (2 * lam2) + 0.5) * Di_c; phs = (1 / (2 * lam1) + 0.5) * Xi_s + (1 / (2 * lam2) + 0.5) * Di_s
    mu, eI = 0.08, 3e-5
    K0 = 3.903829755203727e-21 * 1.0065      # F1 cert (v35) x [est] quasilinear + nonlinear + modulation additions
    FL0, FL2 = 1.64e-20, 8.0e-20
    s = mu * ((FL0 - K0) + eI * (FL2 - K2pt)); Cp = mu * eI * phc; Sp = mu * eI * phs
    print('---', tag); print(' int N3^2 <= %.3e + %.3e T ; int N4^2 <= %.3e + %.3e T (bootstrap-doubled exports)' % (n3c, n3s, n4c, n4s))
    print(' Dpt %.3e Xpt %.3e  K2pt %.3e  epsI*K2pt %.3e ; K0 %.3e ; combined margin x%.3f' % (Dpt, Xpt, K2pt, eI * K2pt, K0, (FL0 + eI * FL2) / (K0 + eI * K2pt)))
    print(' slack s %.3e ; C\' %.3e (2C\' %.3e vs 8.2155e-23: x%.0f) ; S\' %.3e (s/S\' x%.2e)' % (s, Cp, 2 * Cp, 8.2155e-23 / (2 * Cp), Sp, s / Sp))
