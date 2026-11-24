#!/usr/bin/env python3
"""fix9_review/recompute.py -- independent spot recomputation for the review of flux lemma v2 (arb, outward rounding).
Own code for: W(57); sup_OR of zeta^(i), D2v coefficients, KAP2, level-0 Gamma (hand form with the (1,1) pair, design sec. 4);
Nbar_l from Ehat (v36); direct (non-affine) Minkowski window norms at sample T; Kbar0, K2pt, C', S', margins; C0-C6 from the rate tables.
Inputs read (read-only): v36 asm_design_N2ref4.json (E, exports), fix8 chart4_shell.json (chart/d data, ~1e-36..1e-4), c2_t3d.json sub [57,60],
rate tables fin13v4b, v36 hjw_rep.json.  Nothing is written except recompute.json next to this file.
Run: PYTHONHASHSEED=0 nice -n 19 python3 -B recompute.py"""
import json, math, os, sys, hashlib
from fractions import Fraction
from flint import arb, acb, ctx
ctx.prec = 160
B = '/home/claude/s1spec/writeup/lemmaU/bulk'
F8 = B + '/fix8_flux2_impl/out'
RTD = '/home/claude/s1spec/step3nl_v4b/runs/P'
INS = {}
def load(p):
    INS[p] = hashlib.sha256(open(p, 'rb').read()).hexdigest(); return json.load(open(p))
def up(v): return arb(math.nextafter(float(v), math.inf))
def dn(v): return arb(math.nextafter(float(v), -math.inf))
def fup(a): return math.nextafter(float(arb(a).upper()), math.inf)
def fdn(a): return math.nextafter(float(arb(a).lower()), -math.inf)
def amax(a, b): return a if a.upper() >= b.upper() else b
OUT = {}

# ---------------------------------------------------------------- A. W(57), verifier weight (S7 blend on [30.25, 38.25], W = (1+x^2)^-2 before)
def integrand(z, _):
    t = (z - acb('30.25')) / 8
    s = t**4 * (35 - 84*t + 70*t*t - 20*t**3)
    q = 1 + z*z
    return (1 - s) * 4*z/q + s * z/3 * (q.log() * (-acb(1)/3)).exp()
I = acb.integral(integrand, acb('30.25'), acb('38.25')).real
phi = lambda x: ((1 + x*x).log() * (arb(2)/3)).exp() / 4          # phi' = x/(3 (1+x^2)^(1/3))
psi57 = 2*(1 + arb('30.25')**2).log() + I + phi(arb(57)) - phi(arb('38.25'))
W57 = arb(fup((-psi57).exp()))
OUT['W57'] = fup(W57)

# ---------------------------------------------------------------- B. sups on OR
cz = load('/home/claude/s1spec/step3nl_v4b/v4cert/out/c2_t3d.json')['tables']['0.00025']['ZB']['sub']['[57,60]']
CO = {k: up(v) for k, v in cz.items() if isinstance(v, float)}
env = load(RTD + '/rate_table_fin13v4b_a0.00015.json')['args']['env']
S0V = arb(env['S0V']); S1V = arb(env['S1V'])
STIP = arb('-2.2490776515195').union(arb('-2.2490776505754'))       # v4cert/c1_cut_sups.py TIP
SMIN = arb('54.75'); s1 = arb(fup(S1V / SMIN)); ABAR = arb('1e-6')
CBAR = CO['c']; C1OR = CO['c_a'] + CO['c_b'] * s1
JB = arb(1) + arb(0, fup(CBAR)); JP = arb(0, fup(C1OR)); jbar = 1 + CBAR
ch = load(F8 + '/chart4_shell.json')['res']
X2c = load(F8 + '/flux2_coef.json')['res']
kb = up(X2c['profile']['kb']['up']); qs = up(X2c['profile']['qs']['up']); kU = CO['kU']
def S7d(u, i):
    c = [0, 0, 0, 0, 35, -84, 70, -20]
    for _ in range(i): c = [k * c[k] for k in range(1, len(c))]
    r = arb(0)
    for a in reversed(c): r = r * u + a
    return r
NC = int(os.environ.get('NCELL', 24000))
keys = ['z1', 'z2', 'z3', 'z4', 'a0', 'a1', 'a2', 'a3', 'KAP2', 'G00', 'G10', 'G20', 'G11', 'X0lin', 'psi1']
sup = {k: arb(0) for k in keys}; infp = None
for i in range(NC):
    x = (arb(57) + arb(3) * i / NC).union(arb(57) + arb(3) * (i + 1) / NC)
    u = (60 - x) / 3
    # mean-value form in u around the cell midpoint for the polynomial factors
    um = arb(u.mid()); du = u - um
    z = [S7d(um, k) + S7d(u, k + 1) * du for k in range(5)]
    z = [z[k] * (arb(-1) / 3) ** k for k in range(5)]
    z0 = z[0].intersection(arb('0.5').union(arb(0)).union(arb(1))) if z[0].overlaps(arb('0.5') + arb(0, 0.5)) else z[0]
    q13 = ((1 + x*x).log() * (-arb(1)/3)).exp()
    p1 = x / 3 * q13
    p1d = q13 / 3 - arb(2) / 9 * x * x * q13 / (1 + x*x)
    p2 = p1 * p1 - p1d
    sig = x + STIP; C = (sig + arb(0, fup(qs))) / 4 + arb(0, fup(7 * kb * up(X2c['profile']['mb']['up'])))
    Ceff = C + arb(0, fup(ABAR * ((abs(sig) + qs + 1) / (1 - kU))))
    k2 = arb(0, fup(kb * kb))
    def U(key, v):
        v = arb(fup(abs(v))) if not isinstance(v, float) else arb(v)
        if v.upper() > sup[key].upper(): sup[key] = v
    for k in range(1, 5): U('z%d' % k, z[k])
    U('a3', 4*z[1] - 2*p1*z0); U('a2', 6*z[2] - 6*p1*z[1] + p2*z0); U('a1', 4*z[3] - 6*p1*z[2] + 2*p2*z[1]); U('a0', z[4] - 2*p1*z[3] + p2*z[2])
    U('KAP2', 2*z0*z[1] - p1*z0*z0); U('psi1', p1)
    infp = p1 if infp is None or p1.lower() < infp.lower() else infp
    kapW = 4 * ((z[1]*z[1] + z0*z[2]) * JB + z0*z[1]*JP - p1*z0*z[1]*JB)          # kappa'/W, kappa = 4 W zeta zeta' j
    c20 = -kapW + 6*z0*z[2]*JB + 3*z0*z[1]*JP
    c11 = kapW / 2
    c10 = 4*z0*z[3]*JB + 3*z0*z[2]*JP + 2*k2*z0*z[1]
    c00 = z0*z[4]*JB + z0*z[3]*JP + k2*z0*z[2] + Ceff*z0*z[1]
    U('G20', c20); U('G10', c10)
    for key, c in (('G11', c11), ('G00', c00)):
        v = arb(max(fup(c), 0.0))
        if v.upper() > sup[key].upper(): sup[key] = v
    U('X0lin', jbar * abs(z[4]) + abs(Ceff * z[1]))
for k in keys: OUT['sup_' + k] = fup(sup[k])
OUT['inf_psi1'] = fdn(infp)
omW = (W57 / arb(fdn(infp))).sqrt()

# ---------------------------------------------------------------- C. pointwise data
asm = load(B + '/fix6_v36/runs/prod2/t2_A/asm_design_N2ref4.json') if os.path.exists(B + '/fix6_v36/runs/prod2/t2_A/asm_design_N2ref4.json') else load(B + '/fix6_v36/t2_A/out/asm_design_N2ref4.json')
a = asm['res']['tight']
def la(j):
    v = j['arb']
    if isinstance(v, str): return arb(v)
    m, e, r = v; x = arb(int(m)) * arb(2) ** int(e); return x + arb(0, r) if r else x
E = [arb(la(e).upper()) for e in a['E']]
EX = {k: (arb(la(a['export'][k]['C0']).upper()), arb(la(a['export'][k]['C1']).upper())) for k in ('1', '2')}
OUT['E'] = [fup(e) for e in E]; OUT['export'] = {k: [fup(v[0]), fup(v[1])] for k, v in EX.items()}
YMAX = arb('57.76') + up(ch['e']['up']); SMAX = arb('57.76')
RHO = up(ch['rho_meas']['up']); SQR = (1 + RHO).sqrt()
Ab = [up(r['up']) for r in ch['Abar']]; Bb = {k: up(r['up']) for k, r in ch['Bbar'].items()}; B422 = up(ch['B422']['up'])
Pint = [up(r['up']) for r in ch['P_int_ds']]; P2sq = up(ch['P2sq_int_ds']['up']); Dint = [up(r['up']) for r in ch['Dd_int_ds']]; Sds = up(ch['S_ds']['up'])
BETA = arb(2)
zpt = lambda k: SQR * ((YMAX / 28) * BETA * E[k]).sqrt()
y0 = zpt(0) + Dint[0].sqrt() + Ab[0] * Sds.sqrt()
y1 = zpt(1) + Dint[1].sqrt() + Ab[1] * Sds.sqrt()
r2 = zpt(2) + Dint[2].sqrt() + Ab[2] * Sds.sqrt() + Bb['22'] * Pint[2].sqrt()
M0 = SMAX.sqrt() * y0; M1 = y1 / SMIN.sqrt()
f2 = SMIN ** arb(-1.5); b22 = Bb['22'] / SMIN.sqrt()
M2 = f2 * (r2 + y1) / (1 - f2 * b22); y2 = r2 + b22 * M2
M = [arb(fup(v)) for v in (M0, M1, M2)]
N = [arb(fup((W57 * m * m).sqrt())) for m in M]
OUT['Mbar'] = [fup(m) for m in M]; OUT['Nbar'] = [fup(n) for n in N]

# ---------------------------------------------------------------- D. level 0
G = sup
K0lin_direct = G['G00'] * N[0]**2 + G['G10'] * N[0] * N[1] + G['G20'] * N[0] * N[2] + G['G11'] * N[1]**2
d = CO
cN2 = d['d010'] + d['d110'] * s1; cN1 = d['d100'] + d['d200'] * s1 + d['d300'] * s1 * s1
K0nl = (N[1] + G['KAP2'] * N[0]) * (cN2 * N[2] + cN1 * N[1] + d['A1'] * N[0] + d['A0'] * arb('0.5') * omW)
c0p = up(X2c['sup_c0p']['up'])
K0w = N[0] * arb('0.5') * omW * (c0p + ABAR * d['EA1'])
K0 = arb(fup(K0lin_direct + K0nl + K0w))
OUT['K0_direct'] = dict(total=fup(K0), lin=fup(K0lin_direct), nl=fup(K0nl), w=fup(K0w))
# the implementation's kappa form, with ITS kappa (for the comparison of the N's only)
kap = [up(r['up']) for r in X2c['kappa0']]
K0_impl_form = arb(fup(sum((kap[l] * N[l]**2 for l in range(3)), arb(0)) + K0nl + K0w))
OUT['K0_kappa_form_with_my_N'] = fup(K0_impl_form)

# ---------------------------------------------------------------- E. level 2 pointwise
z1, z2, z3, z4 = (G['z1'], G['z2'], G['z3'], G['z4'])
XN2 = jbar*6*z2 + C1OR*3*z1 + (d['d010_a'] + d['d010_b']*s1 + d['d100']) + 2*d['d200']*s1 + (d['d110_a'] + d['d110_b']*s1)*s1 + 3*d['d300']*s1*s1
XN1 = jbar*4*z3 + C1OR*3*z2 + fup(kb*kb)*2*z1 + (d['d100_a'] + d['d100_b']*s1) + (d['d200_a'] + d['d200_b']*s1)*s1 + (d['d300_a'] + d['d300_b']*s1)*s1*s1 + d['A1']
XN0 = G['X0lin'] + C1OR*z3 + fup(kb*kb)*z2 + d['B1']
Xw = d['B0'] + c0p + ABAR * d['EA1']
Dpt = G['a2']*N[2] + G['a1']*N[1] + G['a0']*N[0]
Xpt = XN2*N[2] + XN1*N[1] + XN0*N[0] + Xw*arb('0.5')*omW
K2pt = arb(fup(Dpt*Xpt + Dpt**2/2 + Xpt**2/2))
OUT['K2pt'] = fup(K2pt); OUT['X2_pt_coef'] = dict(N2=fup(XN2), N1=fup(XN1), N0=fup(XN0))
x3 = 4*z1*jbar + d['d010'] + d['d110']*s1
OUT['x3'] = fup(x3)

# ---------------------------------------------------------------- F. windows: direct Minkowski at sample T, and the implementation's affine constants
X4 = load(F8 + '/flux2_data.json')['res']
aff = {k: up(X4[k]['up']) for k in ('m3c', 'm3s', 'm4c', 'm4s', 'sig2c', 'sig2s')}
f3 = SMIN ** arb(-2.5); f4 = SMIN ** arb(-3.5)
b33 = Bb['33'] / SMIN.sqrt(); b43 = Bb['43'] / SMIN.sqrt(); b44 = Bb['44'] / SMIN.sqrt()
c3 = Dint[3].sqrt() + Ab[3]*Sds.sqrt() + Bb['32']*(M[2]/SMIN.sqrt() + Pint[2].sqrt()) + Bb['33']*Pint[3].sqrt()
c4 = Dint[4].sqrt() + Ab[4]*Sds.sqrt() + Bb['42']*(M[2]/SMIN.sqrt() + Pint[2].sqrt()) + Bb['43']*Pint[3].sqrt() + Bb['44']*Pint[4].sqrt() + 2*B422*P2sq.sqrt()
def direct(T):
    T = arb(T); rt = T.sqrt()
    zw = lambda k: SQR * (YMAX**5 / 28 * 2 * (EX[str(k-2)][0] + EX[str(k-2)][1] * T)).sqrt()
    Y1 = rt * y1; Y2 = rt * y2
    mu3 = f3 * (zw(3) + rt*c3 + 3*Y2 + 2*Y1) / (1 - f3*b33)
    Y3 = zw(3) + rt*c3 + b33*mu3
    s2int = M[2]**2 * T / 3 + 2 * M[2] * rt * mu3                      # int s2^2 <= Mbar_2^2 T/3 + 2 Mbar_2 sqrt(T) mu_3 (Cauchy-Schwarz in tau)
    R4 = zw(4) + rt*c4 + b43*mu3 + 2*B422*M[2]/SMIN.sqrt()*s2int.sqrt()
    mu4 = f4 * (R4 + 6*Y3 + 11*Y2 + 6*Y1) / (1 - f4*b44)
    return mu3**2, mu4**2, s2int
A3 = G['a3']
def Cphi(m3, m4, s2):                                                  # int_I phi2 with lam1 = lam2 = 1
    return 1 * (2*x3**2*W57*m3 + 2*(d['d110']*N[2])**2*s2) + 1 * (2*W57*m4 + 2*A3**2*W57*m3)
MU = arb('0.08'); EPSI = arb('3e-5'); FLUX0 = arb('1.64e-20'); FLUX2 = arb('8.0e-20'); F0 = MU*(FLUX0 + EPSI*FLUX2)
Cpp = arb(fup(Cphi(aff['m3c'], aff['m4c'], aff['sig2c']))); Spp = arb(fup(Cphi(aff['m3s'], aff['m4s'], aff['sig2s'])))
Cp = arb(fup(MU*EPSI*Cpp)); Sp = arb(fup(MU*EPSI*Spp))
rows = []
for T in ('0', '1e-3', '1', '100', '1e4', '1e6', '1e9'):
    m3, m4, s2 = direct(T); Td = arb(T)
    dirv = MU*EPSI*Cphi(m3, m4, s2); affv = Cp + Sp*Td
    rows.append(dict(T=T, direct_int_phi=fup(dirv), affine=fup(affv), ok=bool(dirv.upper() <= affv.upper()),
                     m3_direct=fup(m3), m3_aff=fup(aff['m3c'] + aff['m3s']*Td), m4_direct=fup(m4), m4_aff=fup(aff['m4c'] + aff['m4s']*Td),
                     m3_ok=bool(m3.upper() <= (aff['m3c'] + aff['m3s']*Td).upper()), m4_ok=bool(m4.upper() <= (aff['m4c'] + aff['m4s']*Td).upper()),
                     s2_ok=bool(s2.upper() <= (aff['sig2c'] + aff['sig2s']*Td).upper())))
OUT['windows'] = rows
OUT['Cpp'] = fup(Cpp); OUT['Spp'] = fup(Spp); OUT['Cp'] = fup(Cp); OUT['Sp'] = fup(Sp)
s_ = arb(fdn(F0 - MU*K0 - MU*EPSI*K2pt))
OUT['slack_s'] = fdn(s_); OUT['F0'] = fup(F0)
OUT['margin_combined'] = fdn((FLUX0 + EPSI*FLUX2) / (K0 + EPSI*K2pt)); OUT['margin_K0'] = fdn(FLUX0 / K0)
OUT['margin_combined_kappa_form'] = fdn((FLUX0 + EPSI*FLUX2) / (K0_impl_form + EPSI*K2pt))

# ---------------------------------------------------------------- G. conditions on the rate tables (own loop)
two = 2 * Cp; cond = {}
for al in ('0.00015', '0.000125'):
    dd = load(RTD + '/rate_table_fin13v4b_a%s.json' % al); st = dd['steps']
    Lprev = up(dd['Ly1']) * (1 + arb('1e-12')) * (1 + arb('1e-15'))
    c2 = c3m = None; nslack = 0; npos = 0; contig = True
    orow = dd['oracle_rows']
    for k, sk in enumerate(st):
        Ln = up(sk['sqrtL_next'])**2
        Lmx = amax(Lprev, Ln); Lmn = Lprev if Lprev.lower() <= Ln.lower() else Ln
        hi = arb(sk['Ly_hi']); lo = arb(sk['Ly_lo'])
        assert lo.upper() <= dn(float(Lmn.lower())).upper() or lo.upper() <= Lmn.lower(), ('barrier leaves cell below', k)
        m2 = (hi - Lmx) / two
        c2 = m2 if c2 is None or m2.lower() < c2.lower() else c2
        Uk = sk['rate_U']
        if Uk > 0:
            npos += 1
            padf = 1e-9 * max(1.0, abs(Uk))                         # as closure.py:130, float; rr - Rb may lose one ulp of Rb
            pad = arb(padf) - arb(2.3e-16) * abs(arb(Uk)) - arb(2.3e-16) * arb(padf)
            m3 = pad * Lmn.intersection(Lmn) / (arb(Uk) * two)        # barrier L(t) >= min(L_k, L_{k+1}) on the step (U>0: = L_k)
            m3b = pad * lo / (arb(Uk) * two)
            c3m = m3b if c3m is None or m3b.lower() < c3m.lower() else c3m
        Lprev = Ln
    # slack path used?  (rate_U != oracle U of the last accepted cell)
    Us = {r['U'] for r in orow}
    nslack = sum(1 for sk in st if sk['rate_U'] not in Us)
    tl = dd['tail']; ETA = arb(env['ETA'])
    rho_hi = up(abs(tl['rate_s_U'])); LF = arb(fdn(2 * F0 / rho_hi))
    ac = tl['alpha_condition_whole_tail']
    last = st[-1]
    cond[al] = dict(C2=fdn(c2), C3=fdn(c3m) if c3m is not None else None, n_steps=len(st), n_Upos=npos, n_slack_steps=nslack,
                    C4=fdn(ETA * LF / two), C5=fdn(dn(ac['lhs']) / (up(ac['rhs_whole_tail']) * (1 + two / LF))),
                    tail_in_last_cell=fdn((arb(last['Ly_hi']) - up(tl['L_star'])) / two),
                    a_star=tl['a_star'], alpha_last=last['alpha_next'], alpha_up_last=last['alpha_up'], a_star_le_alpha_up_last=bool(tl['a_star'] <= last['alpha_up']),
                    Lstar_le_last_next=bool(tl['L_star'] <= math.nextafter(last['sqrtL_next']**2, math.inf)), C1=fdn(s_ / Sp), C0=fdn(s_ / F0))
hj = load(B + '/fix6_v36/t2_A/out/hjw_rep.json')
Xc = dn(load(RTD + '/rate_table_fin13v4b_a0.00015.json')['unc']['X_const'])
res = dn(hj['kappa_J']) - up(hj['kappa_J_certified'])
cond['C6'] = dict(margin=fdn(res * Xc / Cp), margin_if_2Cp=fdn(res * Xc / two), reserve=fdn(res))
OUT['conditions'] = cond
# comparison with the implementation's certificate
cert = load(F8 + '/flux2_cert.json')['res']
cmpv = {}
for nm, mine, key in (('K0bar', K0, 'K0bar'), ('K2pt', K2pt, 'K2pt'), ('Cp', Cp, 'Cp'), ('Sp', Sp, 'Sp')):
    cmpv[nm] = dict(mine=fup(mine), impl=cert[key]['up'], ratio_mine_over_impl=float((mine / arb(cert[key]['up'])).mid()))
cmpv['margin_combined'] = dict(mine=OUT['margin_combined'], impl=cert['margin_combined']['dn'])
cmpv['W57'] = dict(mine=OUT['W57'], impl=X2c['W57']['up'])
for k_, kk in (('z1', 'z.1'), ('z2', 'z.2'), ('z3', 'z.3'), ('z4', 'z.4'), ('a0', 'frak_a.0'), ('a1', 'frak_a.1'), ('a2', 'frak_a.2'), ('a3', 'frak_a.3'), ('KAP2', 'KAP2')):
    a_, b_ = kk.split('.') if '.' in kk else (kk, None)
    iv = X2c[a_][b_]['up'] if b_ else X2c[a_]['up']
    cmpv[k_] = dict(mine=fup(sup[k_]), impl=iv, impl_ge_mine_minus_1e_3rel=bool(iv >= fup(sup[k_]) * (1 - 1e-3)))
for l in range(3): cmpv['Nbar%d' % l] = dict(mine=fup(N[l]), impl=X4['Nbar'][str(l)]['up'])
OUT['compare'] = cmpv
OUT['args'] = dict(argv=sys.argv, env={k: os.environ.get(k) for k in ('PYTHONHASHSEED', 'NCELL')}, inputs=INS, prec=ctx.prec, ncell=NC)
json.dump(OUT, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'recompute.json'), 'w'), indent=1, sort_keys=True)
print(json.dumps({k: OUT[k] for k in OUT if k not in ('args',)}, indent=1))
