"""Lemma S v3.2: certified contribution of the OUTER cutoff ramp [190, 379.5] lambda_1 to the strip energies (levels 0, 1, 2), python-flint arb,
outward rounding.  Inputs: ramp_out_190_379.5.json (commutator coefficient sups c_jk, arb), ../lemmaO_v32/lemmaO_v32_consts.json (Lemma O),
../lemmaO_v32/lemmaR_consts.json (Lemma R), the frozen Bbox run (sup |f| at tau_1).  Formulas: step3_proof_v32.tex, Sec. 3.7 (Lemma R paragraph).

Strip energy (deviation v = f - U_1, weight frozen at entry, cutoff chi = 1 on [50, 190], ramp on [190, 379.5]):
   Delta_j := total ramp contribution over [tau_1, oo)  <=  q sum_{k <= j+1} c_jk N_k^2,   q = 1/(1 - 4 abar),
   N_k^2 >= sup_t int_ramp z_k(v)^2 ds (k <= 2);   N_3^2 >= time average over [t_1, T) of int_ramp z_3(v)^2 ds
   (eps(tau) dtau = y_E^-4 dt / lambda_1^4, so the k = j+1 = 3 term needs only the time average).
Pieces (each inside the conclusion region |x| <= 1 of a cone member of unit L):  [190, 232.71] L = 190/3;  [232.71, 285] L = 77.57;  [285, 379.5] L = 95.
Member units, tau_m = (50/L)^2:
   kap2 = m2 + C2 sqrt(t) tau_m >= sup|h_xx|,   kap3 = n3 + sqrt(C3 t) tau_m >= ||h_xxx||_L2(-1,1),   kap4 = sqrt(2) M + sqrt(R4) >= rms_t ||h_xxxx||_L2(-1,1),
   p = tan(atan s0 + C1 sqrt(t) tau_m) >= member-frame slope,   Fb = tan(atan F_w + C1 sqrt(t) tau_m) >= |f| (ray chart), s_p = 1+p^2, s_f = 1+Fb^2.
Frames (k curvature, ' = d/dy on the ray chart, s arclength; identities verified in lemmaR_interval.check_identity):
   f' = s_f^(3/2) k,  f'' = s_f^2 (k_s + 3 f k^2),  f''' = s_f^(5/2) (k_ss + 10 f k k_s + (3 + 15 f^2) k^3),  and the same with (h, p) in the member frame.
   |k| <= kap2/L;  ||k_s||_L2(dx_m) <= Ks := kap3 + 3 p kap2^2 sqrt(len);  ||k_ss|| <= kap4 + 10 p kap2 Ks + (3 + 15 p^2) kap2^3 sqrt(len);  dy <= sqrt(s_p) dx_m.
   On a piece [a, b] (len = (b-a)/L, ell = ln(b/a)), with z_1 = y f', z_2 = y^2 f'' + y f', z_3 = y^3 f''' + 3 y^2 f'' + y f' for f and z_k(U_1) = p_U^k U_1:
   D1  = s_f^(3/2) (b/L) kap2                                         >= sup y|f'|
   N2f = (b/L)^(3/2) s_f^2 [ s_p^(1/4) Ks + 3 Fb kap2^2 sqrt(len) ]   >= ||y^2 f''||_L2(ds)
   N3f = (b/L)^(5/2) s_f^(5/2) s_p^(1/4) [ kap4 + 10 (p + Fb) kap2 Ks + (6 + 15 (p^2 + Fb^2)) kap2^3 sqrt(len) ]   >= rms_t ||y^3 f'''||_L2(ds)
   n_0 = (Fb + U_b) sqrt(ell),  n_1 = (D1 + p_U U_b) sqrt(ell),  n_2 = N2f + (D1 + p_U^2 U_b) sqrt(ell),  n_3 = N3f + 3 N2f + (D1 + p_U^3 U_b) sqrt(ell)."""
import json, math, os
from flint import arb, ctx
ctx.prec = 128
H = '/home/claude/s1spec/writeup'
A = lambda v: arb(repr(v)) if isinstance(v, float) else arb(v)
U = lambda v: math.nextafter(float(arb(v).upper()), math.inf)
J = json.load(open(f'{H}/lemmaO_v32/lemmaO_v32_consts.json'))['cases'][0]; d = J['deviation_upper_at_needed_time']
C1, C2, G3 = A(d['C1']), A(d['C2']), A(d['sqrtC3t'])
RC = {c['L_m']: c for c in json.load(open(f'{H}/lemmaO_v32/lemmaR_consts.json'))['cases']}
W = json.load(open(f'{H}/lemmaS_v32/lemmaS_Bbox_1e-05_a0.00015_dv0.15_frozen.json'))
c = json.load(open(f'{H}/lemmaS_v32/ramp_out_190_379.5.json'))['c']
s0, m2, n3, M = arb('0.08'), arb('0.03'), arb('0.004'), arb('0.12'); pU = arb('1.0246'); al = arb('1.5e-4'); abar = arb('2.19103e-7')
Fw = max(A(w['derived'][0]) for w in W)*arb('1.000001')            # |f(tau_1)| <= 0.0928 on the hand-off windows up to y = 379.7
U1 = lambda y: arb('1.3056')*al*arb(y)**pU
pieces = [(arb(190), arb('232.71'), 190/3), (arb('232.71'), arb(285), 77.57), (arb(285), arb('379.5'), 95.0)]
N = [arb(0)]*4; rows = []
for a, b, Lf in pieces:
    r = RC[Lf]; L = arb(190)/3 if Lf == 190/3 else A(Lf)
    assert abs(float(r['t_needed'])/(1.6001e-7*(50/Lf)**4) - 1) < 1e-9 and (3*L).upper() <= a.upper() + 1e-9 and b.upper() <= (5*L).upper()
    tm = (50/L)**2; kap2 = m2 + C2*tm; kap3 = n3 + G3*tm; kap4 = arb(2).sqrt()*M + A(r['R4_upper']).sqrt()
    p = (s0.atan() + C1*tm).tan(); Fb = (Fw.atan() + C1*tm).tan(); sp = 1 + p*p; sf = 1 + Fb*Fb
    ln = ((b - a)/L).sqrt(); ell = (b/a).log(); Ub = U1(b)
    Ks = kap3 + 3*p*kap2**2*ln
    D1 = sf**arb('1.5')*(b/L)*kap2
    N2f = (b/L)**arb('1.5')*sf**2*(sp**arb('0.25')*Ks + 3*Fb*kap2**2*ln)
    N3f = (b/L)**arb('2.5')*sf**arb('2.5')*sp**arb('0.25')*(kap4 + 10*(p + Fb)*kap2*Ks + (6 + 15*(p*p + Fb*Fb))*kap2**3*ln)
    n = [(Fb + Ub)*ell.sqrt(), (D1 + pU*Ub)*ell.sqrt(), N2f + (D1 + pU**2*Ub)*ell.sqrt(), N3f + 3*N2f + (D1 + pU**3*Ub)*ell.sqrt()]
    for k in range(4): N[k] = N[k] + n[k]**2
    rows.append(dict(piece=[U(a), U(b)], L=Lf, kap2=U(kap2), kap3=U(kap3), kap4=U(kap4), R4=r['R4_upper'], p=U(p), Fb=U(Fb), D1=U(D1), N2f=U(N2f), N3f=U(N3f), n=[U(x) for x in n]))
q = 1/(1 - 4*abar); entry = [arb('2e-7'), arb('1e-5'), arb('1e-5')]
Delta = [q*sum((A(c[str(j)][str(k)])*N[k] for k in range(j + 2)), arb(0)) for j in range(3)]
out = dict(ramp=[190, 379.5], pieces=rows, N_sq=[U(x) for x in N], c=c, Delta=[U(x) for x in Delta],
           Delta_over_window_entry=[U(Delta[j]/entry[j]) for j in range(3)],
           level2_split=dict(k_le_2=U(q*sum((A(c['2'][str(k)])*N[k] for k in range(3)), arb(0))/entry[2]), k_eq_3=U(q*A(c['2']['3'])*N[3]/entry[2])),
           note='Delta_j: upper bound for the total outer-ramp contribution to the level-j strip energy (omega = 28/y normalisation at tau_1). '
                'The inner edge of the strip (r = 50 lambda_1) is NOT covered here.')
json.dump(out, open(f'{H}/lemmaS_v32/ramp_cert.json', 'w'), indent=1)
for r_ in rows: print(r_)
print('N_k^2 <=', out['N_sq']); print('Delta_j <=', out['Delta'], ' / window entry (2e-7, 1e-5, 1e-5):', out['Delta_over_window_entry'], out['level2_split'])
