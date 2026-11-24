"""F4: absolute CFLUX constants from F1 energies; N3 via time-integrated level-1 dissipation.
K0_rate, K2_rate: per unit tau (N0..N2 parts);  I3 = kappa3 * int N3^2 dtau over horizon T (time-integrated; NOT a rate).
Dissipation: dE1/dtau <= -2 D1 + sources,  D1 = int A2 eps omega z_3^2 ds (core).  Hence
  int 2 D1 <= E1(tau1) + int(inflow+ramp densities) + 2 c_F sqrt(E1max) int A dtau + eps_core*Ghat11*E1max*T   [assumes F6: int sqrtLam <= 2 sqrtLam1/rho]
On the window y in [54.75,57.75]: z_3^2 integrated ds <= D1/(amin*eps(57.75)*omega(57.75)) (omega=(y/28)^-1).
"""
import json, math
d = json.load(open('lemmaS_F1.json')); c = json.load(open('lemmaS_cflux.json'))
W = c['W57_upper']; k0 = c['kappa_I0']; k2 = c['kappa_I2']
rho = 2.6e-3; sqL = 2.49e-7; intA = 2*sqL/rho
al = 5e-4; p = 1.0246; Rout = 350.; y0 = 14.; gam = 1.0
e = 2*p - gam; Iint = 28.**gam*((Rout**e - y0**e)/e)
out = {}
for T in (1e3, 1e4):
  for nm in ('lin_gam1.0', 'quasi_box1_gam1.0', 'quasi_box10_gam1.0'):
    r = d[nm]['result']; Eu = [float(x) for x in r['Eu']]; E1 = float(r['E'][1]); pa = r['parts'][1]
    cF = 1.4*(0.00244 + 63.5*al + 1.02*al*7.57)*p*math.sqrt(Iint)       # F = cF * A
    amin = (1 + 0.1**2)**-2; eps_c = 28.**-4
    G11 = float(d[nm]['consts']['Ghat_1']['1'])
    bud = float(pa['handoff']) + float(pa['inflow']) + float(pa['ramp_comm']) + 2*cF*math.sqrt(E1)*intA + eps_c*G11*E1*T
    yw = 57.75; denom = amin*yw**-4*(28./yw)
    Eu3_int = (bud/2)/denom * (yw/28.)**gam       # int dtau of unweighted int z_3^2 ds over window  (omega^-1 conv)
    # N_j^2 as in F4 plan
    N0 = W*57.75*Eu[0]; N1 = W*Eu[1]/54.75; N2 = W*2*(Eu[2] + Eu[1])/54.75**3
    N3int = W*3*(Eu3_int + T*(9*Eu[2] + 4*Eu[1]))/54.75**5      # int N3^2 dtau
    K0 = k0[0]*N0 + k0[1]*N1 + k0[2]*N2
    K2r = k2[0]*N0 + k2[1]*N1 + k2[2]*N2
    I3 = k2[3]*N3int
    out[f'{nm}_T{T:g}'] = dict(K0_rate=K0, K2_rate=K2r, I3=I3, K2_total_over_T=K2r*T + I3, budget=bud, Eu3_int=Eu3_int)
    print(f'T={T:g} {nm}: K0={K0:.2e}/tau  K2={K2r:.2e}/tau  I3={I3:.2e}  (K2*T+I3={K2r*T+I3:.2e})  budget={bud:.2e}')
json.dump(out, open('lemmaS_F4.json', 'w'), indent=1)
