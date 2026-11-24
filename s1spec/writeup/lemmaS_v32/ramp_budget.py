"""Lemma S v3.2, open item (Lemma R): budget for the outer cutoff ramp of the strip (B) energies.  FLOAT ESTIMATE, not a certificate
(the coefficients c_jk are arb sups from scratch_lemmaS/ramp_out.py; the ramp norms below use a blanket frame factor FR = 1.1).
Level-j energy (frozen weight, deviation v = f - U_1):   Delta_j := total contribution of the ramp over [tau_1, oo)
     <= (1 - 4 abar)^-1 sum_{k <= j+1} c_jk * N_k^2,      N_k^2 = sup_t int_ramp z_k^2 ds   (k = j+1: time average over [t_1, T], since
     eps(tau) dtau = dt / lambda_1^4 y_entry^4).
Ramp norms from the Lemma O conclusions on the ramp (member r_m = 4r/3, L_m = r/3 for r <= 285, else L_m = 95; deviations at the
horizon times (50/L_m)^2) plus the hand-off deviation (window sups / energies of the frozen Bbox run):
  |z_0| <= S0 + 1.01 C1 sqrt t,   |z_1| <= S1 + FR (r/L_m) C2 sqrt t,   ||z_2||_L2(ds) <= sum sqrt(Eu_2) + ||z_1|| + FR (sum_pieces (r/L_m)^3 C3 t)^(1/2).
k = 3 (level 2 only) is NOT available from Lemma O as stated.  Two entries are printed for it:
  B3      = the largest time average of int_ramp z_3^2 ds for which Delta_2 <= X * (entry E_2 = 1e-5), X = 1 and X = 9;
  est_O   = the value suggested by Lemma O's own energy E_3 (dissipation int int g_5^2 <= C3 t / (2 a_min), interpolation
            int int g_4^2 <= (int int g_3^2 int int g_5^2)^(1/2)): time average of int_{|x|<=1} g_4^2 <= C3 sqrt(t) / sqrt(2 a_min) in member units."""
import json, math
H = '/home/claude/s1spec/writeup'
C = json.load(open(f'{H}/lemmaO_v32/lemmaO_v32_consts.json'))['cases'][0]; d = C['deviation_upper_at_needed_time']
C1, C2, G3, C3 = d['C1'], d['C2'], d['sqrtC3t'], C['C3_upper']; tneed = C['needed']
W = json.load(open(f'{H}/lemmaS_v32/lemmaS_Bbox_1e-05_a0.00015_dv0.15_frozen.json'))
abar = 2.19103e-7; FR = 1.1; amin = (16/17)**2; entry = [2e-7, 1e-5, 1e-5]; out = {}
def Lm(r): return r/3 if r <= 285 else 95.0
def pieces(y0, y1):                      # conclusion regions [3 r_m/4, 5 r_m/4] covering the ramp, largest unit first
    P = []; x = y0
    while x < y1 - 1e-9:
        rm = min(4*x/3, 380.0); hi = min(5*rm/4, y1); P.append((x, hi, rm/4)); x = hi
    return P
for y0, y1 in ((150., 175.), (150., 300.), (190., 380.)):
    c = json.load(open(f'{H}/lemmaS_v32/ramp_out_{y0:g}_{y1:g}.json'))['c']; ell = math.log(y1/y0)
    ws = [w for w in W if w['yE'] < y1 and 1.5*w['yE'] > y0]
    S0 = max(w['S'][0] for w in ws); S1 = max(w['S'][1] for w in ws); sE2 = sum(math.sqrt(1.5*w['yE']/28*w['EB'][2]) for w in ws)
    tm = lambda r: (50/Lm(r))**2
    d0 = S0 + 1.01*C1*tm(y0); d1 = S1 + FR*max((r/Lm(r))*C2*tm(r) for r in (y0, min(y1, 285.), y1))
    P = pieces(y0, y1)
    g2 = FR*math.sqrt(sum((hi/L)**3*(G3*(50/L)**2)**2 for lo, hi, L in P))
    N = [d0*d0*ell, d1*d1*ell, (sE2 + d1*math.sqrt(ell) + g2)**2]
    q = 1/(1 - 4*abar); res = {}
    for j in (0, 1, 2):
        part = q*sum(c[str(j)][str(k)]*N[k] for k in range(min(j + 1, 2) + 1))
        res[f'Delta_{j}_over_entry_k_le_2'] = part/entry[j]
    part2 = res['Delta_2_over_entry_k_le_2']*entry[2]; c23 = q*c['2']['3']
    res['B3_for_Delta2_le_1x_entry'] = (1e-5 - part2)/c23; res['B3_for_Delta2_le_9x_entry'] = (9e-5 - part2)/c23
    estO = sum((hi/L)**5*C3*math.sqrt(tneed)*(50/L)**2/math.sqrt(2*amin) for lo, hi, L in P)*FR**2
    res['est_O_time_avg_int_z3sq'] = estO; res['Delta_2_over_entry_with_est_O'] = (part2 + c23*estO)/entry[2]
    s5 = sum((hi/L)**5 for lo, hi, L in P)*FR**2
    res['K4bar_sq_budget_member_units_1x'] = res['B3_for_Delta2_le_1x_entry']/s5; res['K4bar_sq_budget_member_units_9x'] = res['B3_for_Delta2_le_9x_entry']/s5
    res.update(ell=ell, N_sq=N, c=c, pieces=P, eps_edge=y0**-4)
    out[f'[{y0:g},{y1:g}]'] = res
    print(f'ramp [{y0:g},{y1:g}]: Delta_j/entry (k<=2) = ' + ', '.join('%.3g' % res[f'Delta_{j}_over_entry_k_le_2'] for j in (0, 1, 2)) +
          f" | level 2, k = 3: budget B3 = {res['B3_for_Delta2_le_1x_entry']:.3g} (1x), {res['B3_for_Delta2_le_9x_entry']:.3g} (9x); est_O = {estO:.3g} -> Delta_2/entry = {res['Delta_2_over_entry_with_est_O']:.3g}; needed time-avg int_|x|<=1 h_xxxx^2 (member units) <= {res['K4bar_sq_budget_member_units_1x']:.3g} (1x), {res['K4bar_sq_budget_member_units_9x']:.3g} (9x)")
json.dump(out, open(f'{H}/lemmaS_v32/ramp_budget.json', 'w'), indent=1)
