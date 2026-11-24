"""Final assembly: from run_<alpha>_<eta>_n<NC>.json + global constants (ibp_consts) to the final inequalities.   usage: python3 ibp_final.py RUN.json [OUT.json]"""
import sys, json
import numpy as np
from flint import arb, ctx
ctx.prec = 150
import ibp_consts as K
from simplex import simplex_max
R = json.load(open(sys.argv[1])); OUT = sys.argv[2] if len(sys.argv) > 2 else None
AL = R['alpha1']; ETA = R['eta']
WMAX = 1e-3; AMAX = 1e-3; SREF = 0.02            # hypotheses of the feedback (Neumann) step: A, W1 <= 1e-3, |w|,|w'| <= 0.02 on the strip
up = lambda z: float(z.upper()); mag = lambda z: arb(abs(z).upper())
comps = ['a', 'b', 'beta']; CN = ['00', '01', '11', '02', '12', '22']
def forms_of(reg, i):
    g = R[reg][comps[i]]
    Q = {c: g['gam'][c] for c in CN}
    C = {'00': g['M3'] / 3, '01': 2 * g['M3'] / 3, '11': g['M3'] / 3, '02': g['M2'], '12': g['M2'], '22': 2 * g['M1']}      # cubic remainder / eta
    return Q, C
def build(reg):
    out = {}
    for i in range(3):
        Q, C = forms_of(reg, i)
        Ce = {c: ETA * v for c, v in C.items()}
        J = {c: Q[c] + Ce[c] for c in CN}
        out[comps[i]] = dict(Q=Q, Ce=Ce, J=J, q_quad=simplex_max([Q])[0], q_cubic=simplex_max([Ce])[0], q_joint=simplex_max([J])[0],
                             q_simple_sum=simplex_max([Q])[0] + simplex_max([Ce])[0])
    Qa, Qb = out['a']['Q'], out['b']['Q']; Ca, Cb = out['a']['Ce'], out['b']['Ce']; Ja, Jb = out['a']['J'], out['b']['J']
    out['ab'] = dict(q_quad=simplex_max([Qa, Qb])[0], q_cubic=simplex_max([Ca, Cb])[0], q_joint=simplex_max([Ja, Jb])[0])
    return out
res = {}
for reg in ('core', 'near'):
    res[reg] = build(reg)
# ------------------------------------------------------------------ linear part
m_i = [R['m_i'][c] for c in comps]; off = R['offdiag_ab']
lam = (m_i[0]**2 + m_i[1]**2 + np.sqrt((m_i[0]**2 - m_i[1]**2)**2 + 4 * off**2)) / 2
m_ab = float(np.sqrt(lam)) * (1 + 1e-12)
m_j = R['m_j']; c_j = [up(x) for x in K.c_j]
lead = [up(x) for x in K.lead]; lead_ab = up(K.lead_vec); cF = [up(x) for x in K.cF]; cF_ab = up(K.cF_vec)
# ------------------------------------------------------------------ Neumann feedback (E matrix), hypotheses |alpha| <= AL, W1 <= WMAX, A <= AMAX, strip |w|,|w'| <= SREF
Nref = [[K.N0[i][k] * AL + K.N1[i][k] * WMAX + K.N2[i][k] * AL * AMAX + K.N3[i][k] * SREF + K.N4[i][k] * AL * SREF for k in range(3)] for i in range(3)]
Nf = np.array([[up(Nref[i][k]) for k in range(3)] for i in range(3)])
ev, V = np.linalg.eig(Nf); v = np.abs(V[:, np.argmax(ev.real)].real); d = v + 1e-3 * v.max(); d = d / d.max()
dA = [arb(float(x)) for x in d]
nu = max(up(sum((Nref[i][k] * dA[k] for k in range(3)), arb(0))) / float(d[i]) for i in range(3)) * (1 + 1e-12); assert nu < 1
nd = lambda Nm: [sum((Nm[i][k] * dA[k] for k in range(3)), arb(0)) for i in range(3)]
n0, n1, n2, n3, n4 = ([up(x) for x in nd(Nm)] for Nm in (K.N0, K.N1, K.N2, K.N3, K.N4))
fac = 1 / (1 - nu)
def finalize(reg, qkey='q_joint'):
    qi = [res[reg][c][qkey] for c in comps]
    tauA0 = max(lead[k] / d[k] for k in range(3)); tauA1 = max(m_i[k] / d[k] for k in range(3))
    tauQ = max(qi[k] / d[k] for k in range(3)); tauP = max(cF[k] / d[k] for k in range(3))
    o = {}; formsT = {}
    for i, c in enumerate(comps):
        dm = n0[i] * tauA0 * fac + n0[i] * tauA1 * AL * fac                       # alpha*A part from feedback
        tA = tauA0 + tauA1 * AL
        fb_AW1 = n1[i] * tA * fac                                                    # coefficient of W1*A   (E-matrix coupling)
        fb_AA = n2[i] * AL * tA * fac                                                # alpha*A*A
        sQ = (n0[i] * AL + n1[i] * WMAX + n2[i] * AL * AMAX) * fac * tauQ            # (alpha, W1, alpha A) * quadratic form  (scalar multiple of Q)
        J = dict(res[reg][c]['J'])
        Jfb = dict(J); Jfb['01'] += fb_AW1; Jfb['00'] += fb_AA
        for key, mult in (('00', 1), ('01', 2), ('11', 1), ('02', 2), ('12', 2), ('22', 1)): Jfb[key] += mult * sQ
        formsT[c] = Jfb
        q_tot = simplex_max([Jfb])[0]
        dlead = (n3[i] + AL * n4[i]) * SREF * tauA0 * fac                            # strip-sup coupling (tiny)
        gam = (AL * n0[i] + WMAX * n1[i] + AL * AMAX * n2[i] + n3[i] * SREF + AL * n4[i] * SREF) * fac
        o[c] = dict(lead=lead[i], dlead=dlead, m_new=m_i[i], dm_feedback=dm, m_total=m_i[i] + dm, q_ibp_quad=res[reg][c]['q_quad'], q_ibp_cubic=res[reg][c]['q_cubic'],
                    q_ibp_joint=qi[i], q_feedback_forms=dict(W1A=fb_AW1, AA=fb_AA, scalar_Q=sQ), q_total=q_tot,
                    cPhi=cF[i], cPhi_total=cF[i] + gam * tauP, gamma_ref=gam, d=float(d[i]))
    qab_raw = res[reg]['ab'][qkey]; qab_tot = simplex_max([formsT['a'], formsT['b']])[0]
    o['ab'] = dict(lead=lead_ab, dlead=float(np.hypot(o['a']['dlead'], o['b']['dlead'])), m_new=m_ab, dm_feedback=float(np.hypot(o['a']['dm_feedback'], o['b']['dm_feedback'])),
                   m_total=m_ab + float(np.hypot(o['a']['dm_feedback'], o['b']['dm_feedback'])), q_ibp_quad=res[reg]['ab']['q_quad'], q_ibp_cubic=res[reg]['ab']['q_cubic'], q_ibp_joint=qab_raw,
                   q_total=qab_tot, cPhi=cF_ab, cPhi_total=float(np.hypot(o['a']['cPhi_total'], o['b']['cPhi_total'])))
    o['taus'] = dict(tauA0=tauA0, tauA1=tauA1, tauQ=tauQ, tauPhi=tauP)
    return o
final = {reg: finalize(reg) for reg in ('core', 'near')}
# ------------------------------------------------------------------ strip terms (all multiplied by tiny rho), in terms of S_i = sup_{strip}|w^(i)|
strip = {}
for i, c in enumerate(comps):
    g = R['strip_L1'][c]
    lin_old = [sum(float(abs(K.H[i][j].mid())) * K.KAP[j][ip] * 1.0 for j in range(3)) for ip in range(4)]
    strip[c] = dict(linear_alpha_S0=R['L1k_strip_i'][c], commutator_S0_S3=[up(sum((mag(K.H[i][j]) * arb(K.KAP[j][ip]) for j in range(3)), arb(0))) for ip in range(4)],
                    quad_gamma_L1=g['gam'], cubic_M3_L1=g['M3'], cubic_M2_L1=g['M2'], cubic_M1_L1=g['M1'])
out = dict(alpha1=AL, eta=ETA, NC=R.get('nboxes'), hypotheses=dict(alpha_max=AL, A_max=AMAX, W1_max=WMAX, strip_sup_w_wprime=SREF, jet_box_core=ETA),
           per_pairing=dict(c_j=c_j, m_j=m_j, note="|<w~,L*f_j>| <= c_j A (projected, j=1,2; eps1 for j=3); |alpha <w,(d_alpha L)* chi f_j>| <= m_j |alpha| A"),
           neumann=dict(nu=nu, d=[float(x) for x in d], n0=n0, n1=n1, n2=n2, n3=n3, n4=n4, one_over_1_minus_nu=float(fac)),
           quad_cubic=res, final=final, strip=strip, lead=dict(zip(comps, lead), ab=lead_ab), m_new_Gram=dict(zip(comps, m_i), ab=m_ab))
if OUT: json.dump(out, open(OUT, 'w'), indent=1)
def show(reg):
    f = final[reg]
    print(f"--- region {reg}: nu = {nu:.4f}")
    for c in ('ab', 'a', 'b', 'beta'):
        x = f[c]; print(f"{c:5s} lead {x['lead']:.5g} m_new {x['m_new']:.4g} m_fb {x['dm_feedback']:.4g} m_tot {x['m_total']:.4g} | q_ibp(quad/cub/joint) {x['q_ibp_quad']:.4g}/{x['q_ibp_cubic']:.4g}/{x['q_ibp_joint']:.4g} q_tot {x['q_total']:.4g} | cPhi {x['cPhi_total']:.4g}")
    print('   feedback forms:', {c: {k_: round(v_, 4) for k_, v_ in f[c]['q_feedback_forms'].items()} for c in comps})
for reg in ('core', 'near'): show(reg)
print("per-pairing c_j", c_j, "m_j", m_j)
