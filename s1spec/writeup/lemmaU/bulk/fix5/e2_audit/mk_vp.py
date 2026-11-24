"""E2 audit: build vp/rhs_iv.py = the fin13v4 verifier's rhs_iv.py (read only) with the E2 additions switched by the module variable E2X
(E2X = None reproduces fin13v4 exactly).  Every replacement must match exactly once.  usage: python3 -B mk_vp.py"""
import hashlib
V = '/home/claude/s1spec/step3nl_v4/verifier/rhs_iv.py'
s = open(V).read(); sha = hashlib.sha256(s.encode()).hexdigest()
def rep(old, new):
    global s
    assert s.count(old) == 1, (s.count(old), old[:60]); s = s.replace(old, new)
rep("""           + ND * Dq + cl['modE2'] * AA)
    L = L + c.EPS * (fn * Dq)
""", """           + ND * Dq + cl['modE2'] * AA)
    if E2X is not None: R2p = R2p + e2_extra(al, A, W1, B, _q, sr, V0, V1, V2, V3, supw2f, wstrip, Dq, sI1r, sI2r, fn, c)      # [E2 audit]
    L = L + c.EPS * (fn * Dq)
""")
rep("""    out = [(k['f572'] * fn, 'A'), (eps * fn, 'D4'), (eps * fn * k['t22'], 'W3'), (eps * fn * k['t12'], 'B'),      # E0 forcing, E2 forcing (fin12f4 forms)
""", """    fe = (eps * fn * k['fe2']) if 'fe2' in k else (eps * fn)                                                        # [E2 audit] forcing factor of the E2 level
    out = [(k['f572'] * fn, 'A'), (fe, 'D4'), (fe * k['t22'], 'W3'), (fe * k['t12'], 'B'),      # E0 forcing, E2 forcing
""")
rep("""           (mueps * c2, 'V4'), (mueps * c2 * k['two_PSI1'], 'V3'), (mueps * c2 * k['PSI2'], 'V2')]
    return out
""", """           (mueps * c2, 'V4'), (mueps * c2 * k['two_PSI1'], 'V3'), (mueps * c2 * k['PSI2'], 'V2')]
    if 'fe2A' in k: out.append((eps * fn * k['fe2A'], 'A'))                                                         # [E2 audit] (a,b)-forcing, IBP form
    return out
""")
rep("""    kd = dict(v); kd['f572'] = cl['f572']; kd['t22'] = cl['t22']; kd['t12'] = cl['t12']
""", """    kd = dict(v); kd['f572'] = cl['f572']; kd['t22'] = cl['t22']; kd['t12'] = cl['t12']
    if E2X is not None:
        for n_ in ('fe2', 'fe2A'):
            if n_ in E2X: kd[n_] = E2X[n_]
""")
rep("""def _ei_v4(al, A, B, D4, V0, V2, V4, W1, W3, V1, V3, mid, regS, sr, sI1r, sI2r, _q, c, FL, fn):
""", '''E2X = None      # [E2 audit] dict of P constants (set by e2_probe.py); None = fin13v4
def e2_extra(al, A, W1, B, _q, sr, V0, V1, V2, V3, supw2f, wstrip, Dq, sI1r, sI2r, fn, c):
    """[E2 audit] terms missing from the E2 level, inside eps*[...] (E2_AUDIT.md section 5):
       RM*Dq (modulation beyond the Q2 lead, Cauchy-Schwarz against D_rho w~; optional IBP form ndf*RABx*A) + NS*Dq (nonlinear part on the strip)."""
    x = E2X; v = c.v4; cl = c.cl; two = v['two']
    Q = _q * _q
    S0s = sr[0] + c.sqrtPHI34 * V0
    supv2 = gm(V2, c.KW * (c.PSI1W * V2 + two * V3))
    ABs = (v['SF_ab_0'] + v['SF_ab_al'] * al) * S0s + v['SF_ab_1'] * sI1r + v['SF_ab_2'] * sI2r + v['SF_ab_3'] * supv2
    BTs = (v['SF_beta_0'] + v['SF_beta_al'] * al) * S0s + v['SF_beta_1'] * sI1r + v['SF_beta_2'] * sI2r + v['SF_beta_3'] * supv2
    ABx = (v['CAB0'] + v['CAB0_dlead'] + v['MAB'] * al) * A + v['QAB'] * Q + ABs
    RABx = (x['dab'] + v['CAB0_dlead'] + v['MAB'] * al) * A + v['QAB'] * Q + ABs          # |(a,b) - Minv l|, x-part
    BTx = (v['BETA0'] + v['BETA0_dlead'] + v['BETA1'] * al) * A + v['QB'] * Q + BTs
    ABc = v['PHIAB'] * fn
    n3e = x['n3'] + two * al * x['nh2']
    wpart = x['cw0'] * A + x['cw1'] * (W1 + x['xwc'] * S0s)                                # ||chi (A1,B1)(u) - chi (A1,B1)(u0)||
    u0part = al * (x['ngab'] + al * x['r2'])                                               # ||chi (A1,B1)(u0)||
    RM = BTx * n3e + ABx * (u0part + wpart) + ABc * wpart
    t = None
    if 'nf12' in x: RM = RM + RABx * x['nf12']
    if 'ndf' in x: t = x['ndf'] * RABx * A
    NS = x['srw'] * (cl['nf1'] * V3 + cl['nf2'] * V2 + cl['nf3'] * V1 + cl['nf7'] * (V0 + wstrip) + supw2f * (cl['nf4'] * V1 + cl['nf5'] * V2) + cl['nf6'] * (supw2f * supw2f) * V1)
    r = (RM + NS) * Dq
    return r if t is None else r + t
def _ei_v4(al, A, B, D4, V0, V2, V4, W1, W3, V1, V3, mid, regS, sr, sI1r, sI2r, _q, c, FL, fn):
''')
open('vp/rhs_iv.py', 'w').write(s)
open('vp/PROVENANCE.txt', 'w').write(f"vp/rhs_iv.py = mk_vp.py applied to {V}\nsource sha256 {sha}\npatched sha256 {hashlib.sha256(s.encode()).hexdigest()}\n")
print('ok', sha)
