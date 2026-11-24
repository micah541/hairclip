import json, math, sys
sys.path.insert(0,'.')
raw = json.load(open('cflux_raw.json'))
y0 = 57 - 2.2490776515
W57 = raw['W57_upper']; IW = 3*W57
box = raw['box']; Uj = 0.11
def kappas(Gm, nmax):
    k = [0.0]*(nmax+1)
    for key, v in Gm.items():
        a, b = map(int, key.split(','))
        if a == b: k[a] += v
        else: k[a] += v/2; k[b] += v/2
    return k
k0 = kappas(raw['Gamma']['K0'], 3); k2 = kappas(raw['Gamma']['K2'], 3)
def Nsq(j, sharp):
    if sharp and j == 0: m = 0.0129
    elif sharp and j == 1: m = 0.219
    else: m = box['M%d' % j] + Uj
    return IW*(m/y0**j if not (sharp and j == 1) else m/y0)**2 if True else None
out = {'ramp': [57, 60], 'zeta': 'quintic smoothstep 1->0 (assumed)', 'W57_upper': W57, 'intW_ramp_upper': IW,
       'Gamma': raw['Gamma'], 'kappa_I0': k0[:3], 'kappa_I2': k2[:4],
       'CFLUX_in_uncut_ramp_norms': {'I0': max(k0[:3]), 'I2': max(k2[:4])}}
for tag, sharp in (('box', False), ('lemmaS_sharpened', True)):
    a0 = sum(k0[j]*Nsq(j, sharp) for j in range(3)); a2 = sum(k2[j]*Nsq(j, sharp) for j in range(4))
    out['abs_%s' % tag] = {'K0_abs': a0, 'K2_abs': a2, 'K2_abs_times_epsI_over_epsI': a2}
out['outflow_coeff_c_zeta_zetap_max'] = (60 - 2.249)/4*0.398
json.dump(out, open('lemmaS_cflux.json', 'w'), indent=1)
print(json.dumps(out, indent=1))
