import json, numpy as np, glob
L = lambda f: json.load(open(f))
A = L('logs/f3b.json'); A.update({k: v for k, v in L('logs/f3b_extra.json').items()}); F = L('logs/f3bF.json')
ts = sorted([k for k in A if k.startswith('t=')], key=lambda k: float(k[2:]))
print('### F3b main (variant A)')
print('| t | lam | k_tip | axis crossings (sigma) | top 6 eig sym_rho | n(>-0.05) | n(>-0.1) | n(>-0.25) | M(0.1) | variants n(>-0.05): B (no tau) / E (consistent mu,wc) / C (mu=0) | F (T-clock fixed centre) |')
print('|---|---|---|---|---|---|---|---|---|---|---|')
for k in ts:
    a = A[k]; x = a['A']; fk = F.get(k)
    print(f"| {float(k[2:]):g} | {a['lam']:.3f} | {a['kt']:.3f} | {a['axis_crossing_sigma'][:4]} | {x['top'][:6]} | {x['n']['0.05']} | {x['n']['0.1']} | {x['n']['0.25']} | {x['M']['0.1']} | {a['B']['n']['0.05']} / {a['E']['n']['0.05']} / {a['C']['n']['0.05']} | {fk['n']['0.05'] if fk else '-'} (n.1={fk['n']['0.1'] if fk else '-'}) |")
print('mu, wc per t:', {k: (round(A[k]['mu'], 4), round(A[k]['wc'][0], 3)) for k in ts})
print('F details:', {k: (F[k]['lam'], F[k]['tipY'][0], F[k]['top'][:4]) for k in F})
print('modes (A) 4th/5th eigenvalue, location s50/s90:')
for k in ts:
    m = A[k]['A']['modes']; print(k, [(x['eig'], x['s50'], x['s90']) for x in m[3:5]])
