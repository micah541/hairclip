"""print the T7 report tables from closure.json, t8req.json, eicon.json, c1float_clamped_0.1.json, cert_*.json"""
import json, glob, numpy as np
C1 = json.load(open('c1float_clamped_0.1.json')); E = json.load(open('eicon.json')); R = json.load(open('t8req.json')); CL = json.load(open('closure.json'))
print('| k | ' + ' | '.join(f"{r['k']:.0f}" for r in C1) + ' |')
print('| c(theta=0.03) float | ' + ' | '.join(f"{-r['top_0.03'][0]:.4f}" for r in C1) + ' |')
print('| m_c = c - omega_c (omega_c = 0.03) | ' + ' | '.join(f"{-r['top_0.03'][0]-0.03:.4f}" for r in C1) + ' |')
import closure as CLm
for nm in ('P=2', 'P=1.5', 'P=1', 'P=1|p2@60', 'P=1|rho@40'):
    r = R[nm]; b = r['best']; pk = CL[nm]['per_k']
    rho_k = [CLm.rate(p, b['muE'], 30) for p in pk]
    print(f"| {nm} | mu1 {pk[0]['mu1']:.3f}/{pk[-1]['mu1']:.3f} | a_f {r['worst']['af']:.3g} b_f {r['worst']['bf']:.3g} c<= {max(r['worst']['cI']):.2f} | "
          f"Cg0 {pk[0]['Cg0']:.3g}/{pk[-1]['Cg0']:.3g} | mu_e {b['muE']:.2g} [{r['window'][0]:.1e},{r['window'][1]:.1e}] | rho_c " + '/'.join(f'{v:.3f}' for v in rho_k) +
          f" | J2 {b['J'][2]:.2g} J4 {b['J'][4]:.2g} | A3 {b['A'][str(3)] if isinstance(b['A'], dict) and '3' in b['A'] else b['A'][3] if isinstance(b['A'], list) else b['A']['3']:.2g} | req item3 {b['item3']:.2e} K4good {b['K4good']:.2e} | Xi {b['Xi']:.1e} Yres {b['Yres']:.1e} |")
for f in sorted(glob.glob('cert_0_*.json')):
    c = json.load(open(f)); print(f"cert h={c['h']}: n={c['n']} c_cert={c['c_cert']:.6f} loss={0.06887-c['c_cert']:.5f} t_assemble={c['t_assemble']:.0f}s t_stepB={c['t_stepB']:.0f}s t_check={c.get('t_check', float('nan')):.0f}s total={c['t_total']:.0f}s")
