"""Review items M3 and M6 as a standalone check on the PINNED v35 outputs, so that the pinned files stay reproducible
by the scripts that produced them. (The same asserts placed inside assembly_rep.py / assemble_U.py are kept as
*_n35assert.py; a rerun with them gave identical numbers, see fix2/N35_REPORT.md.)
  M3: verifier ramp sups (with d and the chart terms) lie inside the ramp box (0.52, 0.52, 0.05)
      <=> certified lower bounds of the ramp margins 0.52/rM1, 0.52/rz1, 0.05/rM0 are > 1.
  M6: certified Nbar^2 <= 4 N2_ref at every level, so margins quoted from the 4 N2_ref run apply.
usage: python3 -B n35_asserts.py   (exit 0 = PASS)"""
import json, hashlib, os, sys
from flint import arb
B = '/home/claude/s1spec/writeup/lemmaU/bulk'
A = B + '/fix2/n12/t2_A_v35/out/A_interface.json'; U = B + '/t4_cert/out/U_cert_v35.json'
sha = lambda p: hashlib.sha256(open(p, 'rb').read()).hexdigest()
ai = json.load(open(A)); uc = json.load(open(U)); ok = True
for st in ('tight', 'v31'):
    for k in ('ramp M1', 'ramp z1', 'ramp M0'):
        m = ai['res'][st]['core']['margins'][k]; assert m['bound'] == 'lo'
        lo = arb(m['arb']).lower(); good = lo > 1; ok &= good
        print(f'M3 {st:5s} {k}: certified lower bound {float(lo):.6f} > 1: {good}')
ss = uc['supersolution']; r = ss['N2_over_4N2ref']
good = bool(ss['N2_le_4N2ref']) and max(r) < 1 - 1e-9; ok &= good
print('M6 Nbar^2/(4 N2_ref) =', [round(x, 6) for x in r], 'all < 1:', good)
json.dump(dict(pass_=ok, inputs={A: sha(A), U: sha(U), os.path.abspath(__file__): sha(__file__)}),
          open(B + '/fix2/n35_asserts.json', 'w'), indent=1)
print('PASS' if ok else 'FAIL'); sys.exit(0 if ok else 1)
