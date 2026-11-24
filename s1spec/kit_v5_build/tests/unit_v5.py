"""Light unit tests of the v5 Windows-fix code paths that a WSL run does not exercise (CRLF texts, missing `resource`).  Writes only under a temp dir in the scratch area given as argv[1]."""
import os, sys, re, json, shutil, hashlib, importlib.util
TMP = sys.argv[1]; os.makedirs(TMP, exist_ok=True)
K3, K4 = '/home/claude/s1spec/kit_step3_v5', '/home/claude/s1spec/kit_step4b_v5'
R = []
def ok(n, v, d=''): R.append((n, bool(v), d)); print(('PASS ' if v else 'FAIL ') + n, d)
# 1. chart4_shell.py rss() without `resource` (both kits: identical file)
for K in (K3, K4):
    src = open(K + '/writeup/lemmaU/bulk/fix13_v37chain/flux/code/chart4_shell.py').read()
    a = src.index('try:\n    import resource'); a2 = src.index('    resource = None\n', a) + len('    resource = None\n')
    c = src.index('def rss():'); d = src.index('STIP = s1chart.SIG_TIP')
    blk = src[a:a2] + src[c:d]
    sys.modules['resource'] = None            # import resource -> ImportError
    g = {}; exec(blk, g); del sys.modules['resource']
    ok(f'{os.path.basename(K)}: chart4_shell rss() without resource -> -1 on POSIX (no windll)', g['resource'] is None and g['rss']() == -1, str(g['rss']()))
    g2 = {}; exec(blk, g2)
    ok(f'{os.path.basename(K)}: chart4_shell rss() with resource -> MB > 0', g2['resource'] is not None and g2['rss']() > 0, str(g2['rss']()))
# 2. kit3 check.cmp_text: LF identical, CRLF variant, real difference
spec = importlib.util.spec_from_file_location('chk3', K3 + '/check.py'); m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
ref = K3 + '/reference/flux/texts/flux2_lemma_v2.tex'; t = open(ref, 'rb').read()
lf, crlf, bad = (os.path.join(TMP, n) for n in ('lf.tex', 'crlf.tex', 'bad.tex'))
open(lf, 'wb').write(t); open(crlf, 'wb').write(t.replace(b'\n', b'\r\n')); open(bad, 'wb').write(t.replace(b'\n', b'\r\n') + b'x')
m.rows.clear(); m.cmp_text('lf', lf, ref); m.cmp_text('crlf', crlf, ref); m.cmp_text('bad', bad, ref); m.cmp_text('missing', lf + '.no', ref)
ok('kit3 check.cmp_text: LF PASS, CRLF PASS (modulo CRLF), changed FAIL, missing FAIL', [r[0] for r in m.rows] == ['PASS', 'PASS', 'FAIL', 'FAIL'], str(m.rows))
# 3. kit4b verify_v33.compare with a CRLF copy of the v38 reference texts against the LF rerun (test copy of the WSL run)
T = '/home/claude/s1spec/kit4bv5_test_wsl1'; sys.path.insert(0, T)
import verify_v33 as V
ok('kit4b verify_v33 under test is the v5 file', hashlib.sha256(open(V.__file__, 'rb').read()).hexdigest() == hashlib.sha256(open(K4 + '/verify_v33.py', 'rb').read()).hexdigest(), V.__file__)
RD = os.path.join(TMP, 'v38_out_crlf'); shutil.copytree(T + '/reference/v38_out', RD)
conv = []
for dp, _, fs in os.walk(RD):
    for f in fs:
        if not f.endswith('.json'):
            p = os.path.join(dp, f); b_ = open(p, 'rb').read()
            if b'\n' in b_ and b'\r\n' not in b_: open(p, 'wb').write(b_.replace(b'\n', b'\r\n')); conv.append(os.path.relpath(p, RD))
r_lf = {x: (t_, n) for x, t_, n in V.compare(T + '/reference/v38_out', 'FINAL_manifest_v38.json')}
r_cr = {x: (t_, n) for x, t_, n in V.compare(RD, 'FINAL_manifest_v38.json')}
ok('kit4b verify_v33.compare: %d text files converted to CRLF are EXACT modulo CRLF; every other class unchanged' % len(conv),
   conv and all(r_cr[c.replace(os.sep, '/')] == ('EXACT', 'identical modulo CRLF line endings') for c in conv) and all(r_cr[x] == r_lf[x] for x in r_lf if x not in [c.replace(os.sep, '/') for c in conv]),
   json.dumps(sorted(conv)))
ok('kit4b verify_v33.compare on LF: no DIFFERENT text file', not [x for x, (t_, n) in r_lf.items() if t_ == 'DIFFERENT' and not x.endswith('.json')], str([x for x, (t_, n) in r_lf.items() if t_ == 'DIFFERENT']))
# 4. flux_v38 _lf expressions (copied by regex from the v5 file) on CRLF vs LF
fsrc = open(K4 + '/step4/p2/v38/flux_v38.py').read()
lfdef = re.search(r"^_lf = lambda q: .*$", fsrc, re.M).group(0).split('   #')[0]
g = dict(hashlib=hashlib); exec(lfdef, g)
ok('kit4b flux_v38 _lf: sha(LF) == _lf(CRLF) == _lf(LF); differs for a changed file', g['_lf'](crlf) == g['_lf'](lf) == hashlib.sha256(t).hexdigest() != g['_lf'](bad))
print('UNIT SUMMARY: %d PASS, %d FAIL' % (sum(r[1] for r in R), sum(not r[1] for r in R)))
json.dump(R, open(os.path.join(TMP, 'unit_v5.json'), 'w'), indent=1)
sys.exit(0 if all(r[1] for r in R) else 1)
