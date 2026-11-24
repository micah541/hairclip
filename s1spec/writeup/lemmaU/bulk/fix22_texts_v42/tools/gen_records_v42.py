# gen_records_v42.py -- reads the re-pinned step-4 records and the v4 kits FROM DISK (read-only) and writes tools/records_v42.json:
# (v42: copy of fix19_texts_v40/tools/gen_records_v40.py with output path v42; logic unchanged.)
# full sha256 of the manifests, driver, report, exception records and kits; stage list; pin counts; the result items quoted by the
# v40 texts (v40: copy of fix17_texts_v39/tools/gen_records_v39.py plus the normalised distinct-file count, review fix18 R-1); the check.py SUMMARY lines of the two WSL kit runs.  Everything the builder prints about records comes from this file.
# Usage: PYTHONHASHSEED=0 nice -n 19 python3 -B gen_records_v42.py
import json, hashlib, os, re, sys, platform
S = '/home/claude/s1spec/'; P = S + 'step4/p2/'
D = S + 'writeup/lemmaU/bulk/fix22_texts_v42/'
def sha(p): return hashlib.sha256(open(p, 'rb').read()).hexdigest()
F = {'manifest_v38': P + 'out_v38/FINAL_manifest_v38.json', 'driver_v38': P + 'v38/final_tail_v38.py', 'repin_report': P + 'out_v38/REPIN_REPORT.md',
     'v38_log': P + 'v38/V38_LOG.md', 'manifest_v33': P + 'out_v33/FINAL_manifest_v33.json', 'manifest_v32': P + 'out/FINAL_manifest.json',
     'exc_B8': P + 'v33/EXCEPTION_v32_B8_check_outA5.json', 'exc_v4c': P + 'v38/EXCEPTION_v4c_PA_reg_e2v4b.json',
     'verify_v32_exc': P + 'v33/verify_v32_with_exceptions.py', 'verify_v4c_exc': P + 'v38/verify_v4c_runs_with_exceptions.py',
     'cn34': S + 'writeup/lemmaU/bulk/fix1/lemmaU_rev/constants_needed_v34.json', 'cn34_pre': S + 'writeup/lemmaU/bulk/fix1/lemmaU_rev/constants_needed_v34_pre_f34.json',
     'f34_report': S + 'writeup/lemmaU/bulk/fix2/F34_REPORT.md',
     'kit3': S + 'kit_step3_v4.zip', 'kit4b': S + 'kit_step4b_v4.zip',
     'kit3_notes': S + 'kit3v4_run_wsl5/RUN_NOTES.md', 'kit4b_notes': S + 'kit4bv4_run_wsl1/RUN_NOTES.md',
     'n1_v38': P + 'out_v38/n1_chart.json', 'n1_chain': S + 'writeup/lemmaU/bulk/fix13_v37chain/n1/out/n1_chart_fin13v4c.json',
     'C1_ke_v38': P + 'out_v38/iface/C1_ke.json', 'GS_tau1_v38': P + 'out_v38/iface/GS_tau1.json', 'GS_boot_v38': P + 'out_v38/iface/GS_boot_v4c.json',
     'G1_v38': P + 'out_v38/iface/G1_check.json', 'U1_v38': P + 'out_v38/iface/U1_sym.json', 'iface_cmp': P + 'out_v38/iface/iface_compare_v38.json',
     'flux_v38': P + 'out_v38/flux/flux_v38.json', 'f12_v38': P + 'out_v38/f12/f12_sliver.json', 'pins_v4c': P + 'out_v38/v4c_runs_pins.json',
     'step4_proof': S + 'writeup/step4_proof/step4_proof.tex'}
R = {'args': {'argv': sys.argv, 'python': platform.python_version(), 'inputs': {k: {'path': v, 'sha256': sha(v), 'bytes': os.path.getsize(v)} for k, v in F.items()}}}
m = json.load(open(F['manifest_v38']))
R['stages'] = [s['stage'] for s in m['stages']]
ins = sum(len(s['inputs']) for s in m['stages']); outs = sum(len(s['outputs']) for s in m['stages'])
R['pins'] = {'inputs': ins, 'outputs': outs, 'scripts': len(m['script_sha256']), 'entries': ins + outs + len(m['script_sha256']),
             'hashes_verified_at_end': m['checks']['hashes_verified_at_end'], 'all_items_ok': m['checks']['all_items_ok']}
R['checks'] = {k: m['checks'][k] for k in ('v32_manifest_at_end', 'v33_manifest_at_end')}
R['results'] = m['results']
# every pinned (path, sha256) of manifest_v38 (registry for the hash check)
pins = {}
for s in m['stages']:
    for x in s['inputs'] + s['outputs']:
        pins[x['path']] = x['sha256']
for k, v in m['script_sha256'].items(): pins[k] = v
R['pinned'] = pins
# v40 (review fix18 R-1): distinct files after normalising relative (to step4/) and absolute spellings of the pinned paths
def _norm(q): return q if q.startswith('/') else os.path.normpath(S + 'step4/' + q)
allp = [x['path'] for s in m['stages'] for x in s['inputs'] + s['outputs']] + list(m['script_sha256'])
R['pins']['distinct_files_normalised'] = len(set(_norm(q) for q in allp))
R['pins']['distinct_spellings'] = len(set(allp))
_sp = set(x['path'] for s in m['stages'] for x in s['inputs'] + s['outputs']); _sc = set(m['script_sha256'])
R['pins']['stage_paths_distinct'] = len(_sp); R['pins']['scripts_also_stage_paths'] = len(_sp & _sc)
assert R['pins']['distinct_files_normalised'] == len(_sp | _sc) == R['pins']['stage_paths_distinct'] + R['pins']['scripts'] - R['pins']['scripts_also_stage_paths']
def summ(p, pat):
    t = open(p).read(); mm = re.search(pat, t); return mm.group(0) if mm else None
R['kit3_summary'] = summ(F['kit3_notes'], r'SUMMARY: \d+ PASS, \d+ FAIL, \d+ WARN, \d+ SKIP')
R['kit4b_summary'] = summ(F['kit4b_notes'], r'TOTAL PASS\s+\[\d+ PASS, \d+ FAIL, \d+ WARN\]')
assert R['kit3_summary'] == 'SUMMARY: 80 PASS, 0 FAIL, 1 WARN, 0 SKIP', R['kit3_summary']
assert R['kit4b_summary'].split('[')[1] == '159 PASS, 0 FAIL, 1 WARN]', R['kit4b_summary']
assert len(m['stages']) == 29 and R['pins']['all_items_ok'] is True
ip = m['results']['interface_proofs']
R['derived'] = {'C1_ke_margin_half': 0.5 / ip['C1_ke']['sup_ke'], 'GS_tau1_margin_ku_half': 0.5 / ip['GS_tau1']['sup_ku'],
                'GS_tau1_margin_uprime_one': 1.0 / ip['GS_tau1']['sup_uprime'],
                'note': 'margins = threshold / certified sup, in float (printed rounded down to 3 digits, which absorbs the float error)'}
json.dump(R, open(D + 'tools/records_v42.json', 'w'), indent=1, sort_keys=True)
print('records_v42.json: manifest_v38', R['args']['inputs']['manifest_v38']['sha256'], 'stages', len(R['stages']), 'pins', R['pins']['entries'],
      'kit3', R['args']['inputs']['kit3']['sha256'][:8], 'kit4b', R['args']['inputs']['kit4b']['sha256'][:8])
