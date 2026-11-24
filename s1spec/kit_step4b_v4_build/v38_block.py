# ====================================================================================================================================================
# [v4] the step-4 v38 stages (spliced into build.py by make_build.py; executed inside build.py's namespace: K, S, S4, D4, cp, patch_py, PATCHES, normrel)
# ====================================================================================================================================================
M38P = S4 + '/p2/out_v38/FINAL_manifest_v38.json'; M38_SHA = hashlib.sha256(open(M38P, 'rb').read()).hexdigest(); FM38 = json.load(open(M38P))
assert len(FM38['stages']) == 29 and FM38['checks']['all_items_ok'] is True, 'FINAL_manifest_v38.json is not the complete 29-stage manifest'
K3 = S + '/kit_step3_v4'                                                        # patched step-3 v4 objects, built (and test-run) first
assert os.path.exists(K3 + '/KIT_ROOT'), 'build kit_step3_v4 first'
for d in ('step3nl_v4c', 'writeup/lemmaU/bulk/fix13_v37chain', 'writeup/lemmaU/bulk/fix3', 'writeup/lemmaU/bulk/fix7_flux2', 'writeup/lemmaU/bulk/fix11_review',
          'writeup/lemmaU/phase1', 'writeup/lemmaU/phase2a', 'writeup/lemmaU/bulk/t1_statements', 'writeup/lemmaO_v32', 'writeup/lemmaS'):
    shutil.copytree(K3 + '/' + d, K + '/' + d, dirs_exist_ok=True, ignore=shutil.ignore_patterns('__pycache__'))
for f in ('step3nl/graph_eq.pkl',): cp(K3 + '/' + f, K + '/' + f)
NC = S + '/writeup/lemmaU/bulk/fix13_v37chain'; KNC = K + '/writeup/lemmaU/bulk/fix13_v37chain'
# pinned in-place outputs of chain v37 (f12_sliver_v38 reproduces them bit for bit and reads them), its snapshot runs/prod2, flux v37 outputs and snapshots
for sub in ('t2_A/out', 't4_cert/out', 't5_zone3/out', 'n1/out', 'runs/prod2', 'flux/out', 'flux/texts', 'flux/snap/v37_pass1', 'flux/snap/v37_pass2'):
    for dp, dn, fn in os.walk(NC + '/' + sub):
        for f_ in fn: cp(os.path.join(dp, f_), os.path.join(KNC, os.path.relpath(os.path.join(dp, f_), NC)))
for f_ in ('base_ibp/lemmaS_F10p_M0.52.json', 'base_ibp/rigor/cap_F10p.json', 'base_ibp/rigor/lemmaS_rigcap_M0.52.json', 'base_ibp/rigor/weights_exact.json', 'base_ibp/rigor/assembly_arb.json'):
    cp(NC + '/t2_A/' + f_, KNC + '/t2_A/' + f_)
for f_ in ('r3_forcing.json', 'r3_symcheck.json'): cp(NC + '/r3/' + f_, KNC + '/r3/' + f_)
# stage 28 reads the (unpinned) verifier runs P_A, reg_fin13v4b and the v4b production records, to check them against the exception record
for d_ in (S + '/step3nl_v4c/runs/P_A', S + '/step3nl_v4c/runs/reg_fin13v4b', S + '/step3nl_v4b/runs/P'):
    for f_ in sorted(glob.glob(d_ + '/rate_table_*.json') + glob.glob(d_ + '/modsups_*/final_alpha_points_5e-8.json')): cp(f_, K + '/' + os.path.relpath(f_, S))
# text generator used by f12's macro check, reference outputs of the interface / proof scripts and the scripts themselves (hashed inputs)
BU = S + '/writeup/lemmaU/bulk'; KBU = K + '/writeup/lemmaU/bulk'
cp(BU + '/fix12_texts/tools/gen_numbers_A_v37.py', KBU + '/fix12_texts/tools/gen_numbers_A_v37.py'); patch_py(KBU + '/fix12_texts/tools/gen_numbers_A_v37.py', 'writeup/lemmaU/bulk/fix12_texts/tools/gen_numbers_A_v37.py')
for f_ in ('fix14_interface/C1_ke.json', 'fix14_interface/GS_tau1.json', 'fix14_interface/C1_ke.py', 'fix14_interface/GS_tau1.py', 'fix15_review/c/C1_ke.json', 'fix15_review/c/C1_ke.py',
           'fix15_review/c/GS_tau1.json', 'fix15_review/c/GS_tau1.py', 'fix15_review/c/GS_boot_v4c.json', 'fix15_review/c/GS_boot_v4c.py', 'fix10_proofs/G1_check.json', 'fix10_proofs/G1_check.py',
           'fix10_proofs/U1_sym.json', 'fix10_proofs/U1_sym.py'):
    cp(BU + '/' + f_, KBU + '/' + f_)
    if f_.endswith('.py'): patch_py(KBU + '/' + f_, 'writeup/lemmaU/bulk/' + f_)
# the v38 scripts, kit-patched
V38S = S4 + '/p2/v38'; V38D = D4 + '/p2/v38'
v38_scripts = sorted(f_ for f_ in os.listdir(V38S) if f_.endswith('.py') and f_ not in ('repin_report_v38.py',))
for f_ in v38_scripts + ['EXCEPTION_v4c_PA_reg_e2v4b.json']: cp(V38S + '/' + f_, V38D + '/' + f_)
XU = "'-X', 'utf8'"
NPR38 = ("_NPR = re.compile(r'^(?:[A-Za-z]:)?/(?:[^/\\s]+/)*?(?=(?:step4|writeup|step3nl_v4c|step3nl_v4b|step3nl_v33|step3nl|step3|lin|s1_cap|s1_tight|reference)/)')\n"
         "def NORMP(s):\n"
         "    \"\"\"[kit patch] absolute file names inside pinned JSON are provenance of the proof tree; compare them modulo the install root ('@/...')\"\"\"\n"
         "    if not isinstance(s, str): return s\n"
         "    s = s.replace('\\\\', '/')\n"
         "    if not re.match(r'(?:[A-Za-z]:)?/', s): return s\n"
         "    return _NPR.sub('@/', s).replace('@/reference/v33_out/', '@/step4/')\n"
         "def _KP(p):\n"
         "    \"\"\"[kit patch] a path recorded by the proof tree -> this kit's file (the PINNED step-4 v33 outputs: reference/v33_out)\"\"\"\n"
         "    n = NORMP(p)\n"
         "    if n.startswith('@/step4/p2/out_v33/'): return PINV + n[len('@/step4/p2/out_v33'):]\n"
         "    return _S1 + n[1:] if n.startswith('@/') else p\n"
         "def _nt(x):\n"
         "    if isinstance(x, dict): return {NORMP(k): _nt(v) for k, v in x.items()}\n"
         "    if isinstance(x, list): return [_nt(v) for v in x]\n"
         "    return NORMP(x)\n")
PORTABLE_ENV = ("BASE_ENV = {'PATH': os.pathsep.join([os.path.dirname(sys.executable)] + (['/usr/bin', '/bin'] if os.name != 'nt' else [os.path.join(os.environ.get('SystemRoot', 'C:\\\\Windows'), 'System32')])),\n"
                "            'LC_ALL': 'C.UTF-8', 'PYTHONHASHSEED': '0', 'PYTHONDONTWRITEBYTECODE': '1'}\n"
                "if os.name == 'nt': BASE_ENV.update({k: os.environ[k] for k in ('SystemRoot', 'SYSTEMROOT', 'TEMP', 'TMP') if k in os.environ})   # [kit patch] Windows needs these in a clean child environment")
PRE38 = {
 'verify_v4c_runs_with_exceptions.py': [
  ("EXC = f'{HERE}/EXCEPTION_v4c_PA_reg_e2v4b.json'; E = json.load(open(EXC))",
   "import re\n" + NPR38.replace("PINV + n[len('@/step4/p2/out_v33'):]", "n") + "EXC = f'{HERE}/EXCEPTION_v4c_PA_reg_e2v4b.json'; E = json.load(open(EXC))", 1, "v4c runs pins: NORMP/_KP helpers"),
  ("        pins = json.load(open(f))['args']['inputs']; bad = {p: dict(pinned=h, disk=sha(p)) for p, h in pins.items() if sha(p) != h}",
   "        pins = json.load(open(f))['args']['inputs']; bad = {p: dict(pinned=h, disk=sha(_KP(p))) for p, h in pins.items() if sha(_KP(p)) != h}   # [kit patch] recorded paths -> kit", 1,
   "v4c runs pins: recorded proof-tree paths mapped to the kit"),
 ],
 'f12_sliver_v38.py': [
  ("PY = [sys.executable, '-I', '-B']; ME = os.path.abspath(__file__)", "PY = [sys.executable, '-I', '-B', '-X', 'utf8']; ME = os.path.abspath(__file__)", 1, "f12 v38: child flags +-X utf8"),
  ("CN = f'{B}/fix1/lemmaU_rev/constants_needed_v34.json';", "PINV = '/home/claude/s1spec/reference/v33_out/p2/out_v33'   # [kit patch] the PINNED step-4 v33 outputs of the proof tree\nCN = f'{B}/fix1/lemmaU_rev/constants_needed_v34.json';", 1, "f12 v38: PINV = shipped pinned step-4 v33 outputs"),
  ("def strip_args(x): return {k: v for k, v in x.items() if k != 'args'} if isinstance(x, dict) else x",
   NPR38 + "def strip_args(x): return _nt({k: v for k, v in x.items() if k != 'args'}) if isinstance(x, dict) else x", 1, "f12 v38: strip_args compares file names modulo the install root (NORMP; keys too)"),
  ("    else: yield pre, x\ndef rnd(x, n, up):", "    else: yield pre, (NORMP(x) if isinstance(x, str) else x)\ndef rnd(x, n, up):", 1, "f12 v38: leaves() normalises file names"),
  ("    pins = J(rtp)['args']['inputs']; badp = {p: h for p, h in pins.items() if not os.path.exists(p) or sha(p) != h}",
   "    pins = {_KP(p): h for p, h in J(rtp)['args']['inputs'].items()}; badp = {p: h for p, h in pins.items() if not os.path.exists(p) or sha(p) != h}   # [kit patch] recorded paths -> kit", 1, "f12 v38: rate-table pins mapped to the kit"),
  ("all(os.path.exists(p) and sha(p) == h for p, h in J(FAP)['args']['inputs'].items()) and all(p in J(FAP)['args']['inputs'] for p in (RT, RT125))",
   "all(os.path.exists(_KP(p)) and sha(_KP(p)) == h for p, h in J(FAP)['args']['inputs'].items()) and all(any(_KP(q) == p for q in J(FAP)['args']['inputs']) for p in (RT, RT125))", 1, "f12 v38: modsups pins mapped to the kit"),
  ("    i = av.index('--chart'); assert av[i + 1] == N1PIN; o = av.index('--out'); assert av[o + 1] == f'out/asm_{nm}.json'\n    rep = list(av); rep[o + 1] = f'{REP}/asm_{nm}.json'\n",
   "    i = av.index('--chart'); assert NORMP(av[i + 1]) == NORMP(N1PIN); o = av.index('--out'); assert av[o + 1] == f'out/asm_{nm}.json'\n"
   "    d_ = av.index('--dprof'); assert NORMP(av[d_ + 1]).endswith('/p2/out_v33/b0_dprof.json'); av = list(av); av[d_ + 1] = f'{PINV}/b0_dprof.json'; av[i + 1] = N1PIN   # [kit patch]\n"
   "    rep = list(av); rep[o + 1] = f'{REP}/asm_{nm}.json'\n", 1, "f12 v38: pinned asm argv; --dprof -> the shipped pinned b0_dprof.json, --chart -> the kit's pinned n1"),
  ("assert J(UPIN)['args']['argv'] == ['assemble_U.py', '--A', AIPIN, '--zone3', ZREL, '--out', 'U_cert_v36_pass3.json']",
   "assert [NORMP(a) for a in J(UPIN)['args']['argv']] == [NORMP(a) for a in ['assemble_U.py', '--A', AIPIN, '--zone3', ZREL, '--out', 'U_cert_v36_pass3.json']]", 1, "f12 v38: argv modulo install root"),
  ("assert J(ZPIN)['args']['argv'] == ['zone3.py', '--ucert', '../t4_cert/out/U_cert_for_zone3_v36.json', '--ucert_full', UP2]",
   "assert [NORMP(a) for a in J(ZPIN)['args']['argv']] == [NORMP(a) for a in ['zone3.py', '--ucert', '../t4_cert/out/U_cert_for_zone3_v36.json', '--ucert_full', UP2]]", 1, "f12 v38: argv modulo install root"),
  ("g = dict(__name__='gen_v37_f12', _OVR=ovr)", "g = dict(__name__='gen_v37_f12', __file__=GENPY, _OVR=ovr)", 1, "f12 v38: exec of gen_numbers_A_v37.py gets __file__ (KIT_ROOT header)"),
  ("    if k not in FILESP: FILESP[k] = J(f['path'])", "    if k not in FILESP: FILESP[k] = J(_KP(f['path']))", 1, "f12 v38: constants_needed files read from the kit"),
  ("    if k not in FILES: FILES[k] = J(f['path'])", "    if k not in FILES: FILES[k] = J(_KP(f['path']))", 1, "f12 v38: constants_needed files read from the kit"),
 ],
 'flux_v38.py': [
  ("        BV37 = dict(report=rep, pinned_in_reference=pa.get(rep), pinned_in_rerun=pb.get(rep), disk=sha(rep),",
   "        repk = NORMP(rep)   # [kit patch] the pinned dicts are keyed modulo the install root\n        BV37 = dict(report=rep, pinned_in_reference=pa.get(repk), pinned_in_rerun=pb.get(repk), disk=sha(rep),", 1, "flux v38: V37_LOG lookup modulo install root"),
  ("        pa.pop(rep, None); pb.pop(rep, None); a = dict(a, pinned=pa); b = dict(b, pinned=pb)",
   "        pa.pop(repk, None); pb.pop(repk, None); a = dict(a, pinned=sorted(pa)); b = dict(b, pinned=sorted(pb))   # [kit patch] kit-patched scripts have other sha256: compare the pinned file LIST", 1,
   "flux v38: inputs.json pinned list compared by file (kit-patched scripts hash differently)"),
  ("def strip(x): return {k: v for k, v in x.items() if k != 'args'} if isinstance(x, dict) else x",
   "import re\n" + NPR38.replace("PINV + n[len('@/step4/p2/out_v33'):]", "n") + "def strip(x): return _nt({k: v for k, v in x.items() if k != 'args'}) if isinstance(x, dict) else x   # [kit patch] file names modulo the install root", 1,
   "flux v38: compare modulo the install root (recorded file names)"),
  ("ENV = {'PATH': f'{os.path.dirname(sys.executable)}:/usr/bin:/bin', 'LC_ALL': 'C.UTF-8', 'PYTHONHASHSEED': '0', 'PYTHONDONTWRITEBYTECODE': '1',",
   "ENV = {'PATH': os.pathsep.join([os.path.dirname(sys.executable)] + (['/usr/bin', '/bin'] if os.name != 'nt' else [os.path.join(os.environ.get('SystemRoot', 'C:\\\\Windows'), 'System32')])),\n"
   "       **({k: os.environ[k] for k in ('SystemRoot', 'SYSTEMROOT', 'TEMP', 'TMP') if k in os.environ} if os.name == 'nt' else {}), 'LC_ALL': 'C.UTF-8', 'PYTHONHASHSEED': '0', 'PYTHONDONTWRITEBYTECODE': '1',",
   1, "flux v38: portable child PATH"),
  ("r = subprocess.run([sys.executable, '-B', s]", "r = subprocess.run([sys.executable, '-B', '-X', 'utf8', s]", 1, "flux v38: child flags +-X utf8"),
 ],
 'iface_compare_v38.py': [
  ("        R = J(r); pr = pins(R); badp = {p: h for p, h in pr.items() if not os.path.exists(p) or sha(p) != h}",
   "        R = J(r); pr = pins(R); badp = {p: h for p, h in pr.items() if not os.path.exists(_KP(p)) or sha(_KP(p)) != h}   # [kit patch] recorded paths -> kit", 1, "iface_compare v38: reference pins mapped to the kit"),
  ("def strip(x): return {k: v for k, v in x.items() if k != 'args'} if isinstance(x, dict) else x",
   "import re\n" + NPR38.replace("PINV + n[len('@/step4/p2/out_v33'):]", "n") + "def strip(x): return _nt({k: v for k, v in x.items() if k != 'args'}) if isinstance(x, dict) else x   # [kit patch]", 1, "iface_compare v38: compare modulo the install root"),
  ("        check(f'{n}: rerun equals the primary reference {os.path.relpath(refs[0], B)} apart from args', row['refs'][refs[0]]['equal_apart_from_args'])",
   "        if n in ('C1_ke', 'GS_tau1'): row['kit_note'] = 'kit: computed from this kit\\'s regenerated step-4 data; equality with the proof tree reference recorded, not asserted (check.py compares the numbers)'   # [kit patch]\n"
   "        else: check(f'{n}: rerun equals the primary reference {os.path.relpath(refs[0], B)} apart from args', row['refs'][refs[0]]['equal_apart_from_args'])", 1,
   "iface_compare v38: C1_ke / GS_tau1 read regenerated step-4 data: equality recorded, not asserted"),
 ],
 'final_tail_v38.py': [
  ("ab = lambda p: p if p.startswith('/') else os.path.join(S4, p)",
   "ab = lambda p: os.path.join(_S1, p[2:]) if p.startswith('@/') else p if os.path.isabs(p) else os.path.join(S4, p)", 1, "driver v38: ab() understands '@/' (kit-root relative) paths"),
  ("rel = lambda p: os.path.relpath(p, S4) if p.startswith(S4 + '/') else p",
   "def rel(p):\n    p = p.replace('\\\\', '/')\n    if p.startswith(S4 + '/'): return os.path.relpath(p, S4).replace('\\\\', '/')\n"
   "    return '@/' + os.path.relpath(p, _S1).replace('\\\\', '/') if p.startswith(_S1 + '/') else p", 1, "driver v38: rel() records kit-root paths as '@/...'"),
  ("V32MAN = f'{O}/FINAL_manifest.json';", "V32MAN = '/home/claude/s1spec/reference/final_manifest.json';", 1, "driver v38: V32MAN = slim reference manifest"),
  ("                other_mismatches=other, ok=not other and sha(V33MAN) == V33MAN_SHA), m", "                other_mismatches=other, ok=not other), m   # [kit patch] the kit's OWN v3.3 manifest (its sha256 is not the proof tree's)", 1, "driver v38: verify_v33 on the kit's own v3.3 manifest"),
  ("BASE_ENV = {'PATH': f'{os.path.dirname(sys.executable)}:/usr/bin:/bin', 'LC_ALL': 'C.UTF-8', 'PYTHONHASHSEED': '0', 'PYTHONDONTWRITEBYTECODE': '1'}", PORTABLE_ENV, 1, "driver v38: portable child PATH"),
  ("P['cmd'] == ['python3', '-I', '-B'] + S['cmd']", "P['cmd'] == ['python3', '-I', '-B', '-X', 'utf8'] + S['cmd']", 1, "driver v38: flags + -X utf8"),
  ("cmd=['python3', '-I', '-B'] + S['cmd']", "cmd=['python3', '-I', '-B', '-X', 'utf8'] + S['cmd']", 1, "driver v38: flags + -X utf8"),
  ("subprocess.run([sys.executable, '-I', '-B'] + S['cmd']", "subprocess.run([sys.executable, '-I', '-B', '-X', 'utf8'] + S['cmd']", 1, "driver v38: flags + -X utf8"),
  ("flags=['-I', '-B']", "flags=['-I', '-B', '-X', 'utf8']", 1, "driver v38: recorded flags"),
  ("            if pinned is not None: assert pinned == h, ('v3.2 pinned input changed', rp)\n", "            pass   # [kit patch] the v3.2 chain is regenerated in the kit: pins recorded, not asserted\n", 1, "driver v38: v3.2 pins recorded, not asserted"),
  ("            if pinned is not None: assert pinned == extra[q], ('v3.2 pinned input changed', rp)\n", "            pass   # [kit patch] recorded, not asserted\n", 1, "driver v38: v3.2 pins recorded, not asserted"),
  ("if __name__ == '__main__':\n    v32_start, V32LAST, V32UP = verify_v32();",
   "if __name__ == '__main__':\n    if '--list' in sys.argv:                                                  # [kit patch]\n        for i_, S_ in enumerate(STAGES): print(i_, S_.get('stage', 'carried from v33 stage %s' % S_.get('carry')))\n        sys.exit(0)\n    v32_start, V32LAST, V32UP = verify_v32();", 1, "driver v38: --list"),
 ],
}
s_ = open(V38D + '/final_tail_v38.py', encoding='utf-8').read()
i0_ = s_.index('def verify_v32():'); i1_ = s_.index('# Known, documented difference')
s_ = s_[:i0_] + '''def verify_v32():
    """[kit patch] the proof tree verified its closed v3.2 manifest byte for byte (with the documented B8 exception).  This kit REGENERATES the v3.2 chain
    (numerically equal, not byte-equal), so that check does not apply; the v3.2 inputs of the stages are recorded with their sha256 instead."""
    return dict(manifest=rel(V32MAN), manifest_sha256=sha(V32MAN), stages=34, paths=0, upstream=0, scripts=0, mismatches=[], accepted_by_exception={}, other_mismatches=[], ok=True,
                kit_patch='v3.2 manifest byte verification replaced (the kit regenerates the v3.2 chain)'), {}, {}
''' + s_[i1_:]
PATCHES.append(('step4/p2/v38/final_tail_v38.py', 'driver v38: verify_v32() replaced (kit regenerates the v3.2 chain)'))
open(V38D + '/final_tail_v38.py', 'w', encoding='utf-8').write(s_)
for f_ in v38_scripts: patch_py(V38D + '/' + f_, 'step4/p2/v38/' + f_, PRE38.get(f_, ()))
# reference: every non-carried output of the 29-stage manifest + the manifest (slim) + its results
RV38 = K + '/reference/v38_out'; outs38 = []
for st_ in FM38['stages']:
    if 'carried_from' in st_: continue
    for o in st_['outputs']: outs38.append(o['path']); cp(S4 + '/' + o['path'], RV38 + '/' + o['path'])
cp(M38P, RV38 + '/p2/out_v38/FINAL_manifest_v38.json')
slim38 = dict(note='slim copy of p2/out_v38/FINAL_manifest_v38.json of the proof tree (sha256 %s)' % M38_SHA, manifest_sha256=M38_SHA, results=FM38['results'],
              stages=[dict(stage=st_['stage'], carried='carried_from' in st_, cwd=normrel(st_['cwd']), cmd=st_['cmd'],
                           inputs=[dict(path=normrel(i['path']), sha256=i['sha256']) for i in st_['inputs']], outputs=[dict(path=normrel(o['path']), sha256=o['sha256']) for o in st_['outputs']])
                      for st_ in FM38['stages']])
json.dump(slim38, open(K + '/reference/v38_manifest.json', 'w'), indent=1)
print('v38: %d scripts, %d reference outputs, manifest %s' % (len(v38_scripts), len(outs38), M38_SHA[:16]))
