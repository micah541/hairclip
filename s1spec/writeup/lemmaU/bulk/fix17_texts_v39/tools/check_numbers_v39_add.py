# ======== v39 additions (inserted into check_numbers_v39.py before the verdict) ========
# (8) every number in text written NEW for v39 (tools/v39_inserts.json, recorded by the builder) must match the manifest_v38 data
#     (all JSON under step4/p2/out_v38 plus tools/records_v39.json, which holds the manifest results and the record hashes read from
#     disk): NEW-M38.  NEW-CURRENT: matches only another current certified file (reported).  NEW-WHITELIST: tools/whitelist_v39_new.json
#     with a reason.  NEW-MISS: FAIL.  (Integers and 1-digit numbers are not checked.)
m38 = sorted([t for t in cur if t[1] in M38FILES], key=lambda t: t[0]); M38V = [t[0] for t in m38]
WL39 = json.load(open(D + 'tools/whitelist_v39_new.json'))
ins = json.load(open(D + 'tools/v39_inserts.json'))
new_rows = []; newst = {'NEW-M38': 0, 'NEW-CURRENT': 0, 'NEW-WHITELIST': 0, 'NEW-MISS': 0}
for it in ins:
    seen = set()
    for line in it['text'].split('\n'):
        if line.lstrip().startswith('%'): continue
        for m in NUM.finditer(strip_text(line)):
            mant, e1, e2 = m.group(1), m.group(2), m.group(3); k = sig(mant)
            if k < 2 or ('.' not in mant and e1 is None and e2 is None): continue
            p = Dc(mant) * (Dc(10) ** int(e1 or e2 or 0)); txt = m.group(0)
            if (it['file'], txt) in seen: continue
            seen.add((it['file'], txt))
            c = matches(p, k, M38V, m38)
            if c: st = 'NEW-M38'; mt = c[0][1].replace(S, '') + ' ' + c[0][2]
            elif txt in WL39: st = 'NEW-WHITELIST'; mt = WL39[txt]
            else:
                cc = matches(p, k, CV, cur)
                if cc: st = 'NEW-CURRENT'; mt = cc[0][1].replace(S, '') + ' ' + cc[0][2]
                else: st = 'NEW-MISS'; mt = None
            newst[st] += 1; new_rows.append({'file': it['file'], 'num': txt, 'status': st, 'match': mt})
report['v39_new_numbers'] = new_rows
# (9) every sha256 (64 hex) and every abbreviated hash (8-16 hex followed by '...') printed in the v39 texts is looked up in a registry:
#     all hashes in manifest_v38 / manifest_v33 / the closed v3.2 manifest, the on-disk sha256 of every record file (records_v39.json),
#     the chain SHA256SUMS files, and the numbers_A source pins.  In NEW text an unknown hash is a FAIL; in inherited text it is reported.
HEX = re.compile(r'(?<![0-9a-f])([0-9a-f]{64}|[0-9a-f]{8,16}(?=\.\.\.))(?![0-9a-f])')
reg = {}
for f in [S + 'step4/p2/out_v38/FINAL_manifest_v38.json', S + 'step4/p2/out_v33/FINAL_manifest_v33.json', S + 'step4/p2/out/FINAL_manifest.json']:
    for h in set(re.findall(r'[0-9a-f]{64}', open(f).read())): reg.setdefault(h, os.path.basename(f))
R39 = json.load(open(D + 'tools/records_v39.json'))
for k_, x in R39['args']['inputs'].items(): reg[x['sha256']] = 'disk:' + k_
for f in G(B + 'fix13_v37chain/**/SHA256SUMS*') + G(S + 'step3nl_v4c/**/SHA256SUMS*') + G(B + 'fix6_v36/**/SHA256SUMS*') + G(B + 'fix2/**/SHA256SUMS*'):
    try:
        for h in re.findall(r'[0-9a-f]{64}', open(f).read()): reg.setdefault(h, os.path.relpath(f, S))
    except Exception: pass
for x in mj['args']['inputs'].values(): reg.setdefault(x['sha256'], 'numbers_A source')
REGK = sorted(reg)
def lookup(h):
    if len(h) == 64: return reg.get(h)
    i = bisect.bisect_left(REGK, h)
    return reg[REGK[i]] if i < len(REGK) and REGK[i].startswith(h) else None
newtxt = {}
for it in ins: newtxt.setdefault(it['file'], []).append(it['text'])
hrows = []; hst = {'KNOWN': 0, 'INHERITED-UNVERIFIED': 0, 'NEW-UNKNOWN': 0}
for f in TEXTS:
    if not os.path.exists(D + f) or f == 'erratum_v39.tex': continue
    txt = open(D + f).read(); nt = '\n'.join(newtxt.get(f, []))
    for h in sorted(set(HEX.findall(txt))):
        src = lookup(h)
        if src: st = 'KNOWN'
        elif h in nt: st = 'NEW-UNKNOWN'
        else: st = 'INHERITED-UNVERIFIED'
        hst[st] += 1; hrows.append({'file': f, 'hash': h, 'status': st, 'source': src})
report['v39_hashes'] = hrows
# (10) the stale lemmaU_v34 quote 0.1475 must not occur in any v39 text (the erratum and records files are history)
q1475 = [f for f in TEXTS if os.path.exists(D + f) and not f.startswith('erratum') and '0.1475' in open(D + f).read()]
v39ok = newst['NEW-MISS'] == 0 and hst['NEW-UNKNOWN'] == 0 and not q1475
report['v39_summary'] = {'new_numbers': newst, 'hashes': hst, 'files_with_0.1475': q1475, 'm38_ledger_values': len(m38), 'hash_registry': len(reg)}
for r_ in new_rows:
    if r_['status'] == 'NEW-MISS': print('NEW-MISS', r_)
for r_ in hrows:
    if r_['status'] == 'NEW-UNKNOWN': print('NEW-UNKNOWN', r_)
