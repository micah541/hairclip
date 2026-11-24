# make_numsub_v38.py -- turns the STALE/SUSPECT suggestions of check_numbers_v38.py (tools/stale_suggestions_v38.json) into replacement
# rows for the builders (tools/numsub_v38.json), keeping the printed style (\cdot10^{e} or decimal), the number of significant digits and
# the rounding direction of the old literal.  Rows already in numsub_v38.json are kept; manual rows (key 'manual') are never overwritten.
# Every row records the JSON file and key of the new value.  The output is reviewed by hand before the builders are rerun.
# Usage: python3 -B make_numsub_v38.py
import json, re, os
from decimal import Decimal as Dc, ROUND_CEILING, ROUND_FLOOR
D = '/home/claude/s1spec/writeup/lemmaU/bulk/fix16_texts_v38/'
SUG = json.load(open(D + 'tools/stale_suggestions_v38.json'))
NS = json.load(open(D + 'tools/numsub_v38.json'))
SRCNAME = {'continuation_v38.tex': 'continuation_v38.tex'}
def fmt(old, newraw, up):
    m = re.match(r'(\d+\.?\d*)(\s*\\cdot\s*10\^\{?(-?\d+)\}?)?$', old)
    mant, e = m.group(1), int(m.group(3) or 0)
    k = len(mant.replace('.', '').lstrip('0'))
    v = Dc(newraw)
    if m.group(2):
        q = Dc(1).scaleb(v.adjusted() - k + 1); r = v.quantize(q, rounding=ROUND_CEILING if up else ROUND_FLOOR)
        ex = r.adjusted(); digs = ''.join(map(str, r.as_tuple().digits)).ljust(k, '0')[:k]
        return (digs[0] + ('.' + digs[1:] if k > 1 else '')) + '\\cdot10^{%d}' % ex
    dec = len(mant.split('.')[1]) if '.' in mant else 0
    q = Dc(1).scaleb(-dec); r = v.quantize(q, rounding=ROUND_CEILING if up else ROUND_FLOOR)
    return format(r, 'f')
SKIP = {('continuation_v38.tex', '7.12\\cdot10^{-3}'), ('main_v3_v38.tex', '10.7'), ('step3_v38_main.tex', '10.7'), ('lemmaS_A_v38.tex', '6.25'),
        ('main_v3_v38.tex', '2.69\\cdot10^{-3}'), ('main_v3_v38.tex', '9.872\\cdot10^{-10}'), ('main_v3_v38.tex', '6.867\\cdot10^{-7}')}   # false key matches / records (whitelisted)
added = 0
for s in SUG:
    if s.get('none') or not s.get('new_raw'): continue
    f = s['file']; old = s['old_printed']
    if (f, old) in SKIP: continue
    new = fmt(old, s['new_raw'], s['dir'] == 'up')
    if new == old: continue
    rows = NS.setdefault(f, [])
    if any(r['old'] == old for r in rows): continue
    txt = open(D + f).read()
    rows.append({'old': old, 'new': new, 'all': txt.count(old) > 1, 'key': s['source'], 'dir': s['dir'], 'old_raw': s['old_raw'], 'new_raw': s['new_raw']})
    added += 1
json.dump(NS, open(D + 'tools/numsub_v38.json', 'w'), indent=1)
print('added', added)
for f, rows in NS.items():
    for r in rows: print(f, '|', r['old'], '->', r['new'], '| all' if r.get('all') else '', '|', r.get('key', r.get('manual', ''))[-70:])
