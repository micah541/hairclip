# b43_common.py -- shared helpers of the text-v43 builder (v43: copy of fix22_texts_v42/tools/b42_common.py with v43 names and new
# TITLE_REPS; logic unchanged) (fix24_texts_v43/tools).  The builder reads the text-v42 files of fix22_texts_v42/ (never edits them),
# retags the v42 FILE NAMES of the text set to their v43 names, applies asserted replacements and writes NEW files into fix24_texts_v43/.
# The word "v42" itself is not retagged (sentences about what text v42 did keep it).
#  rep(o, n): o must occur exactly once (count asserted).  rep_all(o, n, k): exactly k >= 1 times.  Inserted text is recorded in
#  tools/v43_inserts.json, so that check_numbers_v43.py can require every number in new text to match certified data.
import os, re, sys, json
S = '/home/claude/s1spec/'
B = S + 'writeup/lemmaU/bulk/'
D9 = B + 'fix22_texts_v42/'
D = B + 'fix24_texts_v43/'
# stems of the text-set files; <stem>_v42 -> <stem>_v43 (also with \_ escapes); step3_v42_main -> step3_v43_main
STEMS = ['bootstrap', 'C1_identification', 'E0_LEMMA_v4c', 'ERRATUM', 'G1_export_k2', 'GS_lemma', 'RECORDS_ERRATA', 'U1_regularity',
         'b12_preamble', 'chart_lemma_body', 'continuation', 'ei_constants', 'ei_lemma_v4c', 'erratum', 'flux2_lemma', 'lemmaS_A', 'lemmaU',
         'main_v3', 'numbers_A', 'step4_handoff', 'STEP4_PATCH_NOTE']
def v43name(f):
    if f == 'step3_v42_main.tex': return 'step3_v43_main.tex'
    return re.sub(r'_v42(\.\w+)$', r'_v43\1', f)
def retag(t):
    for st in sorted(STEMS, key=len, reverse=True):
        for s in (st, st.replace('_', r'\_')):
            sep = r'\_' if '\\' in s else '_'
            t = re.sub(r'(?<![A-Za-z0-9])' + re.escape(s + sep + 'v42') + r'(?![0-9])', lambda m: m.group(0)[:-3] + 'v43', t)
    t = re.sub(r'(?<![A-Za-z0-9])step3_v42_main(?![0-9])', 'step3_v43_main', t)
    t = t.replace(r'step3\_v42\_main', r'step3\_v43\_main')
    # paths into fix22 that name a file of the set now point to its v43 successor in fix24
    def path(m):
        base = m.group(1)
        return 'fix24_texts_v43/' + base if re.match(r'(%s|step3)_v43' % '|'.join(map(re.escape, STEMS)), base) else m.group(0)
    t = re.sub(r'fix22_texts_v42/([\w.\\-]+)', path, t)
    return t
# v43: applied to \title / \date lines only.
TITLE_REPS = [('23 November 2026 (text v42;', '23 November 2026 (text v43; text v42;'),
              ('23 November 2026 (v42: Lean replaces DKS (review W7);', '23 November 2026 (v43: delta check fix23 S-4; v42: Lean replaces DKS (review W7);'),
              ('Step 3, text v42:', 'Step 3, text v43:'), ('(step 3, text v42)', '(step 3, text v43)'), ('(text v42)', '(text v43)'),
              ('Erratum v42', 'Erratum v43'), ('v42, November 23, 2026', 'v43, November 23, 2026'),
              ('proof text, revision v42', 'proof text, revision v43'), ('Draft v3, revision v42 of November 23, 2026', 'Draft v3, revision v43 of November 23, 2026'),
              ('Step 3 of the hair-clip proof, text v42 (', 'Step 3 of the hair-clip proof, text v43 (')]
def retitle(t):
    out = []
    for line in t.split('\n'):
        if re.match(r'\s*\\(title|date|author\{\}\\date)', line):
            for a, b in TITLE_REPS: line = line.replace(a, b)
        out.append(line)
    return '\n'.join(out)
INSERTS = []
class Builder:
    def __init__(self, name, retag_names=True, header=None):
        self.src = D9 + name; self.name = v43name(name); self.dst = D + self.name
        self.t = open(self.src).read(); self.n = 0
        if retag_names: self.t = retag(self.t)
        self.t = retitle(self.t)
        if header: self.t = header.rstrip('\n') + '\n' + self.t
    def rep(self, o, n, new=True):
        c = self.t.count(o)
        if c != 1: sys.exit('%s: replacement count %d for: %r' % (self.name, c, o[:200]))
        self.t = self.t.replace(o, n); self.n += 1
        if new: INSERTS.append({'file': self.name, 'text': n})
    def rep_all(self, o, n, k, new=True):
        c = self.t.count(o)
        if c != k: sys.exit('%s: expected %d occurrences, found %d: %r' % (self.name, k, c, o[:200]))
        self.t = self.t.replace(o, n); self.n += 1
        if new: INSERTS.append({'file': self.name, 'text': n})
    def write(self):
        open(self.dst, 'w').write(self.t)
        print('wrote %s (%d replacements)' % (self.name, self.n))
        return self.dst
def save_inserts():
    json.dump(INSERTS, open(D + 'tools/v43_inserts.json', 'w'), indent=1)
