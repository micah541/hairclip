# b42_common.py -- shared helpers of the text-v42 builder (v42: copy of fix21_texts_v41/tools/b41_common.py with v42 names and new
# TITLE_REPS (incl. the v41 leftover title of step3_main); logic unchanged) (fix22_texts_v42/tools).  The builder reads the text-v41 files of
# fix21_texts_v41/ (never edits them), retags the v41 FILE NAMES of the text set to their v42 names, applies asserted replacements and
# writes NEW files into fix22_texts_v42/.  The word "v41" itself is not retagged (sentences about what text v41 did keep it).
#  rep(o, n): o must occur exactly once (count asserted).  rep_all(o, n, k): exactly k >= 1 times.  Inserted text is recorded in
#  tools/v42_inserts.json, so that check_numbers_v42.py can require every number in new text to match certified data.
import os, re, sys, json
S = '/home/claude/s1spec/'
B = S + 'writeup/lemmaU/bulk/'
D9 = B + 'fix21_texts_v41/'
D = B + 'fix22_texts_v42/'
# stems of the text-set files; <stem>_v41 -> <stem>_v42 (also with \_ escapes); step3_v41_main -> step3_v42_main
STEMS = ['bootstrap', 'C1_identification', 'E0_LEMMA_v4c', 'ERRATUM', 'G1_export_k2', 'GS_lemma', 'RECORDS_ERRATA', 'U1_regularity',
         'b12_preamble', 'chart_lemma_body', 'continuation', 'ei_constants', 'ei_lemma_v4c', 'erratum', 'flux2_lemma', 'lemmaS_A', 'lemmaU',
         'main_v3', 'numbers_A', 'step4_handoff', 'STEP4_PATCH_NOTE']
def v42name(f):
    if f == 'step3_v41_main.tex': return 'step3_v42_main.tex'
    return re.sub(r'_v41(\.\w+)$', r'_v42\1', f)
def retag(t):
    for st in sorted(STEMS, key=len, reverse=True):
        for s in (st, st.replace('_', r'\_')):
            sep = r'\_' if '\\' in s else '_'
            t = re.sub(r'(?<![A-Za-z0-9])' + re.escape(s + sep + 'v41') + r'(?![0-9])', lambda m: m.group(0)[:-3] + 'v42', t)
    t = re.sub(r'(?<![A-Za-z0-9])step3_v41_main(?![0-9])', 'step3_v42_main', t)
    t = t.replace(r'step3\_v41\_main', r'step3\_v42\_main')
    # paths into fix19 that name a file of the set now point to its v42 successor in fix21
    def path(m):
        base = m.group(1)
        return 'fix22_texts_v42/' + base if re.match(r'(%s|step3)_v42' % '|'.join(map(re.escape, STEMS)), base) else m.group(0)
    t = re.sub(r'fix21_texts_v41/([\w.\\-]+)', path, t)
    return t
# v42: applied to \title / \date lines only.  Dates of text v41 become "(text v42; text v41 ...", titles name text v42.
TITLE_REPS = [('21 November 2026 (text v41', '23 November 2026 (text v42; text v41'), ('21 November 2026 (v41:', '23 November 2026 (v42: Lean replaces DKS (review W7); v41:'),
              ('Step 3, text v41:', 'Step 3, text v42:'), ('(step 3, text v41)', '(step 3, text v42)'), ('(text v41)', '(text v42)'),
              ('Erratum v41', 'Erratum v42'), ('v41, November 21, 2026', 'v42, November 23, 2026'),
              ('proof text, revision v41', 'proof text, revision v42'), ('Draft v3, revision v41 of November 21, 2026', 'Draft v3, revision v42 of November 23, 2026'),
              # v42: the v41 leftover (TEXTS_REPORT_v41 section 4): the title of step3_v41_main.tex still said "text v39"
              ('Step 3 of the hair-clip proof, text v39 (', 'Step 3 of the hair-clip proof, text v42 (')]
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
        self.src = D9 + name; self.name = v42name(name); self.dst = D + self.name
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
    json.dump(INSERTS, open(D + 'tools/v42_inserts.json', 'w'), indent=1)
