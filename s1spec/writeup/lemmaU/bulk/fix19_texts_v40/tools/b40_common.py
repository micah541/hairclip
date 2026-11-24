# b40_common.py -- shared helpers of the text-v40 builder (fix19_texts_v40/tools).  The builder reads the text-v39 files of
# fix17_texts_v39/ (never edits them), retags the v39 FILE NAMES of the text set to their v40 names, applies asserted replacements and
# writes NEW files into fix19_texts_v40/.  The word "v39" itself is not retagged (sentences about what text v39 did keep it).
#  rep(o, n): o must occur exactly once (count asserted).  rep_all(o, n, k): exactly k >= 1 times.  Inserted text is recorded in
#  tools/v40_inserts.json, so that check_numbers_v40.py can require every number in new text to match certified data.
import os, re, sys, json
S = '/home/claude/s1spec/'
B = S + 'writeup/lemmaU/bulk/'
D9 = B + 'fix17_texts_v39/'
D = B + 'fix19_texts_v40/'
# stems of the text-set files; <stem>_v39 -> <stem>_v40 (also with \_ escapes); step3_v39_main -> step3_v40_main
STEMS = ['bootstrap', 'C1_identification', 'E0_LEMMA_v4c', 'ERRATUM', 'G1_export_k2', 'GS_lemma', 'RECORDS_ERRATA', 'U1_regularity',
         'b12_preamble', 'chart_lemma_body', 'continuation', 'ei_constants', 'ei_lemma_v4c', 'erratum', 'flux2_lemma', 'lemmaS_A', 'lemmaU',
         'main_v3', 'numbers_A', 'step4_handoff', 'STEP4_PATCH_NOTE']
def v40name(f):
    if f == 'step3_v39_main.tex': return 'step3_v40_main.tex'
    return re.sub(r'_v39(\.\w+)$', r'_v40\1', f)
def retag(t):
    for st in sorted(STEMS, key=len, reverse=True):
        for s in (st, st.replace('_', r'\_')):
            sep = r'\_' if '\\' in s else '_'
            t = re.sub(r'(?<![A-Za-z0-9])' + re.escape(s + sep + 'v39') + r'(?![0-9])', lambda m: m.group(0)[:-3] + 'v40', t)
    t = re.sub(r'(?<![A-Za-z0-9])step3_v39_main(?![0-9])', 'step3_v40_main', t)
    t = t.replace(r'step3\_v39\_main', r'step3\_v40\_main')
    # paths into fix17 that name a file of the set now point to its v40 successor in fix19
    def path(m):
        base = m.group(1)
        return 'fix19_texts_v40/' + base if re.match(r'(%s|step3)_v40' % '|'.join(map(re.escape, STEMS)), base) else m.group(0)
    t = re.sub(r'fix17_texts_v39/([\w.\\-]+)', path, t)
    return t
TITLE_REPS = [('20 November 2026 (text v39', '21 November 2026 (text v40; text v39'), ('text v39:', 'text v40:'), ('(text v39)', '(text v40)'),
              ('revision v39', 'revision v40'), ('Erratum v39', 'Erratum v40'), ('v39, November 20, 2026', 'v40, November 21, 2026'),
              ('text v39,', 'text v40,'), ('Draft v3, revision v40 of November 20, 2026', 'Draft v3, revision v40 of November 21, 2026')]
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
        self.src = D9 + name; self.name = v40name(name); self.dst = D + self.name
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
    json.dump(INSERTS, open(D + 'tools/v40_inserts.json', 'w'), indent=1)
