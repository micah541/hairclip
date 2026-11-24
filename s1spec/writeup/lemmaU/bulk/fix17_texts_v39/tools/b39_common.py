# b39_common.py -- shared helpers of the text-v39 builder (fix17_texts_v39/tools).  The builder reads the text-v38 files of
# fix16_texts_v38/ (never edits them), retags the v38 FILE NAMES of the text set to their v39 names, applies asserted replacements and
# writes NEW files into fix17_texts_v39/.  Sentences that say what was new in v38 keep "v38" (no generic retag of the word v38).
#  rep(o, n): o must occur exactly once (count asserted).  rep_all(o, n): at least once.  ins(...) records inserted v39 text
#  (tools/v39_inserts.json) so that check_numbers_v39.py can require every number in new text to match manifest_v38 data.
import hashlib, os, re, sys, json
S = '/home/claude/s1spec/'
B = S + 'writeup/lemmaU/bulk/'
D8 = B + 'fix16_texts_v38/'
D = B + 'fix17_texts_v39/'
def sha(p): return hashlib.sha256(open(p, 'rb').read()).hexdigest()
# v38 file -> v39 file
FMAP = {'bootstrap_v38.tex': 'bootstrap_v39.tex', 'C1_identification_v38.tex': 'C1_identification_v39.tex',
        'E0_LEMMA_v4c_v38.tex': 'E0_LEMMA_v4c_v39.tex', 'ERRATUM_v38.md': 'ERRATUM_v39.md', 'G1_export_k2_v38.tex': 'G1_export_k2_v39.tex',
        'GS_lemma_v38.tex': 'GS_lemma_v39.tex', 'RECORDS_ERRATA_v38.md': 'RECORDS_ERRATA_v39.md', 'U1_regularity_v38.tex': 'U1_regularity_v39.tex',
        'b12_preamble_v38.tex': 'b12_preamble_v39.tex', 'chart_lemma_body_v38.tex': 'chart_lemma_body_v39.tex',
        'continuation_v38.tex': 'continuation_v39.tex', 'ei_constants_v38.tex': 'ei_constants_v39.tex', 'ei_constants_v38.json': 'ei_constants_v39.json',
        'ei_lemma_v4c.tex': 'ei_lemma_v4c_v39.tex', 'erratum_v38.tex': 'erratum_v39.tex', 'flux2_lemma_v38.tex': 'flux2_lemma_v39.tex',
        'lemmaS_A_v38.tex': 'lemmaS_A_v39.tex', 'lemmaU_v38.tex': 'lemmaU_v39.tex', 'main_v3_v38.tex': 'main_v3_v39.tex',
        'numbers_A_v38.tex': 'numbers_A_v39.tex', 'numbers_A_v38.json': 'numbers_A_v39.json', 'step3_v38_main.tex': 'step3_v39_main.tex'}
def _stems():
    out = []
    for a, b in FMAP.items():
        sa, sb = os.path.splitext(a)[0], os.path.splitext(b)[0]
        out.append((sa, sb))
    out.sort(key=lambda t: -len(t[0]))
    return out
STEMS = _stems()
def retag(t):
    """v38 file names of the text set -> v39 names (with and without \\_ escapes); fix16_texts_v38/<mapped> -> fix17_texts_v39/<mapped>."""
    ph = {}
    # paths into fix16 that are NOT files of the set (reports, tools) stay as they are
    def path(m):
        rest = m.group(1)
        base = rest.split('/')[0]
        if base in FMAP or os.path.splitext(base)[0] in dict(STEMS):
            return 'fix17_texts_v39/' + rest
        k = '\x00F%04d\x00' % len(ph); ph[k] = m.group(0); return k
    t = re.sub(r'fix16_texts_v38/([\w./\\-]*)', path, t)
    for sa, sb in STEMS:
        if sa == 'ei_lemma_v4c':
            t = re.sub(r'ei_lemma_v4c(?=\.tex|\})', 'ei_lemma_v4c_v39', t)
            t = re.sub(r'ei\\_lemma\\_v4c(?=\.tex|\})', r'ei\\_lemma\\_v4c\\_v39', t)
            continue
        t = re.sub(r'(?<![\w])' + re.escape(sa) + r'(?![\w])', sb, t)
        ea, eb = sa.replace('_', r'\_'), sb.replace('_', r'\_')
        t = t.replace(ea, eb)
    for k, v in ph.items(): t = t.replace(k, v)
    return t
INSERTS = []
class Builder:
    def __init__(self, name, retag_names=True, src=None):
        self.name = name; self.dst = D + FMAP.get(name, name)
        self.src = src or (D8 + name)
        self.t = open(self.src).read(); self.n = 0
        if retag_names: self.t = retag(self.t)
    def rep(self, o, n, new=True):
        c = self.t.count(o)
        if c != 1: sys.exit('%s: replacement count %d for: %r' % (os.path.basename(self.dst), c, o[:160]))
        self.t = self.t.replace(o, n); self.n += 1
        if new: INSERTS.append({'file': os.path.basename(self.dst), 'text': n})
    def rep_all(self, o, n, new=True):
        c = self.t.count(o)
        if c < 1: sys.exit('%s: no occurrence of: %r' % (os.path.basename(self.dst), o[:160]))
        self.t = self.t.replace(o, n); self.n += 1
        if new: INSERTS.append({'file': os.path.basename(self.dst), 'text': n})
    def between(self, a, b):
        i = self.t.index(a); j = self.t.index(b, i); return self.t[i:j]
    def head(self, line):
        self.t = line.rstrip('\n') + '\n' + self.t
    def write(self):
        open(self.dst, 'w').write(self.t)
        left = sorted(set(m.group(0) for m in re.finditer(r'[\w\\/.]*v38[\w\\/.]*', self.t)))
        print('wrote %s (%d replacements); remaining v38 tokens: %s' % (os.path.basename(self.dst), self.n, left))
def save_inserts():
    json.dump(INSERTS, open(D + 'tools/v39_inserts.json', 'w'), indent=1)
