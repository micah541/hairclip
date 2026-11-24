# b38_common.py -- shared helpers of the text-v38 builders (fix16_texts_v38/tools).  Every builder reads an existing file (never edits
# it), applies asserted replacements and writes a NEW file into fix16_texts_v38/.
#  rep(o, n)      : o must occur exactly once in the current text (count asserted).
#  rep_all(o, n)  : o must occur at least once; all occurrences replaced.
#  GEN            : generic retag map for files whose references are to the CURRENT objects (applied after the specific
#                   replacements; substrings listed in `protect` are shielded).  Historical files (erratum) do not use GEN.
#  \u00ab...\u00bb       : protection markers for intentional historical references in new text (removed on write).
#  HASHES         : old sha256 (fix6_v36 / fix8 objects) -> new sha256 (fix13_v37chain objects), computed from disk at build time.
import hashlib, os, re, sys
S = '/home/claude/s1spec/'
B = S + 'writeup/lemmaU/bulk/'
D = B + 'fix16_texts_v38/'
def sha(p): return hashlib.sha256(open(p, 'rb').read()).hexdigest()
PAIRS = [('fix6_v36/t2_A/out/A_interface.json', 'fix13_v37chain/t2_A/out/A_interface.json'),
         ('fix6_v36/t2_A/out/asm_design_N2ref4.json', 'fix13_v37chain/t2_A/out/asm_design_N2ref4.json'),
         ('fix6_v36/t4_cert/out/U_cert_v36_pass2.json', 'fix13_v37chain/t4_cert/out/U_cert_v36_pass2.json'),
         ('fix6_v36/t4_cert/out/U_cert_v36_pass3.json', 'fix13_v37chain/t4_cert/out/U_cert_v36_pass3.json'),
         ('fix6_v36/t4_cert/out/u1_ghat.json', 'fix13_v37chain/t4_cert/out/u1_ghat.json'),
         ('fix6_v36/t4_cert/out/u34_primary.json', 'fix13_v37chain/t4_cert/out/u34_primary.json'),
         ('fix6_v36/t4_cert/out/U_cert_for_zone3_v36.json', 'fix13_v37chain/t4_cert/out/U_cert_for_zone3_v36.json'),
         ('fix6_v36/t5_zone3/out/zone3_interface.json', 'fix13_v37chain/t5_zone3/out/zone3_interface.json'),
         ('fix8_flux2_impl/out/flux2_cert.json', 'fix13_v37chain/flux/out/flux2_cert.json')]
HASHES = []
for o, n in PAIRS:
    ho, hn = sha(B + o), sha(B + n)
    if ho != hn: HASHES += [(ho, hn), (ho[:16], hn[:16]), (ho[:8] + '...', hn[:8] + '...')]
GEN = [('step3nl_v4b/V4B_REPORT.md', 'step3nl_v4c/V4C_REPORT.md'), ('step3nl_v4b', 'step3nl_v4c'), ('V4B_REPORT.md', 'V4C_REPORT.md'),
       ('fin13v4b', 'fin13v4c'),
       ('fix8_flux2_impl/out/', 'fix13_v37chain/flux/out/'), ('fix8_flux2_impl/code/', 'fix13_v37chain/flux/code/'),
       ('fix6_v36/', 'fix13_v37chain/'),
       ('numbers_A_v37', 'numbers_A_v38'), ('b12_preamble_v37', 'b12_preamble_v38'), ('chart_lemma_body_v37', 'chart_lemma_body_v38'),
       ('ei_constants_v37', 'ei_constants_v38'), ('lemmaS_A_v37', 'lemmaS_A_v38'), ('flux2_lemma_v37', 'flux2_lemma_v38'),
       ('continuation_v37', 'continuation_v38'), ('lemmaU_v37', 'lemmaU_v38'), ('main_v3_v37', 'main_v3_v38'), ('erratum_v37', 'erratum_v38'),
       ('step3_v37_main', 'step3_v38_main'), ('ei_lemma_v4.tex', 'ei_lemma_v4c.tex'), ('{ei_lemma_v4}', '{ei_lemma_v4c}'),
       ('fix8_flux2_impl/texts/flux2_lemma_v2.tex', 'fix16_texts_v38/flux2_lemma_v38.tex'), ('flux2_lemma_v2.tex', 'flux2_lemma_v38.tex'),
       ('fix10_proofs/U1_regularity.tex', 'fix16_texts_v38/U1_regularity_v38.tex'), ('U1_regularity.tex', 'U1_regularity_v38.tex'),
       ('fix10_proofs/G1_export_k2.tex', 'fix16_texts_v38/G1_export_k2_v38.tex'), ('G1_export_k2.tex', 'G1_export_k2_v38.tex'),
       ('fix14_interface/C1_identification.tex', 'fix16_texts_v38/C1_identification_v38.tex'), ('C1_identification.tex', 'C1_identification_v38.tex'),
       ('fix14_interface/GS_lemma.tex', 'fix16_texts_v38/GS_lemma_v38.tex'), ('GS_lemma.tex', 'GS_lemma_v38.tex'),
       ('step3nl_v4c/e0_lemma/E0_LEMMA_v4b.tex', 'fix16_texts_v38/E0_LEMMA_v4c_v38.tex'), ('E0_LEMMA_v4b.tex', 'E0_LEMMA_v4c_v38.tex'),
       ('fix12_texts/', 'fix16_texts_v38/'), ('text v37', 'text v38'), ('(v37)', '(v38)')]
GLOBAL_PROTECT = ['fix6_v36/V36_REPORT.md', '% Changes against v34: numbers/hashes of chain v36 (verifier fin13v4b, R3 forcing constant, level-2 export r_2);',
                  "% Changes: numbers and citations only (verifier fin13v4b, chain v36, flux lemma v2 / (I3'), c_dU), a v37 status box", '% text version v37 (chain v36 on verifier fin13v4b).']
class Builder:
    def __init__(self, src, dst, generic=True, protect=()):
        self.src = src; self.dst = dst; self.generic = generic; self.protect = list(protect)
        self.t = open(src).read(); self.n = 0
    def rep(self, o, n):
        c = self.t.count(o)
        if c != 1: sys.exit('%s: replacement count %d for: %r' % (os.path.basename(self.dst), c, o[:140]))
        self.t = self.t.replace(o, n); self.n += 1
    def rep_all(self, o, n):
        c = self.t.count(o)
        if c < 1: sys.exit('%s: no occurrence of: %r' % (os.path.basename(self.dst), o[:140]))
        self.t = self.t.replace(o, n); self.n += 1
    def between(self, a, b):
        i = self.t.index(a); j = self.t.index(b, i); return self.t[i:j]
    def write(self):
        t = self.t
        ph = {}
        # text between the markers \u00ab ... \u00bb is protected from GEN/HASHES (intentional historical references); markers are removed
        for i, m in enumerate(re.findall('\u00ab(.*?)\u00bb', t, flags=re.S)):
            k = '\x00K%03d\x00' % i; ph[k] = m; t = t.replace('\u00ab' + m + '\u00bb', k, 1)
        # comment lines are provenance: never retagged
        for i, m in enumerate(re.findall(r'(?m)^[ \t]*%.*$', t)):
            k = '\x00C%04d\x00' % i; ph[k] = m; t = t.replace(m, k, 1)
        if self.generic:
            for i, p in enumerate(self.protect + [x for x in GLOBAL_PROTECT if x in t]):
                if p not in t: sys.exit('%s: protected string absent: %r' % (os.path.basename(self.dst), p[:100]))
                k = '\x00P%03d\x00' % i; ph[k] = p; t = t.replace(p, k)
            for o, n in HASHES + GEN: t = t.replace(o, n)
        for k, p in ph.items(): t = t.replace(k, p)
        assert '\u00ab' not in t and '\u00bb' not in t, 'unbalanced protection marker'
        open(self.dst, 'w').write(t)
        left = sorted(set(m.group(0) for m in re.finditer(r'fin13v4b|fix6_v36|fix8_flux2_impl|step3nl_v4b|_v37\b|v37\.tex|\bv36\b', t)))
        print('wrote %s (%d specific replacements); remaining tags: %s' % (os.path.basename(self.dst), self.n, left))
