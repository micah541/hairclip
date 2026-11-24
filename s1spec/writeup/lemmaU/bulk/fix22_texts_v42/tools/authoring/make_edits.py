import sys, re
sys.path.insert(0, '.')
from new_U1 import *
from new_C1 import *
from new_rest import E as REST
ST = '/tmp/claude-1000/-home-claude/1625b396-7223-44a0-a7db-c079a426b542/scratchpad/stage/'
def v42(f): return f.replace('step3_v41_main', 'step3_v42_main').replace('_v41.', '_v42.')
EDITS = [
 ('U1_regularity_v41.tex', U1_PRE_OLD, U1_PRE_NEW),
 ('U1_regularity_v41.tex', U1_LABEL_OLD, U1_LABEL_NEW),
 ('U1_regularity_v41.tex', '@FROM@' + U1_A_START + '@TO@' + U1_A_END, U1_A_NEW),
 ('U1_regularity_v41.tex', '@FROM@' + U1_B_START + '@TO@' + U1_B_END, U1_B_NEW),
 ('U1_regularity_v41.tex', U1_GAP1_OLD, U1_GAP1_NEW),
 ('U1_regularity_v41.tex', '@FROM@' + U1_GAP2_START + '@TO@' + U1_GAP3_END, U1_GAP23_NEW),
 ('U1_regularity_v41.tex', U1_GAP6_OLD, U1_GAP6_NEW),
 ('C1_identification_v41.tex', C1_PRE_OLD, C1_PRE_NEW),
 ('C1_identification_v41.tex', C1_LAB_OLD, C1_LAB_NEW),
 ('C1_identification_v41.tex', C1_CL_OLD, C1_CL_NEW),
 ('C1_identification_v41.tex', '@FROM@' + C1_I4_START + '@TO@' + C1_I4_END, C1_I4_NEW),
 ('C1_identification_v41.tex', '@FROM@' + C1_I5_START + '@TO@' + C1_I5_END, C1_I5_NEW),
 ('C1_identification_v41.tex', C1_MAXPF_OLD, C1_MAXPF_NEW),
 ('C1_identification_v41.tex', C1_V_OLD, C1_V_NEW),
 ('C1_identification_v41.tex', '@FROM@' + C1_GAP_START + '@TO@' + C1_GAP_END, C1_GAP_NEW),
] + REST
texts = {}
out = {}; order = []
for f, old, new in EDITS:
    if f not in texts: texts[f] = open(ST + v42(f)).read(); order.append(f); out[f] = []
    t = texts[f]
    if old.startswith('@LINE@'):
        pre = old[6:]; lines = [l for l in t.split('\n') if l.startswith(pre)]
        assert len(lines) == 1, (f, pre, len(lines)); old = lines[0]
    elif old.startswith('@FROM@'):
        a, b = old[6:].split('@TO@'); assert t.count(a) == 1, (f, a); i = t.index(a); j = t.index(b, i); assert t.count(b) == 1 or t[i:].count(b) >= 1, (f, b)
        old = t[i:j + len(b)]
    assert t.count(old) == 1, (f, old[:120], t.count(old))
    assert '\n@@ ' not in '\n' + new and '\n@@ ' not in '\n' + old
    texts[f] = t.replace(old, new); out[f].append((old, new))
hdr = """## edits_v42.txt -- the asserted text replacements of text v42, read by build_v42.py (plain text: no escaping).
## Source: writeup/lean_bridge/REVIEW_W7.md section 3 (Lean replaces DKS; items 1-12), the hand bridges S3a/S3b (one macro \\LeanSthree per file),
## the v41 leftovers (TEXTS_REPORT_v41 section 4) and the open-item lists.  No number is computed.
## Block format:  "@@ FILE <v41 source name>" starts a file section; "@@ REP" (exactly once) starts an edit; the old text follows up to
## "@@ NEW"; the new text up to "@@ END".  Old texts are matched AFTER the retag of the set's file names to v42 and the retitle of
## \\title/\\date lines (b42_common.py), and after the earlier edits of the same file (applied in order).
"""
s = hdr
for f in order:
    s += '\n@@ FILE %s\n' % f
    for old, new in out[f]: s += '@@ REP\n' + old + '\n@@ NEW\n' + new + '\n@@ END\n'
open(sys.argv[1], 'w').write(s)
print('edits', sum(len(v) for v in out.values()), 'files', order)
