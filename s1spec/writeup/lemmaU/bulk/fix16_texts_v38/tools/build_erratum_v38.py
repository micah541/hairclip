# build_erratum_v38.py -- builds ERRATUM_v38.md and erratum_v38.tex from fix12_texts/ERRATUM_v37.md and fix12_texts/erratum_v37.tex
# (unchanged; kept verbatim as Parts 1-2) and tools/erratum_part3_v38.md (Part 3, sections K-N; converted to TeX by a small converter:
# headings, lists, `code` -> \fn{}, **bold**, <= -> \le, -> -> \to; other special characters escaped).
# Usage: PYTHONHASHSEED=0 nice -n 19 python3 -B build_erratum_v38.py
import re, sys
sys.path.insert(0, '/home/claude/s1spec/writeup/lemmaU/bulk/fix16_texts_v38/tools')
from b38_common import B, D
T = B + 'fix12_texts/'
p3 = open(D + 'tools/erratum_part3_v38.md').read()

# ---------- Markdown ----------
md = open(T + 'ERRATUM_v37.md').read()
head = ("# Erratum v38 to step 3 v3.1 / v3.2 / v3.3, the step-3 verifiers (fin11, fin12f4, fin13v4, fin13v4b) and step 4\n\n"
        "**v38, Nov 20, 2026** (`writeup/lemmaU/bulk/fix16_texts_v38/`, built by `tools/build_erratum_v38.py`). Parts 1 and 2 below are "
        "`fix12_texts/ERRATUM_v37.md` verbatim (its first line, the v37 title, is kept as a heading); Part 3 (sections K-N) is new. "
        "Where Part 3 settles or supersedes a statement of Part 2 (sections G, I, J), Part 3 governs.\n\n")
md = head + md.replace('# Erratum v37 to step 3', '## (v37 title) Erratum v37 to step 3', 1).rstrip() + '\n' + p3
open(D + 'ERRATUM_v38.md', 'w').write(md)

# ---------- TeX ----------
def esc(t):
    out = []; parts = re.split(r'(`[^`]*`)', t)
    for x in parts:
        if x.startswith('`') and x.endswith('`') and len(x) > 1:
            out.append('\\fn{' + x[1:-1] + '}'); continue
        x = x.replace('\\', '\\textbackslash{}')
        for a, b in [('&', '\\&'), ('%', '\\%'), ('#', '\\#'), ('_', '\\_'), ('^', '\\^{}'), ('$', '\\$'), ('~', '\\~{}')]:
            x = x.replace(a, b)
        x = x.replace('{}{}', '{}')
        x = re.sub(r'\*\*(.+?)\*\*', r'\\textbf{\1}', x)
        x = x.replace('<=', '$\\le$').replace('>=', '$\\ge$').replace('->', '$\\to$')
        out.append(x)
    return ''.join(out)
lines = p3.split('\n'); tex = []; stack = []
def close_lists(level):
    while len(stack) > level: tex.append('\\end{%s}' % stack.pop())
for ln in lines:
    if ln.strip() == '---': close_lists(0); continue
    m = re.match(r'^(#+)\s+(.*)$', ln)
    if m:
        close_lists(0); lvl = len(m.group(1)); txt = esc(m.group(2))
        tex.append({1: '\\part*{%s}', 2: '\\section*{%s}', 3: '\\subsection*{%s}'}.get(lvl, '\\paragraph{%s}') % txt); continue
    m = re.match(r'^(\s*)(- |\d+\. )(.*)$', ln)
    if m:
        depth = len(m.group(1)) // 2 + 1; kind = 'itemize' if m.group(2) == '- ' else 'enumerate'
        close_lists(depth)
        if len(stack) < depth: tex.append('\\begin{%s}' % kind); stack.append(kind)
        tex.append('\\item ' + esc(m.group(3))); continue
    if ln.strip() == '': close_lists(0); tex.append(''); continue
    close_lists(0); tex.append(esc(ln))
close_lists(0)
p3tex = '\n'.join(tex)
t = open(T + 'erratum_v37.tex').read()
for o, n in [("\\title{Erratum v37 to step~3 v3.1/v3.2/v3.3, the step-3 verifiers (fin11, fin12f4, fin13v4) and step~4}",
              "\\title{Erratum v38 to step~3 v3.1/v3.2/v3.3, the step-3 verifiers (fin11, fin12f4, fin13v4, fin13v4b) and step~4}"),
             ("\\date{v37, November 20, 2026 (first version November 13, 2026)}", "\\date{v38, November 20, 2026 (first version November 13, 2026)}"),
             ("\\noindent\\textbf{v37 (November 20, 2026; \\fn{writeup/lemmaU/bulk/fix12_texts/}).}",
              "\\noindent\\textbf{v38 (November 20, 2026; \\fn{writeup/lemmaU/bulk/fix16_texts_v38/}, built by \\fn{tools/build_erratum_v38.py}).} Parts~1 and~2 are \\fn{fix12_texts/erratum_v37.tex} verbatim; Part~3 (sections K--N, at the end) is new and governs where it settles or supersedes sections G, I, J of Part~2.\n\n\\noindent\\textbf{v37 (November 20, 2026; \\fn{writeup/lemmaU/bulk/fix12_texts/}).}"),
             ("\\end{document}", p3tex + "\n\\paragraph{Files (v38).} \\fn{fix16_texts_v38/erratum_v38.tex} (this file) and \\fn{fix16_texts_v38/ERRATUM_v38.md} (same entries); balance checked by script; not compiled.\n\\end{document}")]:
    assert t.count(o) == 1, o[:60]; t = t.replace(o, n)
open(D + 'erratum_v38.tex', 'w').write(t)
print('wrote ERRATUM_v38.md, erratum_v38.tex')
