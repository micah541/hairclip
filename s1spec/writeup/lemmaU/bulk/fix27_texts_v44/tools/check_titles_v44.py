# check_titles_v44.py -- (v44, new; LaTeX compile fix26_latex_v43/LATEX_v43.md defect 2) static check of the moving arguments of the text set.
# For every \part, \chapter, \section, \subsection, \subsubsection and \caption that is NOT starred (starred titles write nothing to the
# .toc/.lof/.aux and take no short title), the title that LaTeX writes to an auxiliary file is the short title [..] if present, else the
# long title {..}.  That written title FAILS if it contains a fragile wrapper (\fn, \key, \pth, \detokenize, \url, \nolinkurl, \verb) or an
# underscore outside math ($...$) that is not escaped as \_.  (In text v43, step3_v43_main.tex:73 had \fn{fix24_texts_v43/lemmaU_v43.tex} in a
# section title; the underscore read back from the .toc stopped pdflatex with "Missing $".)
# Must FAIL on text v43 and PASS on text v44.  Read-only on <text_dir>; writes only <out.json>.
# Usage: PYTHONHASHSEED=0 nice -n 19 python3 -B check_titles_v44.py <text_dir> <out.json>
import re, os, sys, json
DIR, OUT = sys.argv[1].rstrip('/') + '/', sys.argv[2]
CMD = re.compile(r'\\(part|chapter|section|subsection|subsubsection|caption)(\*?)\s*')
FRAG = re.compile(r'\\(fn|key|pth|detokenize|url|nolinkurl|verb)(?![A-Za-z])')
def nocomment(t): return '\n'.join(re.sub(r'(?<!\\)%.*', '', l) for l in t.split('\n'))
def group(t, i, op, cl):
    # t[i] == op; returns (content, index after the closing bracket), brace-aware
    assert t[i] == op
    d = 0; j = i
    while j < len(t):
        c = t[j]
        if c == '\\': j += 2; continue
        if c == '{': d += 1
        elif c == '}': d -= 1
        if op == '[' and c == ']' and d == 0: return t[i + 1:j], j + 1
        if op == '{' and c == '}' and d == 0: return t[i + 1:j], j + 1
        j += 1
    raise ValueError('unclosed')
def fragile(s):
    why = [m.group(0) for m in FRAG.finditer(s)]
    nm = re.sub(r'\$[^$]*\$', ' ', s).replace('\\_', ' ')
    if '_' in nm: why.append('underscore')
    return why
rows = []; nt = 0
for f in sorted(os.listdir(DIR)):
    if not f.endswith('.tex'): continue
    t = nocomment(open(DIR + f, encoding='utf-8', errors='replace').read())
    for m in CMD.finditer(t):
        if t[m.start() - 1:m.start()] == '\\': continue
        i = m.end(); short = None
        if i < len(t) and t[i] == '[': short, i = group(t, i, '[', ']')
        while i < len(t) and t[i] in ' \n': i += 1
        if i >= len(t) or t[i] != '{': continue
        long_, _ = group(t, i, '{', '}')
        nt += 1
        if m.group(2): continue
        written = short if short is not None else long_
        why = fragile(written)
        if why:
            rows.append({'file': f, 'line': t[:m.start()].count('\n') + 1, 'cmd': m.group(1), 'short': short, 'title': long_[:160], 'why': why})
res = {'args': {'argv': sys.argv, 'dir': DIR}, 'titles_seen': nt, 'fail': rows, 'ok': not rows}
json.dump(res, open(OUT, 'w'), indent=1)
for r in rows: print('FAIL title', r['file'], r['line'], r['cmd'], r['why'], r['title'][:100])
print('TITLES', 'PASS' if not rows else 'FAIL', '(%d titles seen)' % nt)
