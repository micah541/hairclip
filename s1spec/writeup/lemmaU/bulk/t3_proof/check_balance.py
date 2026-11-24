#!/usr/bin/env python3
"""Balance check for a LaTeX file (no LaTeX on this machine): braces, environments, inline $, \\[ \\], \\( \\), labels/refs.
usage: python3 -I check_balance.py lemmaU.tex"""
import re, sys
src = open(sys.argv[1]).read()
# strip comments (unescaped %) and verbatim-like \texttt{\detokenize{...}} contents are kept (they are brace-balanced)
lines = [re.sub(r'(?<!\\)%.*$', '', l) for l in src.split('\n')]
txt = '\n'.join(lines)
err = []
# braces
depth = 0
for i, ch in enumerate(txt):
    if ch in '{}' and i > 0 and txt[i - 1] == '\\':
        continue
    if ch == '{': depth += 1
    elif ch == '}':
        depth -= 1
        if depth < 0: err.append('unmatched } at line %d' % (txt[:i].count('\n') + 1)); depth = 0
if depth: err.append('unclosed { : depth %d at end' % depth)
# environments
stack = []
for m in re.finditer(r'\\(begin|end)\{([^}]*)\}', txt):
    ln = txt[:m.start()].count('\n') + 1
    if m.group(1) == 'begin': stack.append((m.group(2), ln))
    else:
        if not stack: err.append('\\end{%s} without begin, line %d' % (m.group(2), ln)); continue
        b, bl = stack.pop()
        if b != m.group(2): err.append('\\begin{%s} (line %d) closed by \\end{%s} (line %d)' % (b, bl, m.group(2), ln))
for b, bl in stack: err.append('unclosed \\begin{%s} line %d' % (b, bl))
# inline math: count unescaped $ per paragraph-independent whole file, after removing \$ ; $$ not used
t2 = txt.replace('\\$', '')
if '$$' in t2: err.append('$$ used')
n = t2.count('$')
if n % 2: err.append('odd number of $: %d' % n)
# per-paragraph $ parity (catches local errors)
for k, para in enumerate(t2.split('\n\n')):
    if para.count('$') % 2:
        err.append('odd $ in paragraph starting line %d' % (t2[:t2.find(para)].count('\n') + 1))
# display math
for a, b in ((r'\[', r'\]'), (r'\(', r'\)')):
    na = len(re.findall(r"(?<!\\)" + re.escape(a), txt)); nb = len(re.findall(r"(?<!\\)" + re.escape(b), txt))
    if na != nb: err.append('%s count %d != %s count %d' % (a, na, b, nb))
# labels / refs
labels = re.findall(r'\\label\{([^}]*)\}', txt)
dup = {l for l in labels if labels.count(l) > 1}
if dup: err.append('duplicate labels: %s' % sorted(dup))
refs = set(re.findall(r'\\(?:ref|eqref)\{([^}]*)\}', txt))
missing = refs - set(labels)
if missing: err.append('undefined refs: %s' % sorted(missing))
print('labels %d, refs %d, $ %d, environments ok=%s' % (len(labels), len(refs), n, not stack))
if err:
    print('FAIL'); [print(' -', e) for e in err]; sys.exit(1)
print('PASS')
