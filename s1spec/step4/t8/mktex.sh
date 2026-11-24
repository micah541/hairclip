#!/bin/sh
# Build LEMMA_F.tex/.pdf from LEMMA_F.md (pandoc + xelatex, FreeSerif); underscores and * escaped against markdown parsing.
set -e
cd "$(dirname "$0")"
python3 - <<'PY'
import re
s = open('LEMMA_F.md').read()
s = s.replace('_', r'\_')
s = re.sub(r'(?<=[xτy])\*', r'\\*', s)
s = s.replace(r'\_*', r'\_\*')
out = []; pat = re.compile(r'^(\d+\.|-) ')
for ln in s.split('\n'):
    if pat.match(ln) and out and out[-1].strip() and not pat.match(out[-1]) and not out[-1].startswith(' '):
        out.append('')
    out.append(ln)
s = '\n'.join(out)
open('/tmp/lemmaF_esc.md', 'w').write(s)
PY
pandoc -f markdown-superscript-subscript-tex_math_dollars /tmp/lemmaF_esc.md -s -o LEMMA_F.tex \
  -V mainfont=FreeSerif -V geometry:margin=1in -V title="Lemma F: far-field energies of the P2 error (step 4, T8)" --pdf-engine=xelatex
for i in 1 2; do xelatex -interaction=nonstopmode -halt-on-error LEMMA_F.tex > /tmp/lemmaF_xelatex.log 2>&1; done
echo "missing glyphs: $(grep -c 'Missing character' /tmp/lemmaF_xelatex.log)"
grep -o "Output written on .*" /tmp/lemmaF_xelatex.log
rm -f LEMMA_F.aux LEMMA_F.out
