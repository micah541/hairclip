# check_macros_v44.py -- (v44: copy of fix24_texts_v43/tools/check_macros_v43.py with one fix, LaTeX compile fix26 defect 1.)
# v44 FIX: the v40-v43 check removed every \fn{...}, \key{...}, \pth{...}, \texttt{...}, ... argument TOGETHER WITH the macro name before
#  collecting the used control words, so an undefined wrapper macro was never seen (\key was used but undefined in GS_lemma, erratum and
#  flux2_lemma of text v43; pdflatex stops there).  Now the names of the wrappers that occur (fn, key, pth, texttt, detokenize, nolinkurl, url)
#  are added to USED before their arguments are removed.  The check stays per standalone document (expanded through its \input's).
#  Must FAIL on text v43 (\key in GS_lemma, erratum, flux2_lemma) and PASS on text v44.
# Usage: PYTHONHASHSEED=0 nice -n 19 python3 -B check_macros_v44.py <text_dir> <tag v39|...|v44> <out.json>
# check_macros_v43.py -- (v43: copy of fix22_texts_v42/tools/check_macros_v42.py; on code lines fix22_texts_v42 -> fix24_texts_v43 and v42 -> v43; comments below are the v42 history, unchanged; v43 additions marked 'v43'.)
# check_macros_v42.py -- (v42: copy of fix19_texts_v40/tools/check_macros_v40.py, logic unchanged) "used subset of defined" macro check for the step-3 text set (review fix18 F-2), plus a valid-name check.
# LaTeX is not installed, so the check is static:
#  (1) VALID NAMES: every name defined by \newcommand, \renewcommand, \providecommand (with or without braces, starred or not),
#      \def, \let, \DeclareMathOperator must be letters only (or a single non-letter).  v39 emitted \providecommand{\BPhiC0oneV}, which
#      TeX reads as \BPhiC followed by the text 0oneV.
#  (2) USED SUBSET OF DEFINED, per standalone document: the document is expanded through its \input's (comments removed; inside
#      step3_*_main the \ifdefined\stepmain\else ... \fi blocks of the \input files are removed, because \stepmain is defined there).
#      USED = every control word \name (letters) outside \fn{}, \key{}, \pth{}, \texttt{}, \detokenize{}, \nolinkurl{}, \url{}.
#      DEFINED = the names defined in that expansion.  A used, undefined name FAILS if
#        (a) it is defined somewhere else in the text set (a project macro missing from this document, e.g. \cdUV in main_v3), or
#        (b) it is not in the list KNOWN of LaTeX-kernel / package control words below.
#      KNOWN is a fixed list (the control words of LaTeX2e, amsmath/amssymb/amsthm, enumitem, booktabs, longtable, geometry, hyperref,
#      tikz, xcolor, mathtools, array, graphicx and fontenc that the set uses); a name is added to it only after checking that it is a
#      genuine kernel/package command.  So every undefined name of the set is either a known command or a FAIL.
# Usage: PYTHONHASHSEED=0 nice -n 19 python3 -B check_macros_v42.py <text_dir> <tag v39|v40|v42> <out.json>
# Read-only on <text_dir>; writes only <out.json>.
import re, os, sys, json
DIR, TAG, OUT = sys.argv[1].rstrip('/') + '/', sys.argv[2], sys.argv[3]
DOCS = [f % TAG for f in ['step3_%s_main.tex', 'lemmaU_%s.tex', 'main_v3_%s.tex', 'U1_regularity_%s.tex', 'C1_identification_%s.tex',
                          'GS_lemma_%s.tex', 'G1_export_k2_%s.tex', 'E0_LEMMA_v4c_%s.tex', 'erratum_%s.tex', 'step4_handoff_%s.tex',
                          'bootstrap_%s.tex', 'continuation_%s.tex', 'ei_lemma_v4c_%s.tex', 'flux2_lemma_%s.tex', 'lemmaS_A_%s.tex']]
KNOWN = set('''
documentclass usepackage begin end title author date maketitle section subsection subsubsection paragraph label ref eqref pageref cite
bibitem item emph textbf textit texttt textsf textup textrm textsc mathrm mathbf mathcal mathfrak mathbb mathsf mathit operatorname
newcommand renewcommand providecommand def let newtheorem theoremstyle DeclareMathOperator input include setcounter addtocounter
newcounter stepcounter arabic roman alph Alph the thesection thetheorem small footnotesize scriptsize tiny normalsize large Large LARGE
huge bf it rm tt sf sc em centering noindent indent par newline linebreak newpage clearpage vspace hspace quad qquad medskip bigskip
smallskip hfill vfill hline toprule midrule bottomrule cmidrule endhead endfirsthead endfoot endlastfoot multicolumn fbox parbox mbox
makebox raisebox linewidth textwidth hsize vsize columnwidth baselineskip ensuremath detokenize nolinkurl url href hypersetup
ldots cdots dots dotsc dotsb vdots ddots cdot times le ge leq geq ne neq approx sim simeq equiv propto subset subseteq supset
supseteq in notin ni cup cap setminus emptyset varnothing forall exists infty partial nabla sum prod int iint oint lim limsup liminf sup
inf max min arg det exp log ln sin cos tan sinh cosh tanh arctan sqrt frac tfrac dfrac binom tbinom left right big Big bigg Bigg bigl
bigr Bigl Bigr biggl biggr Biggl Biggr langle rangle lvert rvert lVert rVert vert Vert mid lfloor rfloor lceil rceil to mapsto
rightarrow leftarrow Rightarrow Leftarrow Leftrightarrow leftrightarrow longrightarrow Longrightarrow longmapsto uparrow downarrow
Uparrow Downarrow nearrow searrow hookrightarrow iff implies impliedby pm mp ast star circ bullet oplus otimes wedge vee neg lnot
alpha beta gamma delta epsilon varepsilon zeta eta theta vartheta iota kappa lambda mu nu xi pi varpi rho varrho sigma varsigma tau
upsilon phi varphi chi psi omega Gamma Delta Theta Lambda Xi Pi Sigma Upsilon Phi Psi Omega ell hbar Re Im prime
hat bar tilde widetilde widehat overline underline vec dot ddot check breve acute grave mathring overbrace underbrace stackrel
overset underset text intertext displaystyle textstyle scriptstyle limits nolimits qedhere qed square blacksquare Box
tag notag nonumber label displaybreak allowdisplaybreaks align gather split cases aligned matrix pmatrix bmatrix
LaTeX TeX today and thanks protect relax ifdefined else fi colon ldotp ss aa AA ae oe o O l L i j S P dag ddag copyright
textbackslash textasciitilde textasciicircum textunderscore textbar textless textgreater textendash textemdash
quad enspace thinspace negthinspace hphantom vphantom phantom smash kern mkern mskip hskip vskip
arraystretch tabcolsep extrarowheight renewenvironment newenvironment enumerate itemize description
setlength addtolength textcolor color colorbox definecolor tikzset draw node fill path coordinate foreach tikzstyle
includegraphics resizebox scalebox rotatebox footnote footnotemark footnotetext
ignorespaces unskip strut null nobreak allowbreak hyphenation sloppy fussy
lesssim gtrsim leqslant geqslant ll gg prec succ preceq succeq perp parallel angle triangle nabla
mathop mathbin mathrel mathord mathpunct mathopen mathclose
DeclarePairedDelimiter coloneqq eqqcolon
numberwithin tableofcontents appendix bibliographystyle bibliography
smallsetminus backslash lbrace rbrace lbrack rbrack
textdegree checkmark
caption refstepcounter setlist thetable varkappa address emergencystretch raggedbottom rule texorpdfstring usetikzlibrary
cong dim bigcup pmod geometry hrule part verb
'''.split())
def read(name):
    p = DIR + (name if name.endswith('.tex') else name + '.tex')
    return open(p, encoding='utf-8', errors='replace').read() if os.path.exists(p) else None
def nocomment(t):
    return '\n'.join(re.sub(r'(?<!\\)%.*', '', l) for l in t.split('\n'))
SKIP = re.compile(r'\\ifdefined\\stepmain\\else(.*?)\\fi', re.S)
missing = []
def expand(name, stepmain, depth=0):
    t = read(name)
    if t is None: missing.append(name); return ''
    t = nocomment(t)
    if stepmain:
        for m in SKIP.finditer(t): assert '\\if' not in m.group(1), (name, 'nested \\if in a stepmain block')
        t = SKIP.sub('', t)
    sm = stepmain or bool(re.search(r'\\def\\stepmain\b', t))
    return re.sub(r'\\input\{([^}]+)\}', lambda m: expand(m.group(1), sm, depth + 1), t)
DEF_PATS = [r'\\(?:re|provide|new)command\*?\s*\{?\s*\\([^\s{}\[\]\\]+)', r'\\(?:re|provide|new)command\*?\s*\{?\s*\\(\\)',
            r'\\def\s*\\([^\s{}#\\]+)', r'\\let\s*\\([^\s{}=\\]+)', r'\\DeclareMathOperator\*?\s*\{\s*\\([^\s{}]+)\}',
            r'\\DeclarePairedDelimiter\s*\{?\s*\\([^\s{}]+)']
def defs(t):
    out = set()
    for p in DEF_PATS: out |= set(re.findall(p, t))
    return out
def valid(n): return bool(re.fullmatch(r'[A-Za-z]+', n)) or (len(n) == 1 and not n.isalpha())
VERB = re.compile(r'\\(?:fn|key|pth|texttt|detokenize|nolinkurl|url)\s*\{(?:[^{}]|\{[^{}]*\})*\}')
VERBI = re.compile(r'\\verb\*?([^A-Za-z\s])(.*?)\1')   # \verb|...| : its content is verbatim
def used(t):
    t = t.replace('\\\\', ' ')                       # \\ (line break) is not the control word of the letters after it
    t = re.sub(r'\\ifdefined\s*\\[A-Za-z]+', ' ', t)     # \ifdefined\X tests X; it does not use it
    t = VERBI.sub(r'\\verb ', t)
    wr = set(re.findall(r'\\(fn|key|pth|texttt|detokenize|nolinkurl|url)\s*\{', t))   # v44: the wrapper names are used too
    t = VERB.sub(' ', t)
    return set(re.findall(r'\\([A-Za-z]+)', t)) | wr
# project macros: every name defined in any .tex file of the set
PROJECT = set(); bad_names = []
for f in sorted(os.listdir(DIR)):
    if not f.endswith('.tex'): continue
    t = nocomment(open(DIR + f, encoding='utf-8', errors='replace').read())
    d = defs(t); PROJECT |= {n for n in d if valid(n)}
    bad_names += [{'file': f, 'name': n} for n in sorted(d) if not valid(n)]
res = {'args': {'argv': sys.argv, 'dir': DIR, 'tag': TAG}, 'invalid_names': bad_names, 'docs': {}}
ok = not bad_names
for doc in DOCS:
    if read(doc) is None: res['docs'][doc] = {'missing': True}; ok = False; continue
    missing.clear()
    t = expand(doc, False)
    d = defs(t); u = used(t)
    und = sorted(u - d)
    proj = [n for n in und if n in PROJECT]
    unknown = [n for n in und if n not in PROJECT and n not in KNOWN]
    r = {'used': len(u), 'defined': len(d), 'undefined_project_macros': proj, 'undefined_unknown': unknown,
         'missing_inputs': sorted(set(missing)), 'ok': not proj and not unknown and not missing}
    res['docs'][doc] = r; ok = ok and r['ok']
res['ok'] = ok
json.dump(res, open(OUT, 'w'), indent=1, sort_keys=True)
for b in bad_names: print('FAIL invalid macro name', b)
for doc, r in res['docs'].items():
    print(('PASS' if r.get('ok') else 'FAIL'), 'macros', doc, {k: v for k, v in r.items() if k in ('undefined_project_macros', 'undefined_unknown', 'missing_inputs', 'missing') and v})
print('MACROS', 'PASS' if ok else 'FAIL')
