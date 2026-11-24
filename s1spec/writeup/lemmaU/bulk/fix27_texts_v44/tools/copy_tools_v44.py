# copy_tools_v44.py -- one-off: copies the generic text-set tools of fix24_texts_v43/tools/ (read-only) to fix27_texts_v44/tools/ with v44
# names.  On code lines (not starting with '#'): 'fix24_texts_v43' -> 'fix27_texts_v44', then 'v43' -> 'v44'.  Comment lines are kept (history),
# and a first line names the source.  Data files keymap/whitelist: the file keys '<stem>_v43.tex' -> '_v44.tex'.
# check_balance: NAMEDIRS also get fix24_texts_v43/, fix25_check_v43/, fix26_latex_v43/ (path roots for the existence check of new paths).
# Not copied (written by hand for v44): b43_common.py, build_v43.py, check_macros_v43.py.
# Usage: PYTHONHASHSEED=0 nice -n 19 python3 -B copy_tools_v44.py
import os, re, json
B = '/home/claude/s1spec/writeup/lemmaU/bulk/'
SRC = B + 'fix24_texts_v43/tools/'; DST = B + 'fix27_texts_v44/tools/'
PY = ['gen_records_v43.py', 'gen_numbers_A_v43.py', 'gen_ei_constants_v43.py', 'check_balance_v43.py', 'check_numbers_v43.py']
for f in PY:
    out = []
    for i, line in enumerate(open(SRC + f).read().split('\n')):
        if not line.lstrip().startswith('#'):
            line = line.replace('fix24_texts_v43', 'fix27_texts_v44').replace('v43', 'v44')
        out.append(line)
    t = '\n'.join(out)
    if f == 'check_balance_v43.py':
        o = "B + 'fix23_delta_check/', "
        assert t.count(o) == 1
        t = t.replace(o, o + "B + 'fix24_texts_v43/', B + 'fix25_check_v43/', B + 'fix26_latex_v43/', ")
        o = '# v44: + fix22, fix23, snap_s3'; assert t.count(o) == 1
        t = t.replace(o, '# v43: + fix22, fix23, snap_s3; v44: + fix24, fix25, fix26')
    if f == 'gen_records_v43.py':
        # v44 additions: the Windows kit runs (Nov 24), the check of text v43 (fix25) and the LaTeX compile of text v43 (fix26), fig_curves.pdf
        o = "     'step4_proof': S + 'writeup/step4_proof/step4_proof.tex'}"; assert t.count(o) == 1
        t = t.replace(o, o[:-1] + ",\n     # v44: Windows runs of the v4 kits (Nov 24), check fix25 of text v43, LaTeX compile fix26 of text v43, the figure of main_v3\n"
                      "     'kit3_win_notes': S + 'kit3v4_run_win1/RUN_NOTES.md', 'kit4b_win_notes': S + 'kit4bv4_run_win1/RUN_NOTES.md',\n"
                      "     'check_v43': S + 'writeup/lemmaU/bulk/fix25_check_v43/CHECK_v43.md', 'latex_v43': S + 'writeup/lemmaU/bulk/fix26_latex_v43/LATEX_v43.md',\n"
                      "     'latex_v43_pdf_sums': S + 'writeup/lemmaU/bulk/fix26_latex_v43/pdf/SHA256SUMS', 'fig_curves': S + 'writeup/main/fig_curves.pdf'}")
        o = "assert len(m['stages']) == 29"; assert t.count(o) == 1
        t = t.replace(o, "# v44: SUMMARY lines of the two Windows runs (run copies with Windows-only fixes; kit v5 pending)\n"
                      "R['kit3_win_summary'] = summ(F['kit3_win_notes'], r'\\d+ PASS, \\d+ FAIL, \\d+ WARN')\n"
                      "R['kit4b_win_summary'] = summ(F['kit4b_win_notes'], r'TOTAL PASS \\[\\d+ PASS, \\d+ FAIL, \\d+ WARN\\]')\n"
                      "assert R['kit3_win_summary'] == '72 PASS, 6 FAIL, 3 WARN', R['kit3_win_summary']\n"
                      "assert R['kit4b_win_summary'] == 'TOTAL PASS [159 PASS, 0 FAIL, 1 WARN]', R['kit4b_win_summary']\n" + o)
    hdr = ('# %s -- (v44: copy of fix24_texts_v43/tools/%s made by tools/copy_tools_v44.py; on code lines fix24_texts_v43 -> fix27_texts_v44 and '
           'v43 -> v44; comments below are the v43 history, unchanged.)\n') % (f.replace('v43', 'v44'), f)
    open(DST + f.replace('v43', 'v44'), 'w').write(hdr + t)
    print('copied', f)
for f in ['keymap_v43.json', 'whitelist_v43.json', 'whitelist_v43_new.json']:
    t = open(SRC + f).read()
    t2 = re.sub(r'([A-Za-z0-9]+)_v43\.(tex|md)', r'\1_v44.\2', t).replace('step3_v43_main', 'step3_v44_main')
    json.loads(t2)
    open(DST + f.replace('v43', 'v44'), 'w').write(t2)
    print('copied', f, t.count('_v43.'), 'file keys renamed')
