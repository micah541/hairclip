# build_v39.py -- builds the step-3 text set v39 in fix17_texts_v39/ from the v38 set in fix16_texts_v38/ (read-only).
# Work order (text v39): (1) every text cites the re-pinned step-4 records (manifest_v38, driver final_tail_v38.py, REPIN_REPORT, the
# exception records, the known F34 constants_needed_v34.json difference) and the v4 kits with their WSL test runs ONLY (no Windows run
# of a v4 kit exists; listed as open); (2) step-4 wording in step4_handoff_v39.tex (make_step4_v39.py); (3) no 0.1475 quote;
# (4) fixes of v38 retag damage found while doing this (listed in TEXTS_REPORT_v39.md).  All record values: tools/records_v39.json.
# Usage: PYTHONHASHSEED=0 nice -n 19 python3 -B build_v39.py   (after gen_records_v39.py, make_step4_v39.py)
import sys, os, json, re, shutil
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import b39_common as C
from b39_common import Builder, D, D8
from v39_values import V
import make_step4_v39  # writes step4_handoff_v39.tex and STEP4_PATCH_NOTE_v39.md (shares C.INSERTS)
v = V
MAN = r'\fn{step4/p2/out_v38/FINAL_manifest_v38.json} (%d stages, %d paths and %d scripts pinned, sha256 \fn{%s})' % (v['nst'], v['paths'], v['scripts'], v['man'])
MAN8 = r'\fn{FINAL_manifest_v38.json} (%d stages, \fn{%s})' % (v['nst'], v['man8'])
KITS = (r'kits \fn{kit_step3_v4.zip} (sha256 \fn{%s}) and \fn{kit_step4b_v4.zip} (sha256 \fn{%s}), test-run in WSL only '
        r'(\fn{kit3v4_run_wsl5}: 80 PASS, 0 FAIL, 1 WARN; \fn{kit4bv4_run_wsl1}: 159 PASS, 0 FAIL, 1 WARN); no Windows run of either v4 kit yet') % (v['kit3'][:8] + '...', v['kit4b'][:8] + '...')
assert v['kit3sum'] == 'SUMMARY: 80 PASS, 0 FAIL, 1 WARN, 0 SKIP' and '159 PASS, 0 FAIL, 1 WARN' in v['kit4bsum']
HDR = '%% %s -- text v39 (fix17_texts_v39, 20 Nov 2026), built by fix17_texts_v39/tools/build_v39.py from fix16_texts_v38/%s (unchanged):\n%% re-pinned step-4 records (manifest_v38 %s), v4 kits (WSL test runs only), step4_handoff_v39.tex; v38 retag fixes.  Not compiled.'
def hdr(b, src): b.head(HDR % (os.path.basename(b.dst), src, v['man8']))
OPEN_STD = (r'the Lean build at commit \fn{d72b744} and the $n=1$ instance; a re-reading of DKS for (D1) and the uniform bounds used in C-1; '
            r'the authors\' confirmations ((P3$\prime$), Lemma~O arc by arc, the single-point statement); the authors\' Windows runs of the two v4 kits; '
            r'LaTeX compilation; the final independent review of the v39 text set').replace(r'\prime', "'").replace("\\'", "'")
out = []
# ---------------------------------------------------------------- copies with retag only
for name in ['b12_preamble_v38.tex', 'E0_LEMMA_v4c_v38.tex']:
    b = Builder(name); hdr(b, name); b.write(); out.append(b.dst)
# ---------------------------------------------------------------- chart lemma body
b = Builder('chart_lemma_body_v38.tex'); hdr(b, 'chart_lemma_body_v38.tex')
b.rep(r'(fin12f4, $\alpha_1$ row:', r'(verifier \fn{fin13v4c}, $\alpha_1$ row (text v39; the v3.3 computation used the \fn{fin12f4} row, which gives the same printed bounds):')
b.rep(r'Then \tC\ (\fn{p2/out_v33/n1_chart.json}):', r'Then \tC\ (\fn{step4/p2/out_v38/n1_chart.json}, stage 19 of %s, on \fn{fin13v4c}; equal apart from args to the chain-v37 file \fn{n1/out/n1_chart_fin13v4c.json}; the v3.3 file \fn{p2/out_v33/n1_chart.json} on \fn{fin12f4} gives the same printed digits):' % MAN8)
b.write(); out.append(b.dst)
# ---------------------------------------------------------------- bootstrap
b = Builder('bootstrap_v38.tex'); hdr(b, 'bootstrap_v38.tex')
b.rep(r'\title{Step 3, text v38: the bootstrap set and its closure}', r'\title{Step 3, text v39: the bootstrap set and its closure}', new=False)
b.rep(r'\author{}\date{20 November 2026 (text v38). \emph{Not compiled.}}', r'\author{}\date{20 November 2026 (text v39; set unchanged from text v38). \emph{Not compiled.}}', new=False)
b.rep(r'This section is the \emph{only} statement of the set in text v38.', r'This section is the \emph{only} statement of the set (text v38, unchanged in text v39). \emph{Records (text v39):} the certificates named in the table below are pinned by the step-4 manifest %s: \fn{C1_ke} stage 22, \fn{GS_tau1} stage 23, \fn{GS_boot_v4c} stage 24, \fn{G1_check} stage 25, \fn{U1_sym} stage 26, flux v37 stage 21, \fn{f12_sliver} stage 20, \fn{runs/P} stage 28; the v38 copies of \fn{C1_ke}, \fn{GS_tau1}, \fn{GS_boot_v4c}, \fn{U1_sym} equal the references cited here apart from args (stage 27); \fn{G1_check} is rerun on chain v37 (stage 25).' % MAN)
b.write(); out.append(b.dst)
# ---------------------------------------------------------------- C-1
b = Builder('C1_identification_v38.tex'); hdr(b, 'C1_identification_v38.tex')
b.rep(r'\date{20 November 2026 (text v38; first version fix round 14). \emph{Not compiled.}}', r'\date{20 November 2026 (text v39; text v38 added (P2-a); first version fix round 14). \emph{Not compiled.}}', new=False)
b.rep(r'(\fn{sup_ke} $=1.7428\cdot10^{-3}$ at slab 180,', r'(\fn{sup_ke} $\le%s$ at slab 180 (text v38 printed $=1.7428\cdot10^{-3}$, the 5-digit truncation; pinned as stage 22 of %s, \fn{out_v38/iface/C1_ke.json}, equal to this file apart from args),' % (v['ke'].replace('0.0017429', r'1.7429\cdot10^{-3}'), MAN8))
b.rep(r'\item \textbf{P2 as an a priori theorem}:', r'\item \textbf{Step-4 wording} (text v39): the hand-off list and the replaced \textup{[A]} sentences are in \fn{step4_handoff_v39.tex} item (3); \fn{step4_proof.tex} (closed) is not edited, the proposed edits are \fn{STEP4_PATCH_NOTE_v39.md}.' + '\n' + r'\item \textbf{P2 as an a priori theorem}:')
b.write(); out.append(b.dst)
# ---------------------------------------------------------------- GS lemma
b = Builder('GS_lemma_v38.tex'); hdr(b, 'GS_lemma_v38.tex')
b.rep(r'\date{20 November 2026 (text v38; first version fix round 14). \emph{Not compiled.}}', r'\date{20 November 2026 (text v39; first version fix round 14). \emph{Not compiled.}}', new=False)
b.rep(r"""\paragraph{Review of \fn{fix17_texts_v39/continuation_v39.tex} l.41 (``Strictness (sketch added in text v38)'').}""",
      r"""\paragraph{Review of \fn{fix12_texts/continuation_v37.tex} l.41 (``Strictness (sketch added in text v37)'').} \emph{(Text v39 restores this reference; text v38 had retagged it to its own file by mistake.)}""")
b.rep(r'\paragraph{Consequential edits owed in fix round 14 (done in text v38, except the step-4 wording).}', r'\paragraph{Consequential edits owed in fix round 14 (done in text v38; the step-4 wording is in \fn{step4_handoff_v39.tex}, text v39).}')
b.rep(r"""Shipping \fn{GS_boot_v4c} as a pinned stage (with $\sqrt{L_{\max}}(1+2C'/L)^{1/2}$, GS-6, and the aligned cell, GS-2) is part of the re-pin.""",
      r"""\fn{GS_boot_v4c} is pinned as stage 24 of %s (\fn{out_v38/iface/GS_boot_v4c.json}, equal to the fix-15 output apart from args), and \fn{GS_tau1} as stage 23. The pinned script is the fix-15 script unchanged: the factor $(1+2C'/L)^{1/2}$ on $\sqrt{L_{\max}}$ (GS-6) and the aligned cell (GS-2) are \emph{not} applied; review fix15 rates both negligible (margins $\times24$ and more). They remain open records items.""" % MAN)
b.write(); out.append(b.dst)
# ---------------------------------------------------------------- U-1
b = Builder('U1_regularity_v38.tex'); hdr(b, 'U1_regularity_v38.tex')
b.rep(r'\date{20 November 2026 (text v38; first version fix round 10', r'\date{20 November 2026 (text v39; first version fix round 10', new=False)
b.rep(r'(review fix15 D-5; text v38 had $\sqrt\Lambda', r'(review fix15 D-5; text v37 had $\sqrt\Lambda', new=False)
b.rep(r'Step~4\'s hand-off list should name GS next to (H1) (step-4 wording, open).'.replace("\\'", "'"),
      r"Step~4's hand-off list with GS next to (H1): \fn{step4_handoff_v39.tex} item (4) (text v39; \fn{step4_proof.tex} is not edited, see \fn{STEP4_PATCH_NOTE_v39.md}); \fn{GS_tau1}, \fn{GS_boot_v4c} and \fn{U1_sym} are stages 23, 24, 26 of %s." % MAN8)
b.write(); out.append(b.dst)
# ---------------------------------------------------------------- G1 export
b = Builder('G1_export_k2_v38.tex'); hdr(b, 'G1_export_k2_v38.tex')
b.rep(r're-checked by \fn{G1_check.py} in this directory (output \fn{G1_check.json}, \fn{all_pass: true});',
      r're-checked by \fn{G1_check.py} (\fn{fix10_proofs}, v36 data) and, on the chain-v37 data of this text, by its v38 copy \fn{step4/p2/v38/G1_check_v38.py}, stage 25 of %s (\fn{out_v38/iface/G1_check.json}, \fn{all_pass: true}; text v39: the rerun owed in text v38 is done);' % MAN)
b.write(); out.append(b.dst)
# ---------------------------------------------------------------- ei lemma
b = Builder('ei_lemma_v4c.tex'); hdr(b, 'ei_lemma_v4c.tex')
b.rep(r'\author{}\date{20 November 2026 (text v38). \emph{Not compiled.}}', r'\author{}\date{20 November 2026 (text v39). \emph{Not compiled.}}', new=False)
b.rep(r'(the text-v38 copy of', r'(the text-v39 copy of', new=False)
b.rep('\x00K000\x00', 'the E_I functions are unchanged from fin13v4b, rhs_iv.py b045b9ff..., only term H8 was added [text v39 restores this comment: the v38 builder left a NUL-delimited placeholder here]', new=False)
b.rep(r'text v38 printed $\times1.679$, which omitted A-1;', r'text v37 printed $\times1.679$, which omitted A-1;', new=False)
b.rep(r'\item Audited since text v38', r'\item Audited since text v37', new=False)
b.write(); out.append(b.dst)
# ---------------------------------------------------------------- flux lemma
b = Builder('flux2_lemma_v38.tex'); hdr(b, 'flux2_lemma_v38.tex')
b.rep(r'\date{20 November 2026 (text v38). \emph{Not compiled.}', r'\date{20 November 2026 (text v39). \emph{Not compiled.}', new=False)
b.rep(r'\noindent\textbf{Status (text v38).}', r'\noindent\textbf{Status (text v39; v38 status follows).}', new=False)
b.rep(r'Open: the outputs \fn{fix13_v37chain/flux/out/*} are in no manifest and no kit (re-pin; first regenerate \fn{SHA256SUMS_v37} and pin a non-log file instead of \fn{V37_LOG.md}, review B-v37-1);',
      r'Records (text v39): the outputs \fn{fix13_v37chain/flux/out/*} and \fn{flux/texts} are pinned by stage 21 of %s, which reruns the 8 flux stages into \fn{out_v38/flux} and reproduces every file (JSON apart from args; pkl and texts byte for byte; verdict PASS, combined margin $\times%s$), and they are in \fn{kit_step3_v4.zip} (WSL test run only). Open: review item B-v37-1 (\fn{inputs.json} pins an older sha256 of \fn{V37_LOG.md}, the one non-args difference of the rerun; \fn{SHA256SUMS_v37} not regenerated);' % (MAN, v['flux']))
b.write(); out.append(b.dst)
# ---------------------------------------------------------------- lemma S(A)
b = Builder('lemmaS_A_v38.tex'); hdr(b, 'lemmaS_A_v38.tex')
b.rep(r'\title{Lemma S(A), repaired (step 3, text v38):', r'\title{Lemma S(A), repaired (step 3, text v39):', new=False)
b.rep(r'\date{20 November 2026 (text v38: chain v37, verifier \fn{fin13v4c}).', r'\date{20 November 2026 (text v39: chain v37, verifier \fn{fin13v4c}).', new=False)
b.rep(r'(generated by \fn{fix16_texts_v38/tools/gen_numbers_A_v38.py}; v37 values,', r'(generated by \fn{fix17_texts_v39/tools/gen_numbers_A_v39.py}; v37 values,', new=False)
b.rep(r'these outputs carry no args/sha256 block, fix11 B1, to be pinned at the re-pin)', r'these outputs carry no args/sha256 block, fix11 B1; text v39: the \fn{fix13_v37chain/t2_A/base_ibp} outputs are now pinned by sha256 as inputs of \fn{FINAL_manifest_v38.json} (stage 20), but they still carry no args block, open)')
b.rep(r'Both are certified for $f-U_1$ in \fn{step4/p2/out_v33/FINAL_manifest_v33.json} (21 stages, sha256 \fn{8bc47e83...}):',
      r'Both are certified for $f-U_1$ in \fn{step4/p2/out_v33/FINAL_manifest_v33.json} (21 stages, sha256 \fn{8bc47e83...}), carried by hash into %s as stages 0--4 and 10--17:' % MAN)
b.rep(r'by the step-4 v33 stage \fn{n1_chart} \tC.', r'by the step-4 v33 stage \fn{n1_chart}, rerun on \fn{fin13v4c} as stage 19 of \fn{FINAL_manifest_v38.json} \tC.')
b.rep(r'\item \textbf{v38 open items (not mathematics of this file):}', r'\item \textbf{Open items (text v39; not mathematics of this file):}', new=False)
b.rep(r'the step-4 re-pin of stages 8, 9, 18, 19, 20 of manifest\_v33 on v37 is not done (\fn{fix13_v37chain/V37_REPORT.md}, re-pin list).',
      r'the step-4 re-pin of stages 8, 9, 18, 19, 20 on v37 is done: %s (\fn{step4/p2/out_v38/REPIN_REPORT.md}).' % MAN)
b.write(); out.append(b.dst)
# ---------------------------------------------------------------- continuation
b = Builder('continuation_v38.tex'); hdr(b, 'continuation_v38.tex')
b.rep(r'\title{Step 3, text v38: continuation,', r'\title{Step 3, text v39: continuation,', new=False)
b.rep(r'\author{}\date{20 November 2026 (text v38). \emph{Not compiled.}}', r'\author{}\date{20 November 2026 (text v39). \emph{Not compiled.}}', new=False)
b.rep(r'text v38 printed $7.12\cdot10^{-3}$', r'text v37 printed $7.12\cdot10^{-3}$', new=False)
b.rep(r'(The sketch of text v38 is withdrawn:', r'(The sketch of text v37 is withdrawn:', new=False)
b.rep(r'(text v38 printed the difference lhs$-$rhs, $177.07$, as a margin)', r'(text v37 printed the difference lhs$-$rhs, $177.07$, as a margin)', new=False)
b.rep(r'(\fn{step4_proof.tex} l.672, [A]) is no longer used.', r'(\fn{step4_proof.tex} l.672, [A]) is no longer used; the step-4 wording is \fn{step4_handoff_v39.tex} item (3) (text v39).')
b.rep(r"step~4's hand-off list (Theorem~4(c)) should name GS next to (H1) (step-4 wording, open). The strictness step is certified (\fn{GS_boot_v4c.json}); shipping \fn{GS_boot_v4c} as a pinned stage is part of the re-pin.",
      r"step~4's hand-off list with GS next to (H1) is \fn{step4_handoff_v39.tex} item (4) (text v39; \fn{step4_proof.tex} is not edited, proposed edits in \fn{STEP4_PATCH_NOTE_v39.md}). The strictness step is certified (\fn{GS_boot_v4c.json}). Both are pinned: stages 23 and 24 of %s." % MAN8)
b.rep(r'Step 4 stage \fn{h4far_v33} (\fn{FINAL\_manifest\_v33.json}, 21 stages, \fn{8bc47e83...}; \fn{step4/p2/v33/M9_REPORT.md})',
      r'Step 4 stage \fn{h4far_v33} (\fn{FINAL\_manifest\_v33.json}, 21 stages, \fn{8bc47e83...}; \fn{step4/p2/v33/M9_REPORT.md}), re-pinned on zone~III of chain v37 as stage 18 of %s (same rows, same margin),' % MAN8)
b.write(); out.append(b.dst)
# ---------------------------------------------------------------- lemma U
b = Builder('lemmaU_v38.tex'); hdr(b, 'lemmaU_v38.tex')
b.rep(r'\large proof text, revision v38 (numbers and hashes of the (A) chain v37 on verifier fin13v4c)}', r'\large proof text, revision v39 (numbers and hashes of the (A) chain v37 on verifier fin13v4c; step-4 records of the re-pin v38)}', new=False)
b.rep(r'\date{20 November 2026 (v38; revision of \fn{fix17_texts_v39/lemmaU_v39.tex}, itself a revision of',
      r'\date{20 November 2026 (v39; revision of \fn{fix16_texts_v38/lemmaU_v38.tex}, a revision of \fn{fix12_texts/lemmaU_v37.tex}, itself a revision of', new=False)
b.rep(r"""the step-4 re-pin on \fn{fin13v4c} and chain v37 (stages 8, 9, 18, 19, 20 of manifest\_v33, \fn{fix13_v37chain/V37_REPORT.md} ``re-pin list'') is pending and needs Micah\'s go; the final independent review of the v38 text set is pending.""".replace("\\'", "'"),
      r"the step-4 re-pin on \fn{fin13v4c} and chain v37 is done (%s; \fn{step4/p2/out_v38/REPIN_REPORT.md}), and the %s; the final independent review of the v39 text set is pending." % (MAN, KITS))
b.rep(r'final after reviews M9, N1/N2 and F1/F2; unchanged and \emph{not} re-pinned on \fn{fin13v4c}/chain v37: the re-pin of stages 8, 9, 18, 19, 20 is pending, \fn{fix13_v37chain/V37_REPORT.md} re-pin list, needs Micah\'s go);'.replace("\\'", "'"),
      r'final after reviews M9, N1/N2 and F1/F2; unchanged), superseded for step~3 by the re-pin %s: 16 stages carried from manifest\_v33 by hash, stages 8, 9, 18, 19, 20 re-pinned on \fn{fin13v4c} and chain v37, stages 21--28 new (flux v37, C1\_ke, GS\_tau1, GS\_boot, G1\_check, U1\_sym, interface comparison, \fn{runs/P} pins); driver \fn{step4/p2/v38/final_tail_v38.py}, report \fn{step4/p2/out_v38/REPIN_REPORT.md};' % MAN)
b.rep(r'changes listed in \fn{fix16_texts_v38/LEMMAU_MAIN_v37_NOTES.md}', r'changes listed in \fn{fix12_texts/LEMMAU_MAIN_v37_NOTES.md} (v37), \fn{fix16_texts_v38/TEXTS_REPORT_v38.md} (v38) and \fn{fix17_texts_v39/TEXTS_REPORT_v39.md} (v39)', new=False)
b.rep(r'proved there source by source (Theorem ``Lemma D, $k=2$\'\'; sha256 \fn{4004a774f5ef41764c9feb32327b8d42c972f70340ec4cdca9cf79b98cf45e05})'.replace("\\'", "'"),
      r'proved there source by source (Theorem ``Lemma D, $k=2$\'\'; sha256 of its first version \fn{fix10_proofs/G1_export_k2.tex}: \fn{4004a774f5ef41764c9feb32327b8d42c972f70340ec4cdca9cf79b98cf45e05}; text v38 paired this hash with its own copy by mistake)'.replace("\\'", "'"))
b.rep(r'(\fn{fix17_texts_v39/U1_regularity_v39.tex}, sha256 \fn{827a5f51ed7af4d7a21c663183b6cdb06a63fffb64c420fa795d49f78ff2e9e7};',
      r'(\fn{fix17_texts_v39/U1_regularity_v39.tex}; sha256 of its first version \fn{fix10_proofs/U1_regularity.tex}: \fn{827a5f51ed7af4d7a21c663183b6cdb06a63fffb64c420fa795d49f78ff2e9e7};')
b.rep(r'\fn{tools/numbers_check_v38.json}', r'\fn{tools/numbers_check_v39.json}', new=False)
b.rep(r'it replaces stage 19 at the re-pin).', r'it replaces stage 19 in manifest\_v38: \fn{out_v38/n1_chart.json}, sha256 \fn{%s}, equal to it apart from args).' % (v['n1v38'][:8] + '...'))
b.rep(r'Nothing in this item is open for the v3.3 chain as pinned; the re-pin of stages 8, 9, 18, 19, 20 of manifest\_v33 on \fn{fin13v4c} and chain v37 is pending (\fn{fix13_v37chain/V37_REPORT.md}, re-pin list), and the hand-off values (H1)--(H6) are unchanged.',
      r'Nothing in this item is open for the chain as pinned: the re-pin on \fn{fin13v4c} and chain v37 is %s (stage 18 reruns the far-family H4 on zone~III of chain v37: minimum margin $\times%s$, unchanged), and the hand-off values (H1)--(H6) are unchanged.' % (MAN8, v['h4far']))
b.rep(r'(re-pin and kits, the Lean build and the instance at $n=1$, DKS, records, LaTeX, the final review;',
      r'(the authors\' Windows runs of the v4 kits, the Lean build and the instance at $n=1$, DKS, records, LaTeX, the final review; the re-pin is done, manifest\_v38;'.replace("\\'", "'"))
b.write(); out.append(b.dst)
# ---------------------------------------------------------------- step-3 main
b = Builder('step3_v38_main.tex'); hdr(b, 'step3_v38_main.tex')
b.rep(r'\title{Step 3 of the hair-clip proof, text v38 (', r'\title{Step 3 of the hair-clip proof, text v39 (', new=False)
b.rep(r'\date{20 November 2026 (text v38). \emph{Not compiled.} v3.1, v3.2, the v3.3 texts and text v38 are unchanged; this is a new file.}', r'\date{20 November 2026 (text v39). \emph{Not compiled.} v3.1, v3.2, the v3.3 texts, text v38 and \fn{step4_proof.tex} are unchanged; this is a new file.}', new=False)
b.rep(r'\textbf{Status (text v38). Step~3 is NOT closed.}', r'\textbf{Status (text v39). Step~3 is NOT closed.} \emph{New in v39:} the step-4 records are the re-pin %s (driver \fn{final_tail_v38.py}); the step-4 hand-off as step~3 reads it, with (H6) by C-1 and (GS), is \fn{step4_handoff_v39.tex}; %s. \emph{Text v38 status:}' % (MAN, KITS))
b.rep(r'\textbf{Open:} the step-4 re-pin on \fn{fin13v4c}/v37 and the kits (the authors\' go); the Lean build and the $n=1$ instance;'.replace("\\'", "'"),
      r"\textbf{Open:} the authors' Windows runs of the v4 kits; the Lean build and the $n=1$ instance;")
b.rep(r"(Text v38: for the clip, (H6) is Theorem C-1 of \fn{C1_identification_v39.tex}, not step~4's [A] sentence.)",
      r"(For the clip, (H6) is Theorem C-1 of \fn{C1_identification_v39.tex}, not step~4's [A] sentence; step-4 wording: \fn{step4_handoff_v39.tex} item (3).)")
b.rep(r'\fn{fix14_interface/GS_tau1.json}; Lemma GS-1 of \fn{GS_lemma_v39.tex}).', r'\fn{fix14_interface/GS_tau1.json}, pinned as stage 23 of %s; Lemma GS-1 of \fn{GS_lemma_v39.tex}; hand-off item (4) of \fn{step4_handoff_v39.tex}).' % MAN8)
b.rep(r'\noindent\emph{Step-4 values} \tC\ (\fn{step4/p2/out/FINAL_manifest.json}, 34 stages, and \fn{step4/p2/out_v33/FINAL_manifest_v33.json}, 21 stages, 366 paths, sha256 \fn{8bc47e83e370fec3e856e9eb61cf0863ec0738320806b294cb2d74772db983ec};',
      r'\noindent\emph{Step-4 values} \tC\ (records of text v39: %s, driver \fn{step4/p2/v38/final_tail_v38.py}, sha256 \fn{%s}, report \fn{step4/p2/out_v38/REPIN_REPORT.md}; it carries the v3.3 stages 0--7, 10--17 by hash, re-pins 8, 9, 18, 19, 20 on \fn{fin13v4c}/chain v37 and adds 21--28; it verifies the unchanged \fn{step4/p2/out/FINAL_manifest.json}, 34 stages, with the one documented exception \fn{EXCEPTION_v32_B8_check_outA5.json}, and the unchanged \fn{step4/p2/out_v33/FINAL_manifest_v33.json}, 21 stages, 366 paths, sha256 \fn{8bc47e83e370fec3e856e9eb61cf0863ec0738320806b294cb2d74772db983ec}, with the one known F34 difference of \fn{constants_needed_v34.json}; details in \fn{step4_handoff_v39.tex} \S1;' % (MAN, v['drv']))
b.rep(r'the stage \fn{n1_chart} (stage 19) and the explicit per-stage environment are N1 and N2 (\fn{fix2/n12/N12_REPORT.md}).',
      r'the stage \fn{n1_chart} (stage 19; rerun on \fn{fin13v4c} in manifest\_v38) and the explicit per-stage environment are N1 and N2 (\fn{fix2/n12/N12_REPORT.md}).')
b.rep(r'(v37 \key{supersolution.M1z_hyp_check}; text v38 printed $0.28578$', r'(v37 \key{supersolution.M1z_hyp_check}; text v37 printed $0.28578$', new=False)
b.rep(r'\item \textbf{Hypotheses of the theorem} \tH: (H1)--(H5), (H2f) (step 4, values unchanged; step~4 not re-pinned on \fn{fin13v4c}/v37); (H6) and (GS) for the clip are now C-1 and GS-1; the hand-off list of step~4 should name GS and cite C-1 for (H6) (step-4 wording).',
      r'\item \textbf{Hypotheses of the theorem} \tH: (H1)--(H5), (H2f) (step 4, values unchanged; re-pinned on \fn{fin13v4c}/v37 in manifest\_v38); (H6) and (GS) for the clip are C-1 and GS-1; the step-4 wording (hand-off list with GS, (H6) by C-1) is \fn{step4_handoff_v39.tex}, and the corresponding edits of the closed \fn{step4_proof.tex} are proposed in \fn{STEP4_PATCH_NOTE_v39.md} (not applied).')
old = b.between(r'\item \textbf{Re-pin and kits} (the authors\' go):'.replace("\\'", "'"), r'\item \textbf{Records}')
b.rep(old, r'\item \textbf{Re-pin and kits}: \emph{re-pin done} (%s; driver \fn{final_tail_v38.py} consults \fn{EXCEPTION_v32_B8_check_outA5.json}; \fn{runs/P_A}, \fn{runs/reg_fin13v4b} through \fn{EXCEPTION_v4c_PA_reg_e2v4b.json}); %s. \textbf{Open:} the authors\' independent Windows runs of both v4 kits.'.replace("\\'", "'") % (MAN8, KITS) + '\n')
b.rep(r'the pins of \fn{runs/P_A} and \fn{runs/reg_fin13v4c}, the checksum files,', r'the pins of \fn{runs/P_A} and \fn{runs/reg_fin13v4b} (now covered by the exception record; a rerun into new directories is still recommended), the checksum files, the GS\_boot items GS-2, GS-6,')
b.rep(r'\item \textbf{LaTeX compile} of every v38 file; \textbf{final independent review} of this text set.', r'\item \textbf{LaTeX compile} of every v39 file; \textbf{final independent review} of this text set.', new=False)
b.write(); out.append(b.dst)
# ---------------------------------------------------------------- main_v3
b = Builder('main_v3_v38.tex'); hdr(b, 'main_v3_v38.tex')
b.rep(r'Draft v3, revision v38 of November 20, 2026', r'Draft v3, revision v39 of November 20, 2026', new=False)
b.rep(r'\textbf{Status (20 November 2026, v38): not compiled; step 3 is NOT closed.}', r'\textbf{Status (20 November 2026, v39): not compiled; step 3 is NOT closed.} \emph{New in v39:} the step-4 re-pin on \fn{fin13v4c}/v37 is done (%s; driver \fn{final_tail_v38.py}; exception records for the B8 file and for \fn{runs/P_A}, \fn{runs/reg_fin13v4b}); the hand-off as step~3 reads it, with (H6) by C-1 and the new item (GS), is \fn{step4_handoff_v39.tex}; %s. \emph{v38 text of this box:}' % (MAN, KITS))
b.rep(r'\textbf{Open:} the step-4 re-pin on \fn{fin13v4c}/v37 (stages 8, 9, 18, 19, 20 of manifest\_v33, a new final\_tail driver reading the EXCEPTION record, a new manifest) and the kits, which need the authors\' go; the Lean build and the $n=1$ instance; a re-reading of DKS for (D1); the \texttt{[read]} inputs and \texttt{[hand]} steps listed in the step-3 texts; LaTeX compilation; the final independent review of the v38 text set.'.replace("\\'", "'"),
      r'\textbf{Open:} %s; the \texttt{[read]} inputs and \texttt{[hand]} steps listed in the step-3 texts.' % OPEN_STD)
b.rep(r'text v38 printed the differences $177.07$, $175.09$ as margins', r'text v37 printed the differences $177.07$, $175.09$ as margins', new=False)
b.rep(r'Proof text: \dref{LU} (\fn{fix17_texts_v39/lemmaU_v39.tex}, revision v38 of', r'Proof text: \dref{LU} (\fn{fix17_texts_v39/lemmaU_v39.tex}, revision v39 of', new=False)
b.rep(r'for v3.3 the addendum \dref{R4v33} (new stages in \fn{step4/p2/v33}, manifest \fn{FINAL_manifest_v33.json}). The step-4 document itself is closed and unchanged.',
      r'for v3.3 the addendum \dref{R4v33} (new stages in \fn{step4/p2/v33}, manifest \fn{FINAL_manifest_v33.json}) and the re-pin \dref{R4v38} (%s, driver \fn{step4/p2/v38/final_tail_v38.py}); the hand-off as step~3 v39 reads it is \dref{S4H}. The step-4 document itself is closed and unchanged.' % MAN)
b.rep(r'satisfies (H1)--(H6) of Theorem~\ref{thm:W} (Table~\ref{tab:handoff}), uniformly for every good time in the slab.',
      r'satisfies (H1)--(H6) of Theorem~\ref{thm:W} (Table~\ref{tab:handoff}), uniformly for every good time in the slab; for step~3 v3.3, (H6) is Theorem C-1 and the hand-off list includes (GS) \dref{S4H}.')
b.rep(r'(a \fn{fin12f4} value: stage \fn{f12\_sliver} still reads the \fn{fin12f4} verifier and must be re-pinned on \fn{fin13v4c} together with \fn{v4bcert/out/strip_cubic_v4c.json}, \fn{fix13_v37chain/V37_REPORT.md} re-pin list; the \fn{fin13v4c} value $s_{I2}\le\sItwoV$ is larger; TODO, needs the authors\' go)'.replace("\\'", "'"),
      r"""(a \fn{fin12f4} value of manifest\_v33; re-pinned in \dref{R4v38} as stage 20 on \fn{fin13v4c}, with the rate-table pins incl.\ \fn{v4bcert/out/strip_cubic_v4c.json}: $sI_2'\le%s$, $M_1'\le%s\le%s$, %d of %d checks pass)""" % (v['sI2'], v['M1p'], v['M1A'], v['f12p'], v['f12n']))
b.rep(r'All passes \textup{[cert]}.' + '\n' + r'\end{itemize}',
      (r'All passes \textup{[cert]}.' + '\n' + r'\item \emph{Re-pin v38 (new, \dref{R4v38}).} %s, written by \fn{step4/p2/v38/final_tail_v38.py}: the 16 v3.3 stages 0--7, 10--17 carried by hash; stages 8 (H7 ratio, margin $\times%s$), 9 (covering notes, H6 task-5 margin $\times%s$), 18 (far H4, $\times%s$), 19 (N1 chart on \fn{fin13v4c}) and 20 (\fn{f12\_sliver}) re-pinned on the verifier \fn{fin13v4c} and the (A) chain v37; stages 21--28 new: flux v37 rerun and pin (combined margin $\times%s$), the interface certificates C1\_ke, GS\_tau1, GS\_boot, G1\_check, U1\_sym and their comparison, and the pins of \fn{step3nl_v4c/runs/P}. The closed v3.2 manifest verifies with the one documented exception \fn{EXCEPTION_v32_B8_check_outA5.json}; manifest\_v33 verifies with the one known F34 difference of \fn{constants_needed_v34.json}; \fn{runs/P_A}, \fn{runs/reg_fin13v4b} are accepted only through \fn{EXCEPTION_v4c_PA_reg_e2v4b.json}. Two passes are byte-identical apart from timing \textup{[cert]}.' + '\n' + r'\end{itemize}') % (MAN, v['h7mar'], v['h6mar'], v['h4far'], v['flux']))
b.rep(r'(H5) fixed-$L$ cone H7 & $t^*\ge$ ratio & ratio $\le1.6000014\cdot10^{-7}$, $t^*\ge3.68334\cdot10^{-7}$ & $\times2.30$\\\bottomrule',
      (r'(H5) fixed-$L$ cone H7 & $t^*\ge$ ratio & ratio $\le1.6000014\cdot10^{-7}$, $t^*\ge3.68334\cdot10^{-7}$ & $\times2.30$\\' + '\n' +
       r'(H6) regularity (v3.3) & smooth immersion & Theorem C-1; $\sup|\tilde\kappa e|\le%s$ (\dref{R4v38} stage 22) & $\times%s$ (on $\frac12$)\\' + '\n' +
       r"""(GS) tip chart (v3.3) & $\sup|ku|<\frac12$, $\sup|u'|<1$ & $%s$, $%s$ on $|x|\le%s$ (\dref{R4v38} stage 23) & $\times%s$, $\times%s$\\\bottomrule""") % (v['ke'], v['keMar'], v['gsKu'], v['gsUp'], v['gsXS'], v['gsKuMar'], v['gsUpMar']))
old = b.between(r'Reproduction kits: \fn{kit_step3.zip}', r'\subsection{Independent reproduction}')
b.rep(old, r'Reproduction kits: \fn{kit_step3.zip} (step~3 as of v3.1), \fn{kit_step4a.zip} ($\Gamma_0$ and P1; a revised checker is in \fn{kit_step4a_v3.zip}), and for the v3.3 line on \fn{fin13v4c}/chain v37 the two v4 kits: \fn{kit_step3_v4.zip} (sha256 \fn{%s}; verifier \fn{fin13v4c}, the certificates of \fn{v4bcert}, the N1 chart, the (A) chain, R3, Lemma~U and zone~III, flux v37, G1, U1, GS\_boot) and \fn{kit_step4b_v4.zip} (sha256 \fn{%s}; base, P2, hand-off, the v3.3 stages and the 29-stage v38 driver). Both write only inside the kit (a guard blocks writes, \texttt{chdir} and child processes outside it).' % (v['kit3'], v['kit4b']) + '\n\n')
b.rep(r'\item \emph{Kit B:} pending.', r'\item \emph{Kit B and the v4 kits:} not yet run on the authors\' machine. Test runs in WSL on the machine that holds the proof tree (so not independent): \fn{kit_step3_v4.zip} in \fn{kit3v4_run_wsl5}, 80 PASS, 0 FAIL, 1 WARN (float noise in an unchanged kit-v2 output); \fn{kit_step4b_v4.zip} in \fn{kit4bv4_run_wsl1}, 159 PASS, 0 FAIL, 1 WARN (two interface comparisons recorded, not asserted, by design), 6~h~20~min. The authors\' Windows runs of both v4 kits are open.'.replace("\\'", "'"))
b.rep(r'Open: (a) the step-4 re-pin on \fn{fin13v4c}/v37 (stages 8, 9, 18, 19, 20 of manifest\_v33, the flux outputs, \fn{C1_ke.json}, \fn{GS_tau1.json}, \fn{GS_boot_v4c.json}, a new final\_tail driver reading the EXCEPTION record, a new manifest) and the kits;',
      r'Done since v38: (a) the step-4 re-pin (\dref{R4v38}, %s) and the two v4 kits with WSL test runs. Open: (a$\prime$) the authors\' Windows runs of both v4 kits;'.replace(r'\prime', "'").replace("\\'", "'") % MAN8)
b.rep(r'\item \textbf{Status v38 (20 November 2026): step 3 is not closed.}', r'\item \textbf{Status v39 (20 November 2026): step 3 is not closed.}', new=False)
b.rep(r'\pth{fix16_texts_v38/TEXTS_REPORT_v38.md};', r'\pth{fix16_texts_v38/TEXTS_REPORT_v38.md}; \pth{fix17_texts_v39/TEXTS_REPORT_v39.md};', new=False)
b.rep(r'(g) the final independent review of the v38 text set.', r'(g) the final independent review of the v39 text set; (h) the authors\' confirmations ((P3$\prime$), Lemma~O arc by arc, the single-point statement).'.replace(r'\prime', "'").replace("\\'", "'"))
b.rep(r'\item \textbf{Reproduction:} kit~B (P2 and the hand-off) has not been run end to end; none of the kits contains the v3.3 additions (\fn{KITS_TODO.md}). The authors decide whether to build and run them.',
      r'\item \textbf{Reproduction:} the v4 kits \fn{kit_step3_v4.zip} and \fn{kit_step4b_v4.zip} contain the v3.3 line on \fn{fin13v4c}/chain v37 and pass their WSL test runs (\S\ref{ss:repro}); the authors\' independent Windows runs are open.'.replace("\\'", "'"))
b.rep(r'\dref{R4v33} & step-4 addendum for v3.3 &', r'\dref{R4v38} & step-4 re-pin v38 (29 stages) & \pth{s1spec/step4/p2/out_v38/FINAL_manifest_v38.json}; \pth{s1spec/step4/p2/out_v38/REPIN_REPORT.md}; \pth{s1spec/step4/p2/v38/final_tail_v38.py}\\' + '\n' + r'\dref{S4H} & step-4 hand-off as step 3 v39 reads it & \pth{s1spec/writeup/lemmaU/bulk/fix17_texts_v39/step4_handoff_v39.tex}; \pth{STEP4_PATCH_NOTE_v39.md}\\' + '\n' + r'\dref{R4v33} & step-4 addendum for v3.3 &')
b.rep(r'\dref{LU} & Lemma U (v3.3), proof text (v38) and certificates (chain v37) &', r'\dref{LU} & Lemma U (v3.3), proof text (v39) and certificates (chain v37) &', new=False)
b.write(); out.append(b.dst)
# ---------------------------------------------------------------- erratum (historical: no retag of the body) and ERRATUM md
part4 = open(D + 'tools/erratum_part4_v39.md').read()
b = Builder('erratum_v38.tex', retag_names=False)
b.head('%% erratum_v39.tex -- text v39: erratum_v38.tex (fix16_texts_v38, unchanged) verbatim as Parts 1-3, plus Part 4 (section O, text v39).  Built by fix17_texts_v39/tools/build_v39.py.  Not compiled.')
b.rep(r'\title{Erratum v38 to step~3', r'\title{Erratum v39 to step~3', new=False)
b.rep(r'\date{v38, November 20, 2026 (first version November 13, 2026)}', r'\date{v39, November 20, 2026 (first version November 13, 2026)}', new=False)
b.rep(r'\noindent\textbf{v38 (November 20, 2026;', r'\noindent\textbf{v39 (November 20, 2026; \fn{writeup/lemmaU/bulk/fix17_texts_v39/}, built by \fn{tools/build_v39.py}).} Parts~1--3 are \fn{fix16_texts_v38/erratum_v38.tex} verbatim; Part~4 (section O, at the end) is new: the records of the step-4 re-pin v38, the step-4 wording, and corrections of text v38.' + '\n\n' + r'\noindent\textbf{v38 (November 20, 2026;', new=False)
tex4 = '\n'.join(l for l in part4.split('\n') if not l.startswith('#md'))
tex4 = tex4.replace('@MAN@', MAN).replace('@KITS@', KITS).replace('@OPEN@', OPEN_STD)
b.rep('\\end{document}', tex4 + '\n\\end{document}')
b.write(); out.append(b.dst)
md = open(D8 + 'ERRATUM_v38.md').read()
open(D + 'ERRATUM_v39.md', 'w').write('# ERRATUM_v39 (text v39): ERRATUM_v38.md verbatim, then Part 4\n\nPart 4 (text v39) is section O of `erratum_v39.tex`; its source is `tools/erratum_part4_v39.md`.\n\n---\n\n' + md)
out.append(D + 'ERRATUM_v39.md')
# ---------------------------------------------------------------- records errata md (copy + section)
rec = open(D8 + 'RECORDS_ERRATA_v38.md').read()
rec9 = open(D + 'tools/records_errata_add_v39.md').read()
open(D + 'RECORDS_ERRATA_v39.md', 'w').write('# RECORDS_ERRATA_v39: RECORDS_ERRATA_v38.md verbatim, then the v39 status section at the end\n\n' + rec + '\n\n' + rec9)
C.INSERTS.append({'file': 'RECORDS_ERRATA_v39.md', 'text': rec9})
out.append(D + 'RECORDS_ERRATA_v39.md')
C.save_inserts()
# sanity: no 0.1475 quote outside the erratum history
for f in out:
    t = open(f).read()
    if os.path.basename(f).startswith(('erratum', 'ERRATUM', 'RECORDS')): continue
    assert '0.1475' not in t, f
print('built %d files' % len(out))
