# main_v2.tex: changes from main.tex
Line numbers are in main_v2.tex (old line in brackets). main.tex untouched. Balance check: braces 0, environments matched (no LaTeX here, not compiled).

(O) wording: HYPOTHESIS_O.md did not exist, so a PLACEHOLDER is used (marked text plus a one-sentence informal description).

- L12 [12] Theorem title: "conditional on Hypotheses (F) and (O)".
- L16 [new] `\newtheorem*{hypO}` for Hypothesis (O).
- L56 [55] Date: "Draft v2 of November 13, 2026".
- L60 [59] Abstract: "two stated hypotheses ... (Hypotheses (F) and (O))".
- L64 [63] Status box: date 13 Nov; adds the Lemma S(A) omissions found by the second review, the repair with corrected constants (erratum in preparation), new (O), "two explicit hypotheses".
- L90-93 [new] Hypothesis (O) environment after (F): PLACEHOLDER `[Hypothesis (O): wording pending, see lemmaU/bulk/t1_statements/HYPOTHESIS_O.md]` + informal sentence (outer cutoff of S(A) on [350,700] adds only small ramp terms).
- L95 [90] Main Theorem: "Assume Hypotheses (F) and (O)".
- L110 [105] "closes under (F) and (O)".
- L113-114 [108-109] Remark title "beyond (F) and (O)"; item (a) independent of (F) or (O).
- L219 [214] Dependency figure node: "Hypotheses (F), (O)" (outer cutoff added).
- L270 [265] "If the bootstrap closes (under (F) and (O), ...)".
- L287-290 [282] Subsection retitled "The gap, and Hypotheses (F) and (O)"; new paragraph "Status of Lemma S(A)": omitted ramp terms (top pair, cutoff x jet, U_j, I_3), repaired with corrected constants (erratum in preparation); at v3.1 thresholds E_2 = 0.2686, core M1 x2.73, flux x4.9 / x5.3, all repaired margins >= x2; outer cutoff on [350,700] = (O); (F) unchanged.
- L295 [287] "under (F) and (O) the bootstrap closes as written, with the corrected Lemma S(A) constants".
- L334 [326] Step-4 review paragraph: "conditional on (F) and (O)".
- L370 [362] Related work: "conditional on Hypotheses (F) and (O)".
- L376-377 [368] Open items: new item for (O); (F) item re-labelled "The far-field gap" (its macro text unchanged).
- L379 [370] "the open step-3 items are (F) and (O)".
- L390 [381] AI-use statement: added sentence on the AI-assisted review (13 Nov) that found and repaired the Lemma S(A) omissions; remaining gaps now (F) and (O).

Not changed (deliberately): the (F) macro text, bibliography, section 4.3 caveats, kit/reproduction text, Step 3 "Documents" line (does not yet cite an erratum document), the stale "31 stages" sentence in section on computer assistance (the CLAUDE.md says 34; not part of this task, flag for PI).
- L351: manifest figures corrected to 34 stages / 97 upstream files (checked against FINAL_manifest.json; main.tex had 31 / 117, stale).

## (F1)+(O) update (13 Nov 2026)
Backup: main_v2_pre_F1O.tex. Line numbers are in the new main_v2.tex. Braces balanced (0), environments matched; not compiled. (F) was already split into (F1) windows / (F2) shells in the macro, so no split was needed.
- L12 Main Theorem title: "(F1) and (O)".
- L60 Abstract: "(Hypotheses (F1) and (O))".
- L64 Status box: "(F1) and (O)", plus clause that (F2) of the former (F) is now a consequence of repaired Lemma S(A) under (O).
- L90-91 New paragraph after the (F) environment: (F) has parts (F1), (F2); (F2) follows from (F1)+(O) via repaired Lemma S(A) (erratum in preparation), so remaining hypotheses are (F1) and (O).
- L93-99 (O) placeholder replaced by HYPOTHESIS_O.md LaTeX statement (O1)-(O3) verbatim, plus provisional-wording remark (may be revised (r2) with the uniformly local Lemma S(B) work; (O) not a corollary of (F1)).
- L103 Main Theorem: "Assume Hypotheses (F1) and (O)".
- L118 "closes under (F1) and (O)"; L121-122 remark title and item (a): (F1).
- L227 Figure node "Hypotheses (F1), (O)".
- L278 "under (F1) and (O)".
- L295 Subsection title "(F1) and (O)"; L296 "Hypothesis (F) is unchanged as a text, but (F2) is now a consequence of repaired S(A) under (O)".
- L303 "Hypothesis (F1) asks exactly for D3 and D4"; "under (F1) and (O) the bootstrap closes".
- L342 review paragraph "(F1) and (O)"; L378 related work; L385 open item relabelled "Hypothesis (F1)"; L387 "open step-3 items are (F1) and (O)"; L398 AI-use statement "(F1) and (O)".
Not changed: (F) macro text and its environment title, comments at L33, the stale "31 stages" flag from before.
