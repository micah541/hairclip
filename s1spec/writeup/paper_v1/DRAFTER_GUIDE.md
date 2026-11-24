# Drafter guide (paper v1): one page

Read first: `OUTLINE.md` (your unit's row in §7.2, its section in §2/§3), then this page, `notation.md`, `namemap.md`. Nothing existing is edited; all output goes to `/home/claude/s1spec/writeup/paper_v1/`.

## Files
| file | owner | use |
|---|---|---|
| `preamble.tex` | coordinator | `\input` by `main.tex` and `supplement.tex` (no `\documentclass` inside); add macros **by request only** |
| `main.tex`, `supplement.tex` | coordinator | masters. The supplement master must say `\newif\ifsupplement\supplementtrue` before `\input{preamble}`. Compile from `paper_v1/` |
| `numbers/numbers_paper.tex`, `numbers/ledger.json` | U0 | number macros; read-only for you |
| `notation.md`, `namemap.md` | U1 | symbols and archive->paper names. A symbol or name you need that is missing: ask, do not invent |
| `partA_conventions.tex` | U1 | supplement Part A |
| `ISSUES.md` | everyone | append-only gap log (below) |
| your unit's files | you | e.g. `sec2_objects.tex`, `partE_shells.tex`; one `.tex` per section or Part, no `\documentclass`, no `\begin{document}` |

Balance check (no LaTeX on this machine): `python3 -I /home/claude/s1spec/writeup/lemmaU/bulk/t3_proof/check_balance.py yourfile.tex`. Expect "undefined refs" only for labels owned by other units and for `#1` in macro definitions.

## Rules (OUTLINE §5, condensed)
1. **Numbers.** Only macros from `numbers_paper.tex`. A literal certified number is an error. Same number in two places: same macro. Rounding: upper bounds up, lower bounds and margins down; margin = threshold/bound (printed through a number macro; `\mgn{<macro>}` typesets the factor sign). If two archive documents disagree, the pinned JSON wins; report the disagreement in `ISSUES.md`.
2. **Trust classes.** Supplement: tag `\tcert`, `\test`, `\thand`, `\tlit` ([cert] [est] [hand] [lit]). Main paper: certified is unmarked; non-certified numbers carry `\test` ("(floating point)") and never appear in a theorem. Every theorem ends with its class: "The proof is computer-assisted; see Section 8." (for hand proofs: say so).
3. **Names.** Use the paper names of `namemap.md`: (H1)-(H7) via `\Hh{n}`, (O1)-(O9) `\Oh{n}`, (L1)-(L3) `\Lh{n}`, (B1)-(B5) `\Bh{n}`, (K1)-(K6) `\Kh{n}`, (W1)-(W5) `\Wh{n}`, Phase I/II `\PhaseI`, `\PhaseII`; Theorems A-D by `\begin{thmlet}`; Main Theorem `\begin{maintheorem}`. **Forbidden outside Parts J and K:** "Hypothesis (F)", "(O)" as a hypothesis, B0, M9, N1, G-A1, fin12f4, v3.1/v3.2/v3.3/v35, review codes, stage and file names, "H2f", "P1'".
4. **Notation.** One meaning per symbol (`notation.md`). Renamed archive symbols are never used in their old meaning: log-radius is `s`, arclength of the flow curve `\arclen` (`\kss`), rate `\rate` ($\varsigma$), scale $\lambda$, spectral parameter $\mu$, slow eigenvalue `\muslow`, energy weight `\muI`, centre `\rc`, Lemma O constants `\MO`, `\KO`, cylinder unit `\cyl`, singular time `\Tsing`, profile curve `\Gsh`, parameter of $\Gamma_0$ `\vth`. Macro names are in `preamble.tex` ("NOTATION MACROS").
5. **Paths.** No file path or JSON key in the main paper (`\fn` warns). In the supplement, only in Parts J/K or in the "Certified values" lines of a lemma: `\fn{path}`, `\key{dotted.key}`.
6. **Consolidation (D5, D7, D15).** No erratum, no version history, no "as in v3.2 §x", no intermediate thresholds in the text of Parts D-I; one proof, final hand-off list only. Rewriting for exposition only: no statement, constant or script changes.
7. **Style.** Plain declarative sentences; no status boxes, dates or "pending"; do not paper over a gap.
8. **Environments.** `theorem proposition lemma corollary definition convention remark` (numbered by section; in the supplement the Part letter replaces the section number); `hyplist` for lettered condition lists: `\begin{hyplist}{H} \item ... \end{hyplist}` gives (H1), (H2), ...

## Labels and cross references
Labels: `<unit>:<kind>-<name>`, unit = `intro obj s3 ... s9` (main paper) or `pA ... pK` (supplement), kind = `thm lem prop cor def eq tab fig sec`; e.g. `pE:lem-chart`, `s6:prop-closure`. Never reuse another unit's label. Same document: `\Lref{label}` (plain `\ref`). Other document (xr-hyper, both directions): `\Xref{pE:lem-chart}`, giving the number, so write "Lemma~\Xref{pE:lem-chart} of the supplement" or `\SuppPart{E}`. Citation keys: only those of `/home/claude/s1spec/writeup/main/main.bib` (the DKS paper is `DziukKuwertSchatzle2002`); a new reference needs a verified entry first.

## Reporting a gap
If a statement, constant, source or name cannot be settled from the sources: (1) stop that paragraph with `\GAP{one sentence: what is missing}` (it prints in red and counts; `\finalpapertrue` turns every `\GAP` into an error); (2) append to `ISSUES.md`: `- [unit, file:label] what; what you checked; what is needed (who can decide)`. Never fill a gap with a guess, a rounded number or a changed constant. Report at the end of your run: the file list, ledger macros used, and the `ISSUES.md` entries added.
