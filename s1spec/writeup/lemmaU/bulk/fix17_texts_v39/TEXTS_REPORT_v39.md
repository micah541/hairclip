# TEXTS_REPORT_v39: step-3 text set v39 (re-pinned step-4 records, step-4 wording, v4 kits)

Nov 20, 2026. In this report `D` = `/home/claude/s1spec/writeup/lemmaU/bulk/fix17_texts_v39/`.

The base is the v38 set in `fix16_texts_v38/`. It was read only; its `SHA256SUMS` still verifies.

**Step 3 is not closed.** Every status box in the set says so, and §4 lists what is still open.

- **What was written.** Only files under `D`. Nothing outside `D` changed (§3.4).
- **What was not done.** `step4_proof.tex` was not edited, and nothing in the proof tree was run.
- **Resources.** One process at a time, `nice -n 19`, `python3 -B`, `PYTHONHASHSEED=0`. The full rebuild with all checks, `tools/build_all_v39.sh`, takes about 4.6 s at about 61 MB peak RSS.
- **Not compiled.** LaTeX is not installed.

## 1. Files

| v39 file | from (unchanged) | change |
|---|---|---|
| `step4_handoff_v39.tex` | new | Step-4 hand-off as step 3 v39 reads it. §1: records of the re-pin. §2: hand-off list with (H1)–(H5) plus the v3.3 items, the (H4) range, **(H6) via C-1** (replaces the [A] sentences), and **(GS) added**. §3: the manifest_v38 stages that step 3 cites. §4: the kits, WSL test runs only |
| `STEP4_PATCH_NOTE_v39.md` | new | 10 proposed edits against the closed `step4_proof.tex` (sha256 6b4c183a…). Each old string is quoted verbatim, and the script asserts that it occurs exactly once. Not applied. |
| `step3_v39_main.tex` | `step3_v38_main.tex` | Status box, (H6) and (GS) items, the "Step-4 values" paragraph (now manifest_v38), and "What remains" (re-pin done; Windows runs open) |
| `main_v3_v39.tex` | `main_v3_v38.tex` | v39 status box, step-4 section: re-pin item, f12_sliver TODO resolved, Theorem 4(c) note, (H6)/(GS) rows in Table handoff, kits paragraph, reproduction section, open items, doc table rows R4v38, S4H |
| `lemmaU_v39.tex` | `lemmaU_v38.tex` | Status box; the S4 item now points to manifest_v38; the n1 stage-19 replacement; the M9 item; the date line; and two hash/path pairings (§2.3) |
| `continuation_v39.tex`, `lemmaS_A_v39.tex`, `flux2_lemma_v39.tex`, `GS_lemma_v39.tex`, `U1_regularity_v39.tex`, `C1_identification_v39.tex`, `G1_export_k2_v39.tex`, `bootstrap_v39.tex`, `chart_lemma_body_v39.tex`, `ei_lemma_v4c_v39.tex` | v38 namesakes | Record citations (manifest_v38 stage numbers) and the fixes of §2.3 |
| `E0_LEMMA_v4c_v39.tex`, `b12_preamble_v39.tex` | v38 | Names and header only |
| `numbers_A_v39.{tex,json}`, `ei_constants_v39.{tex,json}` | regenerated | Generators are copies of the v38 ones. All 122 macros and the EI table are identical to v38 (checked). |
| `erratum_v39.tex` | `erratum_v38.tex` | Parts 1–3 verbatim. New **Part 4, section O**: O1 records, O2 kits, O3 step-4 wording, O4 v38 retag damage, O5 display roundings, O6 the 0.1475 quote, O7 chart body, O8 status |
| `ERRATUM_v39.md`, `RECORDS_ERRATA_v39.md` | v38 | Verbatim, plus a pointer to Part 4 and a section "Status after the step-4 re-pin v38" |

The v38 file names of the set are retagged to the v39 names. The word "v38" itself is **not** retagged anywhere (see §2.3 for why).

### Tools (`D/tools/`)
- `gen_records_v39.py` writes `records_v39.json` from disk. It contains every hash, size, stage list, pin count and manifest result that the texts print, plus derived margins, and it asserts the two kit SUMMARY lines.
- `v39_values.py` holds the outward rounding: bounds up, margins down.
- `make_step4_v39.py` writes the handoff file and the patch note.
- `build_v39.py` does every other edit. Each edit is an asserted unique replacement, and the inserted text is recorded in `v39_inserts.json`.
- `check_balance_v39.py` and `check_numbers_v39.py` are the checks.
- `build_all_v39.sh` runs everything, and writes `SHA256SUMS` (24 files).

## 2. Changes
### 2.1 Records cited (task 1)

| record | value |
|---|---|
| manifest_v38 | `step4/p2/out_v38/FINAL_manifest_v38.json` (sha256 be909bbf42543ab017bb9c81ba58227c80fc4f010cc1b50a66a515b848f90418) |
| manifest_v38 contents | 29 stages; 629 paths and 55 scripts, i.e. 992 pin entries; 16 stages carried, 13 rerun or new |
| driver | `final_tail_v38.py` (b353337a…) |
| REPIN_REPORT | ba53dba3… |
| closed v3.2 manifest | e01cb4da…; the B8 file is accepted only through `EXCEPTION_v32_B8_check_outA5.json` (45091246…) |
| manifest_v33 | 8bc47e83…; one known F34 difference, `constants_needed_v34.json` (pinned 0f3c7455…, disk e70e70f0…) |
| `runs/P_A`, `runs/reg_fin13v4b` | accepted through `EXCEPTION_v4c_PA_reg_e2v4b.json` (ff0151c0…) |
| kit_step3_v4.zip | 0f899eb2de74…ce853 |
| kit_step4b_v4.zip | ab47b5416e8a…df7f0fa |

- **Kit runs (WSL test runs only):**
  - `kit3v4_run_wsl5`: 80 PASS, 0 FAIL, 1 WARN.
  - `kit4bv4_run_wsl1`: 159 PASS, 0 FAIL, 1 WARN.
- **Following your note**, every text says that these WSL runs are not independent, and that **no Windows run of either v4 kit exists**. The Windows runs are listed as open in each open list.
- **The texts now cite these results:**
  - G1_check rerun on chain v37 (stage 25; this was open in v38);
  - flux outputs pinned and reproduced (stage 21);
  - GS_tau1 / GS_boot / C1_ke / U1_sym pinned (stages 22–26);
  - f12_sliver on fin13v4c: 339/339 checks, M1′ ≤ 0.432858 ≤ 0.4329;
  - the n1 chart on fin13v4c (stage 19).
- Every "re-pin pending" or "TODO" sentence of v38 is replaced.

### 2.2 Step-4 wording (task 2)
- `step4_handoff_v39.tex` item (3): (H6) is Theorem C-1 (i), (iv), with C1_ke sup|κ̃e| ≤ 1.7429e-3 (×286 below ½).
- It replaces these step-4 [A] sentences:
  - "smooth for t>0 by parabolic regularity" (l.672), and "(H6) follows from (b)" (l.131);
  - "continuation uses the literature fact (P1)" (l.130, 336), and "continuation fact (P1) … from the literature" (l.370).
- Item (4) adds **(GS)**:
  - sup|ku| ≤ 6.5170e-3 (×76.7) and sup|u′| ≤ 0.054182 (×18.4), certified on |x| ≤ 60.25 (stage 23);
  - the step-3 strict improvement 5.3903e-4 / 0.040936 (stage 24) is given for reference.
- Not removed: l.789/791 (T maximal by (P1)). That is step 3's own criterion (P1′) and keeps its own open items; this is stated in the file and in the patch note.

### 2.3 Fixes of v38 text slips (found while doing this; none changes a number or an argument; erratum O4/O5)
1. **v38 builder retag damage.** Its rule `'text v37'→'text v38'` rewrote historical sentences. "text v38 printed/had/sketch/Audited since" is restored to v37 in these places:
   - continuation ×3;
   - U1;
   - ei_lemma ×2;
   - main_v3;
   - step3 main.
2. The same rule produced three more slips, all restored:
   - `runs/reg_fin13v4c` (that directory does not exist; the correct name is `reg_fin13v4b`) in step3 main;
   - the GS_lemma review heading pointed to its own successor file instead of `fix12_texts/continuation_v37.tex` l.41;
   - the lemmaU date line referred to itself.
3. `lemmaU_v38` printed the sha256 of `fix10_proofs/G1_export_k2.tex` (4004a774…) and of `fix10_proofs/U1_regularity.tex` (827a5f51…) next to the paths of the v38 copies. The pairings are now explicit.
4. `ei_lemma_v4c.tex` (v38) contained a NUL-delimited builder placeholder in a header comment, so `grep` treats the file as binary. The sentence is restored, and the balance check now rejects control characters.
5. C1_identification printed `sup_ke = 1.7428e-3` (a truncation). It now prints ≤ 1.7429e-3.
6. The GS strict-improvement bounds are printed outward as 5.3903e-4 / 0.040936. The 5.3902e-4 / 0.040935 quoted in the task note are truncations, not upper bounds.
7. `chart_lemma_body` named fin12f4 and the v33 n1 file. It now names fin13v4c and the stage-19 output. The printed bounds are unchanged and still hold, as the number check confirms.
8. `lemmaS_A`: the `base_ibp` outputs are now pinned by sha256 (stage-20 inputs) but still carry no args block.

### 2.4 Stale 0.1475 (task 3)
- No v38 or v39 text contains 0.1475; the checker asserts this.
- `lemmaU_v39` prints √2·0.10434 ≤ 0.1476.
- The stale quote survives only in the superseded `fix1/lemmaU_rev/lemmaU_v34.tex`, which is unchanged. Stage 20 of v38 checks 0.1476 instead.

## 3. Checks (all in `tools/`, results reproduced by `build_all_v39.sh`)
### 3.1 Balance (`check_balance_v39.py` → `balance_v39.json`): **23 PASS, 0 FAIL**
- **Per file:** braces, `\begin`/`\end` and `$` parity, plus duplicate labels. All 19 files pass, and so does the `\input`-expanded `step3_v39_main.tex`, where every `\ref` resolves.
- **Cross-document:** every `\ref` of a standalone file resolves to a `\label` somewhere in the set (222 labels), and every `\input` target exists.
- **Control characters:** none.
- **Paths:** all 92 file paths named in the new v39 text exist. The one exception is allowed on purpose: the wrong name `runs/reg_fin13v4c`, which erratum O4 quotes.

### 3.2 Numbers (`check_numbers_v39.py` → `numbers_check_v39.json`): **PASS**
- **Inherited from the v38 checker:** KEYED 386, FILE 128, VALUE 1113, WHITELIST 137, RECORD 15, HISTORICAL 671, STALE 0, SUSPECT 0, UNMATCHED 0.
  - 122 macros pass.
  - The v38 ledger now also includes every JSON under `step4/p2/out_v38` and `records_v39.json`.
  - The checked texts now also include `step4_handoff_v39`, `chart_lemma_body_v39` and the patch note. This added two whitelist entries for the chart body: the formula constant 6.25 and the step-1 constant X0.
- **Pinned rows: 56/56 pass.** These are direction-aware. 26 of them are new and keyed directly to `FINAL_manifest_v38.json` results:
  - H7 ratio, t*, margin;
  - far H4;
  - flux;
  - C1_ke, GS_tau1, GS_boot;
  - f12 sI2′, M1′;
  - (H2) for f−U₁ and (H2f): bounds and margins.
- **(8) Numbers in new v39 text:** 66 match manifest_v38 data, 3 match other current files, 4 are whitelisted with a reason, and **0 miss**.
  - The 3 current-file matches are erratum quotes of old values.
  - The 4 whitelisted numbers are: the kit WARN 5.7e-16, the horizon 27.15, and the erratum quotes 0.28578 and 0.1476.
- **(8b) Direction:** a number after `\times` must be ≤ the value, and a number after `\le` must be ≥ the value. **0 fail.**
- **(9) Hashes:**
  - 89 known: in manifest v38/v33/v32, the disk sha of a record, a SHA256SUMS file, or fix10/fix14/fix15 files.
  - **0 unknown in new text.**
  - 12 inherited ones are unverified: v34/v35-era hashes in lemmaU, and the Lean source hashes in U1.
- **(10)** 0.1475 absent.
- **Negative test**, on a scratch copy that was then deleted. One changed hash digit gave NEW-UNKNOWN; a margin raised to 15.59 gave NEW-DIR-FAIL; a bound lowered gave PINNED-FAIL. All three were caught.
- **Limits:**
  - Value matching proves only that a value of that size exists in the ledger. The direction check can pass when some other ledger value matches. Only the 56 keymap rows are tied to a specific key.
  - Numbers with fewer than 3 significant digits are not checked by the v38 part. The new-text check also covers 2-digit numbers.

### 3.3 Reproducibility
- Two consecutive `build_all_v39.sh` runs gave byte-identical outputs: `SHA256SUMS` matches for all 24 files.
- Numbers macros and the EI table are identical to v38.

### 3.4 Nothing outside `D` changed
- `find s1spec s1_* -newer <start marker>`, excluding `D`, finds 0 files. Only the parent directory's mtime changed, because `D` was created in it.
- `fix16_texts_v38/SHA256SUMS` verifies.
- `step4_proof.tex` is 6b4c183a…, unchanged.
- manifest_v38 is be909bbf…, unchanged.

## 4. Open (step 3 is not closed)
1. **Lean:** `lake build` of CWFlow at commit d72b744, and an n = 1 instantiation file.
2. **DKS re-read:** pp. 1234–35, (D1), and the uniform bounds on ∂ₛᵐκ used in C-1.
3. **LaTeX compile** of every v39 file.
4. **Author confirmations:** (P3′), Lemma O arc by arc, the single-point statement.
5. **Final independent review** of the v39 text set.
6. **Windows runs** of `kit_step3_v4.zip` and `kit_step4b_v4.zip` by the authors. Only the WSL test runs exist.
7. **Step-4 document:** a new step-4 version that applies `STEP4_PATCH_NOTE_v39.md`. This needs the authors' go; `step4_proof.tex` itself stays closed.
8. **Records** (`RECORDS_ERRATA_v39.md`, status section), all minor and not affecting a bound:
   - R2: the third cert pass has no driver record;
   - R3 / B-v37-1: the V37_LOG pin, and `SHA256SUMS_v37` not regenerated;
   - B1: `base_ibp` outputs lack an args block;
   - B2–B4 and C-7;
   - ramp_U's internal `u1_ghat_sha`;
   - GS_boot items GS-2 and GS-6 not applied;
   - optional reruns of `runs/P_A` and `runs/reg_fin13v4b`.
9. **[read]/[hand] items** of the verifier lemmas and C-1, as in v38:
   - DS0–DS3, K_ρ, K_W;
   - W monotone on [31.25, 34.25];
   - the tail monotonicity statement;
   - the D-1a pins;
   - the regularity of the base in τ at slab ends;
   - 98 inherited-unverified step-1/2/4 quotes.
