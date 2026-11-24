# STEP 4 RESULT, v3.3 addendum (new hand-off items for step 3 v3.3)

Addendum to `/home/claude/s1spec/step4/STEP4_RESULT.md` (base v2). That file is **not edited**; `STEP4_RESULT.md` stays authoritative for everything not listed here. Nov 13, 2026. Source: `step4/p2/v33/T6_REPORT.md`; hashes in `V33_INDEX.md` section 1.7.

Status. Step 3 v3.3 (Lemma U, repaired Lemma S(A); `writeup/lemmaU/bulk/`) needs hand-off data that step 4 v2 did not certify: a datum on the far support [60, 700], a tightened shell threshold, and the fixed-L cone cylinders. All of it is certified below. **Closed step 4 is untouched:** `p2/out/FINAL_manifest.json` (sha256 e01cb4da71ba825318b5395f053595d6e980e84c1ac286b9d72ff1e4376a8c94; 34 stages, 523 paths, 97 upstream, 26 scripts) verified before and after the v33 run, and a size/mtime snapshot of all 10,811 files outside `p2/v33/`, `p2/out_v33/`, `p2/out_v33_pass1/` is identical. Everything is arb with outward rounding, no sampling. Labels: **[cert]**; margin = threshold / bound, rounded down. Not yet independently reviewed (task 8 of the bulk plan).

## 1. New hand-off items (uniform for tau* in slab 180; `p2/v33/`, `p2/out_v33/`)
At Ybar = 1e-12 (certified for Ybar <= 1.7e-12, the level of the n3 boxes).

| item | required | base (arb, slab 180, max over 32 t_c balls) | with error (Ybar = 1e-12) | margin | Ybar limit |
|---|---|---|---|---|---|
| **(H2f)** j = 0, G_0 = int_{ln 60}^{ln 700} (28/y)(f - U_1)^2 ds | 1e-6 | 3.1263e-7 | 3.1363e-7 | **x3.19** | 5.7e-8 |
| (H2f) j = 1 | 5e-6 | 1.0031e-6 | 1.0034e-6 | **x4.98** | 1.7e-8 |
| (H2f) j = 2 | 1e-5 | 2.4302e-6 | 2.4351e-6 | **x4.11** | 1.17e-9 |
| **tightened (H2)** shells [14, 60], j = 0 | 1e-8 | 8.010e-11 | 8.013e-11 | x124.8 | |
| tightened (H2) j = 1 | 1e-6 | 1.335e-10 | 4.71e-10 | x2125 | |
| tightened (H2) j = 2 | 1e-3 | 2.936e-10 | 1.559e-6 | x641 | |
| **fixed-L cone** (H5) slope datum | 0.08 | | 0.03029 (worst row) | **x2.64** | |
| fixed-L cone (H5) sup \|h_0''\| = K2 | 0.03 | | 9.876e-3 | x3.04 | |
| fixed-L cone (H5) \|\|h_0'''\|\|_L2 = N3 | 0.004 | | 1.461e-4 | x27.4 | |
| fixed-L cone H3 jets K2 / K3 / K4 | 0.12 | | 9.88e-3 / 5.96e-5 / 8.36e-4 | x12.2 / x2012 / x143 | |
| fixed-L cone Lemma O H4, int (h_0^(5))^2 | 0.0144 | | 1.286e-5 (cone r_c 300-310; chart factor 1.0031) | x1119 | |
| fixed-L cone Lemma O H1 (graph) | sup\|h_0\| <= H - C_0 sqrt(t_1) | | graph on all 8 rows: xi <= -179.1 at arc start, >= +170.0 at arc end, XL = 150; sup\|h_0\| <= 0.2162, H <= 0.2185 | holds | |
| fixed-L cone H7 time ratio | t* >= ratio | | ratio <= 1.6000014e-7 (a-bar = 2.19103e-7); t* >= 3.68334e-7 | **x2.30** | |

Notes.
- (H2f) is the new threshold (5e-7, 5e-6, 1e-5) of the r2 statements, with **(H2f)_0 raised to 1e-6** and beta_A = 2 in v3.3 (so (1e-6, 5e-6, 1e-5)). It is the datum of the repaired Lemma S(A) on [60, 700], together with the tightened (H2) = (1e-8, 1e-6, 1e-3). Parts of the level-0 base value: [60, 379.5] 4.84e-8, [379.5, 700] 2.65e-7. The error term changes G by at most 0.3%.
- Method: win.py's cell machinery on x in [45, 725] (fidelity asserted: identical to `win.energies` on [45, 390]; an unwidened ball reproduces the pinned `win_0.json` bit for bit); 32 IFT balls widened by (x_c, lambda) +- 3.28e-11 and alpha +- 1.23e-12; the e-term is summed in L2 over 51 bins of ratio 1.05 (Minkowski), not maximised over sliding windows (T6 section 1).
- Fixed-L cone family: L = 50, X = 3, centres r_c in [300, 380] in 8 rows of width 10, cylinders [r_c - 150, r_c + 150], in L = 50 units. This replaces the unit-L cone family of v3.2 in step 3 v3.3 (H5 of v3.2 is unchanged and still holds for v3.2).
- The thinnest new margins: fixed-L cone H7 t* **x2.30**, slope **x2.64**, (H2f)_0 **x3.19**. The overall thinnest margins of step 4 are unchanged (|alpha| x1.27, ...; STEP4_RESULT section 1).

## 2. The new manifest
- Driver `p2/v33/final_tail_v33.py` -> **`p2/out_v33/FINAL_manifest_v33.json`** (sha256 52af64e6c5fb923a51429ce7b98bf536502df6338868897cde8f9c1a978e3ca6): 10 stages, 152 paths, 23 scripts. Written atomically; resumable with `--from`; progress in `FINAL_manifest_v33.partial.json`; children run `python3 -I -B` under `nice -n 19`.
- Stages: 0 h2f_ift, 1 h2f_balls -> h2f_win, 2 t9_h2f, 3 h2f_camp, 4 h2t_v33, 5 t9_cone_v33, 6 h4loc_v33, 7 t9_lemmaO_v33 (LO_NS = 48), 8 ratio_v33, 9 notes_v33. Total about 280 s.
- **Relation to FINAL_manifest.json.** The new manifest is additional and does not replace it. Its first action verifies FINAL_manifest.json (34 stages, 523 paths, 97 upstream, 26 scripts); its last action re-verifies it. 52 of its inputs are pinned v3.2 files, each checked against the hash that FINAL_manifest.json records for it. No stage reads a v33 output back into a v3.2 stage. A complete v3.3 step-4 record is therefore the pair (FINAL_manifest.json, FINAL_manifest_v33.json).
- Reproduction (T6 section 6): pass 1 kept in `p2/out_v33_pass1/`; pass 2 is the clean run the manifest records. 56 of 74 JSON outputs are byte-identical between passes; the rest differ only in `cpu`/`seconds` or in recorded hashes of per-ball files that differ only in `cpu`. Pass 1's only failure was the H1 end test, see section 4.
- Fidelity of the copies: with v3.2 geometry, `t9_cone_v33.py` and `h4loc_v33.py` reproduce the pinned v3.2 cone rows bit for bit; `t9_lemmaO_v33.py` reproduces the v3.2 sup|h_0|, H and H8 to about 1e-9 relative, not bitwise (likely python-flint 0.9.0 transcendental enclosures).

## 3. Reused v3.2 outputs (not rerun)
| item | pinned output | status |
|---|---|---|
| far family (centres >= 490 in v3.3; a subset of the certified v3.2 family, centres >= 376) | `t9_far_v32.json`, `h4loc_v32.json`, `t9_lemmaO_far.json` (ok) | reused; bases start at r >= 264 >= 150, inside t9_far's region |
| H6 (v3.2) | `t9_h6_v32` | **superseded** by step 3 v3.3's H6 (`zone3_interface.json`, cone_family.H6, minimum x1.86 at the reference Bhat); its pinned output stays in FINAL_manifest.json and is not a v3.3 input |
| (H1) alpha, sqrt(Lambda), (H3) windows and K4, (H4) graph, (H5) far data, t1 interval, all of STEP4_RESULT sections 1-2 | as in STEP4_RESULT | unchanged |
All of these are recorded in `out_v33/notes_v33.json`.

## 4. Arc-end note (H1 of Lemma O; for review)
The v3.2 script `t9_lemmaO.py` computed the arc-end cells but **never tested |xi| > XL there**. The v33 copy asserts it. With v3.2's arc of 1.1 XL the end test fails by a few units (crude enclosures: the width-10 centre row and slab-time motion), so the arc was lengthened to 1.3 XL; cos(theta - theta_c) > 0 is checked on the whole longer arc, so the test is stronger, not weaker. sup|h_0| of about 0.2 reflects the same enclosure widths; Lemma O puts no other condition on H. H8 is informational only (distance >= 0.72 against the old 0.1 rule). **Open question for the review (not decided here): does closed v3.2's cone H1 need an erratum**, because its script did not check the arc ends? The far family's H1 (`t9_lemmaO_far.py`, STEP4_RESULT section 3) is a separate script and does check |xi(+-G)| >= 226.18 > XL.
