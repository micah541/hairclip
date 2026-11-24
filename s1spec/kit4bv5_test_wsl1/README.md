# Reproduction kit 4B v5: step 4 (v3.2 chain), the 21 v3.3 stages and the 29-stage v38 re-pin chain

**New in v5 (against v4): the Windows fixes** of the first Windows run of v4 (Nov 24 2026, 159 PASS / 0 FAIL / 1 WARN after run-copy patches; `CHANGES_v5.md`): (1) the write guard splits the joined Windows command line (CommandLineToArgvW); (2) `chart4_shell.py` runs without the Unix-only `resource` module (its pin in `inputs/MANIFEST.sha256` is regenerated); (3) the flux texts and their recorded hashes are compared modulo CRLF line endings (`flux_v38.py`, `verify_v33.py`); (4) `verify_v4c_runs_with_exceptions.py` joins relative paths with `/`. No certified computation changed.

**Windows 11 (PowerShell).** Tested: CPython 3.13.7 from python.org, a venv, python-flint 0.9.0, numpy 2.5.4, sympy 1.14.0, mpmath 1.3.0, scipy 1.18.1 (`requirements.txt`). Unzip to a short path outside OneDrive, e.g. `C:\kitruns\kit4bv5`:
```powershell
Expand-Archive -Path .\kit_step4b_v5.zip -DestinationPath C:\kitruns\kit4bv5
cd C:\kitruns\kit4bv5
py -3.13 -m venv kitenv
.\kitenv\Scripts\python.exe -m pip install --upgrade pip
.\kitenv\Scripts\python.exe -m pip install -r requirements.txt
.\kitenv\Scripts\python.exe run_all.py --jobs 4          # 5-7 h; resumable: run the same command again
.\kitenv\Scripts\python.exe check.py
```
Expected: `TOTAL PASS`; 1 WARN (`iface_compare_v38.json` has 2 fewer checks than the proof tree's, by design: C1_ke / GS_tau1 read this run's own step-4 data). **The INFO line "outputs vs proof tree: byte-identical / timing-fields-only / provenance-only / different / missing" is INFO, not a failure**: with other floats (another platform / libm) almost every regenerated v3.2 file differs from the proof tree's in the last digits (Windows v4 run: 0/0/0/498 different; WSL 30/0/0/468); every claim is recomputed on this run's own outputs and checked as an inequality.

**New in v4 (against v3, below):** (1) run_all stage `v38` runs `step4/p2/v38/final_tail_v38.py`: the step-4 hand-off re-pinned on verifier fin13v4c, (A) chain v37,
flux lemma v2 and the interface / proof checks (29 stages; proof-tree manifest `FINAL_manifest_v38.json`, compared in `check.py` section v38). (2) **Write guard**:
`run_all.py` loads `kit_guard.py`, so every process of the kit refuses to write, rename, delete, change directory or start a child outside this folder (the
Nov 15 2026 incident); files read outside the kit are listed by check.py. (3) Checker fix: path-valued keys are normalised (the v3 WARN on f12_sliver.json).
Run it exactly as v3: `run_all.py --jobs 6`, then paste the SUMMARY of `check.py`. Expected time: v3 + about 15 min.

---

# Reproduction kit 4B v3: the rest of step 4 (v3.2 chain) and the 21 stages of the v3.3 hand-off chain

**What is new compared with kit B (kit_step4b.zip, 18a41a7c...).** (1) The six kit bugs of the first end-to-end run (`kit4b_run_wsl1/RUN_NOTES.md`, 76 PASS, 0 FAIL) are fixed in the driver
and the checker, and one more bug found while porting (7, below). (2) The 21 stages of the v3.3 step-4 chain (`step4/p2/v33/`, proof-tree manifest `FINAL_manifest_v33.json`, sha256
8bc47e83...: (H2f) on [60, 700], tightened (H2), the fixed-L cone family, B0, the M9 arb stages, N1 chart terms, f12_sliver) run after the v3.2 chain, and `check.py` verifies their margins
(STEP4_RESULT_v33_addendum, sections 1, 5, 7) and compares every output with the proof tree's (`reference/v33_out`). The kit reuses nothing from earlier runs.

### The kit fixes (all marked `[kit fix N]` / `[kit patch]` in the sources)
1. `t7_b`: `step4/t7/se_widths_secorr.json` is regenerated (`se7.py widths --source se_corr`, renamed) before the t6c-source widths step; it is not shipped. (The proof tree's file differs from the regenerated one at relative ~1e-13.)
2. `h1pts` (the 16 exact H1 points, `h1.py`) is a stage before `h1dt`.
3. `tr()` rewrites the proof tree's path prefixes in a single pass (a kit unpacked below `/home/claude/s1spec` was doubled).
4+5. `tail1` = the pre stages + the `T8 ...` stages of pass A only, then it returns (`n3boxes` and `t9` need the nonlinear closure that `t7_nl` publishes); `tail2` runs the whole tail.
6. The B8 stage `b8` (`check_outA5.py`) is in `STAGES`; `verify_outputs.py` resolves the reference path `../../s1_tight/consistency/check_outA5.json` to the kit's `s1_tight`.
7. (new) The tail's own last stage (B8) has the working directory `../../s1_tight/consistency` in the manifest; the old lookup resolved it to `/home/claude/s1_tight/consistency`, OUTSIDE the kit, so in the WSL test run it ran in (and overwrote the output of) the proof tree. Now it resolves to the kit's `s1_tight`, and `tail_run` asserts that every stage's cwd lies inside the kit.
- `check.py`: the round-1 vs round-2 `n3_b1` comparison drops the `args` block (the numbers were identical; only argv/paths/hashes differ).

Independent of kit A: it starts from our reference P1 checkpoints (slab 95 for the base seed, slab 96 for the switch), shipped in `step4\t3\p1\ckpt\`.
Everything is plain Python 3.13 + python-flint (arb), one single-threaded BLAS per process, no compiler, no network. Kit A and kit B may run at the same time
(kit B uses `--jobs 6` worker processes by default, kit A one core).

## Run it (Windows 11, PowerShell, unzipping to `C:\kit4b`)

```powershell
Expand-Archive -Path .\kit_step4b_v3.zip -DestinationPath C:\kit4b
cd C:\kit4b
dir run_all.py                      # must exist here
py -3.13 -m venv kitenv
.\kitenv\Scripts\python.exe -m pip install --upgrade pip
.\kitenv\Scripts\python.exe -m pip install -r requirements.txt
.\kitenv\Scripts\python.exe run_all.py --jobs 6
.\kitenv\Scripts\python.exe check.py
```
`run_all.py` runs all 39 stages in order (up to `--jobs` worker processes at a time, default 6), then runs `check.py` itself. **Paste the check output back** (the SUMMARY block at least); no zip is needed.
Options: `--jobs N`, `--list` (stages and state), `--only a,b`, `--from NAME`, `--force a,b`, `--status`. There is no quick mode: the stages depend on each other and on the full base.
Unzip to a short path outside OneDrive; plug the PC in and stop it sleeping (`powercfg /change standby-timeout-ac 0`).

## Stages

| stage | what |
|---|---|
| inputs | SHA-256 check of the shipped data inputs |
| t1_seed, t1_float, t1_polish, t1_se | base v2: seed from P1 slab 95 (arb 192), float base (181 slabs), arb-128 polish (**serial**, resumable at slab granularity), (SE) offsets |
| t6, rhoend, se_corr | T6 residual enclosures of the polished base; far-end residual sup; (SE) correction (its widths are used for slabs 0-1 only) |
| t6c | B1: T6's engine run on the (SE)-corrected base itself, slabs 2-180 (residuals.json and jets): the one T6 data set of the final tail |
| t7_a, t7_b | T7 certificates on the polished base (+T8), then on the corrected base (se7 on t6c), published to `step4\t7\constants.json` |
| h1pts | the 16 exact H1 points (`h1.py`), read by the H1 summary at the end of `h1dt` |
| h1dt, winslab | slab-180 hand-off: 32 derivative balls, 32 window/shell balls, with their summaries |
| tail1 | round 1 of the final tail (`reference\final_manifest.json`): the pre stages (switch, (Y), T8 cconst/merge, clam, the 16 H1 points and H1 summary) and the `T8 ...` stages of pass A, which make `n3_b1`; nothing else (the rest needs T7's nonlinear closure) |
| t7_nl | T7 republished with the nonlinear closure (Lemma 4c) on tail1's `n3_b1` (the proof tree's T7 is circular with T8 in this way) |
| tail2 | the final tail again on the republished T7: this is the reference chain (34 stages: t9 bootstrap, passes A and B (camp / n3far iteration, n3boxes, t9), t9_far, h4loc, Lemma O, v3.2 scripts, t1clock, axis crossing, and last the B8 check). `check.py` verifies that round 1 and round 2 give the same `n3_b1` |
| b8 | the B8 hull check `check_outA5.py` once more, standalone (`s1_tight\consistency\check_outA5.json`) |
| v33_00 ... v33_20 | the 21 stages of `step4\p2\v33\final_tail_v33.py` (one run_all stage each, `--from i --until i`), in this order: h2f_ift, h2f_balls (32 balls), t9_h2f, h2f_camp, h2t (tightened (H2)), cone (fixed-L), h4loc, lemmaO_h1, ratio (H7), notes, b0_dprof, h2u, m9_balls (arb, 64 sub-runs), h3_arb, h2t_arb, h2u_arb, h2f_camp_arb, lemmaO_arb, h4far, n1_chart, f12_sliver. They read this kit's own v3.2 outputs (`p2\out\h1slab\ift_*, win_*`, `t9_cone_kch.json`, `n3_b1.json`, `amp_A1_w180.json`, the t6c jets, `t7\constants.json`) and write `step4\p2\out_v33\` incl. `FINAL_manifest_v33.json` |

The console shows `[stage i/39]`, a heartbeat every 3 minutes and an ETA. Expected time: about 3 to 4.5 hours on the i7-10700 (the v3.3 stages add about 15 to 30 minutes; the proof tree's serial v3.3 chain took 13 minutes), about 5 hours at `--jobs 2` on a 2-core WSL sandbox (see the end of this file).

## Interrupted? (Ctrl-C, reboot, crash)

Run the same command again. Finished stages (`kit_out\state\*.done`) are skipped; inside a stage every tool is restartable (the polish continues at the first unstored slab with `--resume`; T6, t6c, T7, rho_end, se_v2, h1, h1dt, win skip finished units). A second `run_all.py` in the same folder refuses to start.
The polish after a restart builds a fresh float Jacobian, so Newton paths after the restart slab differ from an uninterrupted run, within the tolerance |E| < 1e-26.

## The v3.3 stages: what is pinned, what is regenerated

* Regenerated by this run: every output of the 21 stages, from this run's own v3.2 chain (so the v3.3 numbers can differ from the proof tree's in the last digits, exactly as the v3.2 outputs do).
* **f12_sliver is the exception: it recomputes the (A) / U / zone III chain from PINNED inputs only** (the proof tree's `b0_dprof.json`, `n1_chart.json`, `h2u_v33_arb.json` and the closed
  step-3 v3.3 files under `writeup/lemmaU/bulk/`), and asserts that each harness reproduces the pinned file bit for bit (modulo the installation path). It therefore reads the proof tree's three step-4
  files from `reference/v33_out/p2/out_v33/` (SHA-256 in `inputs/MANIFEST.sha256`) and not this run's stages 10, 15 and 19; this run's own `b0_dprof`, `n1_chart` and `h2u_v33_arb` are still compared
  with them by `check.py`. In short: f12_sliver checks the step-3 side of the chain from fixed data; the other 20 stages check step 4's hand-off data on this run's own chain.
* The driver `final_tail_v33.py` of the proof tree verified the closed v3.2 manifest byte for byte before and after. The kit regenerates the v3.2 chain, so it records instead which of the v3.2 files it
  reads differ from the proof tree's (`checks.v32_pins_differing_from_proof_tree` in `FINAL_manifest_v33.json`); a missing input is a hard error.
* Children run `python -I -B -X utf8` with an explicit environment (as in the proof tree; `-X utf8` added for Windows).
* `check.py` section "v3.3": the thresholds of the addendum tables with the addendum margins as WARN-level floors, then the output comparison (EXACT / ARGS / NUMERIC / DIFFERENT / MISSING, `verify_v33.py`).

## What `check.py` compares

1. **Every certified inequality** of the final tables of `STEP4_RESULT.md` (sections 1, 2 and 6: switch and (Y), T7, the T8 limits and N3' boxes, H1-H5, Lemma O H4 and H1, the v3.2 data and closure margins), recomputed from this run's outputs and evaluated strictly (no tolerance), with the proof tree's value printed next to yours. The list is `step4_claims.py`; the proof tree's values are in `reference\ours.json`.
2. **Every output against the proof tree** (`reference\outputs.json`, SHA-256 of about 500 files): byte-identical or different. Files can legitimately differ in bookkeeping fields (cpu time, paths); the list of differing files is printed, and every claim above is recomputed on your own outputs anyway.
3. T1 / T6 norms against the T5 targets.
`kit_out\check_report.txt` is the same text.

## Inputs shipped as data, with SHA-256 (`inputs\MANIFEST.sha256`, verified by the first stage)

These cannot be regenerated on Windows or belong to earlier steps:
* `step4\t3\p1\ckpt\slab_00095.*`, `slab_00096.*`: our P1 reference checkpoints (kit A recomputes its own).
* `s1_tight\prof\prof_tight.txt`, `s1_tight\ext813\ext813_prof.txt` (and `ext813_points.txt`), `s1_tight\ext_tip\prof_tip.txt` (and `ext_tip_check.json`): the CAPD-based step-1 profiles (C++/CAPD, not rerun; used by the S1 chart of H1 and of the windows).
* `step3\outA5.txt` (the Lyapunov CAPD hull boxes behind the 3425 step-3 IBP boxes, regenerated at run time by `ibp_boxes.py` from this file), `step3\hs_uforms.pkl`, `step3\lyap_consts2.json`, `step3nl\agent_mod\lemma\lyap_extra.json`.
* `writeup\lemmaU\bulk\...` (the closed step-3 v3.3 inputs of f12_sliver: `t2_A_v35`, `t4_cert`, `t5_zone3`, the fin12f4 verifier files in `step3nl_v33\`) and `reference\v33_out\` (the proof tree's v3.3 outputs, also the three pinned inputs of f12_sliver).
* `step4\t7\phi_store.json` (the fixed functions phi_1, phi_2: DOFs of the exact midpoints of the v1 slab 187, an earlier base) and `step4\t7\fixed_coarse.json` (lossless cache of the coarse tip mesh, keyed by the tip knots).

## Differences from the reference proof tree (read this)

* **Earlier point enclosures for the H1 derivative balls.** In the proof tree the `h1dt` balls and the window/shell inputs (`win_slab`) were computed from the point enclosures that `h1.py` gave *before* its rerun with the exact u4_cell hull (`p2\out\h1slab\pre_u4`); the points were rerun afterwards. Both enclose the same zero. To reproduce the reference bit for bit the kit ships those earlier enclosures (`reference\h1_pre_u4`, SHA-256 in `inputs\MANIFEST.sha256`) and `h1dt` / `win_slab` read them. `set KIT_DT_POINTS=final` uses this run's own points instead (then those files differ in the 9th digit).
* **Two aggregation scripts were inline in the proof tree** (`h1slab/summary.json`, `win_summary.json`): `p2\t9h\h1slab_summaries.py` reproduces them byte for byte.
* **The tail is run twice** (tail1, t7_nl, tail2) because T7's nonlinear closure reads `n3_b1`, which T8 computes from T7's constants (a fixed point, checked by `check.py`). `t9_cone_kch.json` is bootstrapped the same way (the values t9 reads do not depend on it; pass A and pass B must agree).
* `p2\t7` was a symlink to `t7` in our tree; the one script that used it (`p2\t8se\cconst.py`) now points to `step4\t7` directly.
* All source patches are marked `[kit patch]` (kit-relative paths through the `KIT_ROOT` marker as in the step-3 kit, `os.nice` is not used, `-I` is not used so that `PYTHONHASHSEED=0` applies). No fork, threads, signals or shell calls.
* Stages that cannot run on Windows: none. Not regenerated (C++/CAPD, listed below): the S1 profiles and the step-3 hull boxes.
* The first part of the chain (Gamma0, P1 slabs 0-96) is kit A; kit B trusts the shipped slab 95/96 checkpoints (hashes below).

## Send back

The SUMMARY block of `check.py` (or the whole output). If a stage fails, also the last 40 lines of `kit_out\logs\<stage>.log`.

## Expected runtime (estimate, not measured end to end)

The v3.2 part of this kit (kit B) has now been run end to end once (WSL, 76 PASS, 0 FAIL, `kit4b_run_wsl1/RUN_NOTES.md`); **this v3 build, with the six fixes and the 21 v3.3 stages, has not been run yet**. Kit B (2 cores): the full run is meant to be done on your machine as the independent reproduction. Individual stages and the whole tail were tested piecewise (fixtures, stage by stage, and `check.py` on the proof tree itself: all claims pass there). The sandbox stage timings add up to roughly 4.5 hours at 2 cores; the arb-128 polish (about 40 min) and the T7 certificate (about 75 min) parallelise little. On an i7-10700 (8 cores, 16 threads) with `--jobs 6` expect **about 2.5 to 4 hours**, depending on memory and thermals; the run is resumable (see above), and the console ETA is updated as stages finish. Use `--jobs 8` for a little more speed if the machine is otherwise idle.
