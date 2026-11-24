# Reproduction kit 4B: the rest of step 4 (base v2, T6, direct T6 on the corrected base, T7, T8, switch, (Y), slab-180 hand-off, Lemma O, v3.2 data)

Independent of kit A: it starts from our reference P1 checkpoints (slab 95 for the base seed, slab 96 for the switch), shipped in `step4\t3\p1\ckpt\`.
Everything is plain Python 3.13 + python-flint (arb), one single-threaded BLAS per process, no compiler, no network. Kit A and kit B may run at the same time
(kit B uses `--jobs 6` worker processes by default, kit A one core).

## Run it (Windows 11, PowerShell, unzipping to `C:\kit4b`)

```powershell
Expand-Archive -Path .\kit_step4b.zip -DestinationPath C:\kit4b
cd C:\kit4b
dir run_all.py                      # must exist here
py -3.13 -m venv kitenv
.\kitenv\Scripts\python.exe -m pip install --upgrade pip
.\kitenv\Scripts\python.exe -m pip install -r requirements.txt
.\kitenv\Scripts\python.exe run_all.py --jobs 6
.\kitenv\Scripts\python.exe check.py
```
`run_all.py` runs all 16 stages in order (up to `--jobs` worker processes at a time, default 6), then runs `check.py` itself. **Paste the check output back** (the SUMMARY block at least); no zip is needed.
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
| h1dt, winslab | slab-180 hand-off: 32 derivative balls, 32 window/shell balls, with their summaries |
| tail1 | the final tail (`reference\final_manifest.json`: switch, (Y), T8 cconst/merge, clam, camp / n3far iteration, window camp, t9 ...) on T7 without the nonlinear closure |
| t7_nl | T7 republished with the nonlinear closure (Lemma 4c) on tail1's `n3_b1` (the proof tree's T7 is circular with T8 in this way) |
| tail2 | the final tail again on the republished T7: this is the reference chain (34 stages incl. the 16 exact H1 points and H1 summary (round 1 only), t9_far, h4loc, Lemma O, v3.2 scripts, t1clock, axis crossing, and last the B8 check `check_outA5`). `check.py` verifies that round 1 and round 2 give the same `n3_b1` |

The console shows `[stage i/16]`, a heartbeat every 3 minutes and an ETA. Expected time: about 2.5 to 4 hours (see the end of this file).

## Interrupted? (Ctrl-C, reboot, crash)

Run the same command again. Finished stages (`kit_out\state\*.done`) are skipped; inside a stage every tool is restartable (the polish continues at the first unstored slab with `--resume`; T6, t6c, T7, rho_end, se_v2, h1, h1dt, win skip finished units). A second `run_all.py` in the same folder refuses to start.
The polish after a restart builds a fresh float Jacobian, so Newton paths after the restart slab differ from an uninterrupted run, within the tolerance |E| < 1e-26.

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

This kit has **not been run end to end in our sandbox** (2 cores): the full run is meant to be done on your machine as the independent reproduction. Individual stages and the whole tail were tested piecewise (fixtures, stage by stage, and `check.py` on the proof tree itself: all claims pass there). The sandbox stage timings add up to roughly 4.5 hours at 2 cores; the arb-128 polish (about 40 min) and the T7 certificate (about 75 min) parallelise little. On an i7-10700 (8 cores, 16 threads) with `--jobs 6` expect **about 2.5 to 4 hours**, depending on memory and thermals; the run is resumable (see above), and the console ETA is updated as stages finish. Use `--jobs 8` for a little more speed if the machine is otherwise idle.
