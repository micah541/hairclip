# Kit B first end-to-end run (WSL, Nov 14 2026)
Shipped zip kit_step4b.zip (sha256 18a41a7c205d344e...), unzipped fresh here. `nice -n 19 python3 run_all.py --jobs 2`, Python 3.14.4, flint 0.9.0, numpy 2.5.3.

## Kit bug 1 (t7_b, stage 11): missing step4/t7/se_widths_secorr.json
- Symptom: `se7.py widths` (source t6c) needs `t7/se_widths_secorr.json` for slabs 0-1 (outside P2); the kit neither ships nor makes it. All t7_b steps then fail.
- In the proof tree it is the pre-B1 se_corr-source widths file kept under that name (T7_REPORT.md:318).
- Fix applied in this run dir only: `PYTHONHASHSEED=0 python3 se7.py widths --source se_corr --p2 <kit>/step4/p2`, then `mv se_widths.json se_widths_secorr.json`; removed the failed `step4/t7/run_v2se`; resumed run_all.
- The regenerated file differs from the proof tree's at relative ~1e-13 (this run's `p2/out/se_corr.json` already differs from the proof tree's, the expected platform/polish-path variance); slabs 0-1 values here are slightly smaller.
- Builder fix needed (kit_step4b_build): in `st_t7_b`, run the se_corr widths step and rename before `run_t7` (or ship the file as a pinned input). Not yet changed.

## Kit bug 2 (h1dt, stage 12): the 16 exact H1 points are never computed
- `run_all.py` defines `st_h1pts` (h1.py, 16 points -> p2/out/t9_h1_p*_coh_s180.json) but it is missing from STAGES; `h1slab_summaries.py h1` (end of h1dt) needs those files -> FileNotFoundError. The 32 dt balls themselves finished.
- Fix in this run dir: inserted ('h1pts', st_h1pts, 20) before h1dt in STAGES (backup run_all_pre_fix2.py); resumed. Builder fix needed likewise.

## Monitoring mistake (mine)
- The run stopped at 22:05 Nov 14; my watcher used `pgrep -f "run_all.py --jobs 2"`, which matched its own command line, so no alert fired; ~13 h idle until Nov 15 11:20.

## Kit bug 3 (tail1, switch): doubled path when the kit sits under /home/claude/s1spec
- `run_all.py` `tr()` rewrites '/home/claude/s1spec/step4' -> S4 and then '/home/claude/s1spec' -> KIT with chained str.replace; with KIT = /home/claude/s1spec/kit4b_run_wsl1 the second replace hits the first result (.../kit4b_run_wsl1/kit4b_run_wsl1/step4/...). Only affects a kit unpacked under /home/claude/s1spec (not Windows C:\kit4b), but it is a real bug.
- Fix in this run dir: single-pass re.sub (backup run_all_pre_fix3.py). tail1's cconst/merge steps ran before the failing switch; no stray directory was created.

## Kit bug 4 (tail1, n3boxes [A]): round 1 cannot pass the tip-box assert
- n3boxes asserts T7 meta.tip_box (Lemma 4b, nl closure), which only exists after t7_nl republishes T7; tail1 runs before t7_nl, so the assert always fails there. The n3 box inequalities themselves held in round 1 with the proof tree's margins (x6.43, x4.69, x1.99, x1.70).
- Fix in this run dir (bootstrap pass and pass A): tail1 (round 1) skips the n3boxes stages (pure checks, no outputs); tail2 runs them in full (backup run_all_pre_fix4.py).

## Kit bug 5 (tail1, t9 [A]): round 1 is circular as designed
- After fix 4, t9.py [A] fails too: KeyError 'tip_box' in T7 constants (t9.py:102). So round 1 cannot run any group-A stage beyond the T8 ones.
- Fix in this run dir (supersedes fix 4 for round 1): tail1 = pre stages + the 'T8 ...' stages of pass A (they produce n3_b1.json, all t7_nl needs), then return; tail2 runs the full tail incl. the t9 bootstrap, passes A/B, n3boxes, post (backup run_all_pre_fix5.py). t9_cone_kch.json does not exist yet, so tail2 bootstraps it.

## Kit bug 6: B8 stage (st_b8, check_outA5) never runs
- Like st_h1pts, `st_b8` is defined but missing from STAGES, so check.py WARNs "reference outputs not produced: check_outA5.json". Ran `check_outA5.py` by hand in kit/s1_tight/consistency (9 s): ALL INSIDE, output compared with the proof tree below.

## Result (Nov 15 12:45, after fixes 1-5)
- check.py: 76 PASS, 0 FAIL, 2 WARN; every claim equals the proof tree's value (e.g. rho_c 0.04711 vs 0.04713 is the only visible digit change; H1 sqrtLam 2.284e-8, |alpha| 1.174e-4, H5 x9.76, H4 x4.18).
- WARN "round 1 vs round 2 n3_b1 differ": false alarm, the files differ only in `args` (sha compare); all numbers identical (checked field by field: 0 differences). The T7/T8 fixed point holds.
- outputs: 30 byte-identical, 467 numerically different in low digits (platform/polish path, expected per README), 1 missing (B8, now run by hand).
- B8 output: identical to the proof tree's apart from the outA5 path and timing; check.py still reports it missing because it resolves '../../s1_tight' relative to step4 (one level off when the kit is not at the drive root?) -> checker path bug, minor.
