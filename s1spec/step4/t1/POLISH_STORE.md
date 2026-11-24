# POLISH_STORE: stored arb-128 polished base (P11, 188 slabs)

**Command** (1 process, `nice -n 10`, OMP/OPENBLAS/MKL threads = 1; same args as T1_REPORT plus `--store`; every 10th slab also gets norms):
`python3 -u polish_all.py P11 187 3,5,30,61,90,122,150,183,187 --store`
Wall 1938 s (32 min; `polish_P11.log`/`.json` regenerated, numerics unchanged). Storage: `polished/` = 55 MB (188 `slab_XXX.pkl` + `manifest.json`).

**Format** (`polstore.py`; the only code changes are that module and the `--store` block in `polish_all.py`). Each pkl: every arb of `Yn` ((s+1) x (n+4); rows = relative times 0, c_1..c_8; columns = DOF, then ℓ̃, ln λ, c, t) as exact ints (mid mantissa, mid exponent, rad mantissa, rad exponent); plus cfg, knot parameters `C` (all inserted elements), n, N, P, s, dt, tau0, k, prec = 128, polish record. `load_slab(path)` rebuilds Mesh/Rep/Dyn/Radau and a `SlabPolish` with `Yn` set; no re-polishing. ã, b̃ are derived from the state by `ArbOp` as before. Manifest: slab, tau0..tau1, k0..k1, dof, n_el, maxE, norms (25 sampled slabs). Note: the arb(mid, rad) constructor rounds radii up one mag-ulp; `dec` pre-compensates (verified on 20000 random arb).

**Reload check** (`verify_store.py`, fresh processes, slabs 3, 90, 187): re-encoded Yn identical to stored ints; the recomputed `slab_norms(SP, 6)` are bit-identical (all 9 norms, equal as doubles) to the in-run values for slabs 3, 90, 187. Max residual at the stored state (1e-28 to 6e-28) is below the in-run value (pre-final-step, up to 4.7e-27).

**Residual norms** (max over the 25 sampled slabs; targets from the report):
| | R2 | R1 | Rs | RIs | RF0 | RF1 | RF2 | RF3 | RF4 |
|---|---|---|---|---|---|---|---|---|---|
| target | 2e-14 | 5e-13 | 5e-11 | 1e-10 | 1e-7 | 1e-7 | 1e-7 | 1e-7 | 1e-5 |
| all sampled | 9.6e-6 | 2.0e-6 | 1.5e-6 | 8.6e-9 | 1.0e-3 | 4.0e-3 | 5.3e-2 | 1.1 | 28 |
| slabs >= 5 | 5.4e-16 | 6.3e-15 | 4.2e-14 | 2.9e-18 | 1.9e-11 | 4.7e-11 | 3.3e-10 | 5.4e-9 | 1.5e-7 |

**Misses**: sampled slabs 0 (k = 20.0, all nine norms) and 3 (k = 22.4: R2 = 2.0e-12, RF3 = 2.7e-7). Slab 5 (k = 24.2) onward and all later sampled slabs (10, 20, ..., 187) meet every target. Slabs 1, 2, 4 have no norms in this run (not sampled); slab 1 failed all but RIs in the report, so slabs 0-2 (k < ~22) and probably 4 should be treated as missing. Max polish residual E < 1e-26 everywhere.
