# lean/ : vendored Lean package (subset of CWFlow)

**Author decision needed before anything here is made public.** The 16 hair-clip files
(`CWFlow/Curve/*.lean`) import a large part of the general CWFlow library. Their transitive
import closure is **184 of the library's 198 `.lean` files (44,529 lines)**: almost the whole
library, including `Extension/`, `Flow/`, `Interp/`, `Moser/`, `Weinstein/`, `Bridge/` and `Scalar/`.
Only `CWFlow/AK/` (8 files), 1 `Extension/`, 2 `Flow/`, 2 `Scalar/` and 1 `Weinstein/` file are
left out. If this repository is published, all of those library files become public. The
alternatives are: publish CWFlow on its own first and require it here as a git dependency
(then only `CWFlow/Curve/` would live here), or keep `lean/` out of the public repo.

## Source
- Repository: the author's local CWFlow repo (no public remote recorded), commit
  `3130aef39829f02c126aeb166fbb9583f7670108` ("Curve: n=1 CDF, L2 blow-up criterion,
  maximal flow, shift/glue"). The files were taken with `git show <commit>:<path>`, so the
  bytes are exactly those of that commit.
- The 16 `CWFlow/Curve/*.lean` files are byte-identical to the snapshots
  `s1spec/writeup/lean_bridge/snap_final/` and `snap_s3/` cited by the step-3 texts.
- The closure was computed by parsing the `import` lines, starting from `CWFlow.Curve.All` and
  following only `CWFlow.*` modules (script: `hairclip_repo_build/tools/closure.py`, kept outside the repo).

## Dependencies (not vendored)
- `ShortTime`: the author's package `parabolic-short-time-lean`. The original lakefile uses a
  local path (`../MM`); `lakefile.toml` here requires it from git at commit
  `7b15b01b2d957ffa977396b3d555ad37b27cb515`. That commit was the clean HEAD of `../MM`, and it is on
  `origin/main`. Whether that GitHub repository is public was **not** checked. It pins
  DifferentialGeometry v0.1.3 (qinz1yang/differential-geometry) and Mathlib v4.33.1.
- Toolchain: `leanprover/lean4:v4.33.1` (`lean-toolchain`).
- The original `lakefile.toml` and `lake-manifest.json` are kept as `lakefile.orig.toml` and
  `lake-manifest.orig.json`. They are for reference only, because they point at `../MM`. A fresh `lake update`
  writes a new manifest.
- The library root is `CWFlow.Curve.All` (`roots` in `lakefile.toml`). The original `CWFlow.lean`
  (which imports the whole library) is not vendored.

## Status
- **This vendored package was NOT built here.** A Lean/Mathlib build is too heavy for this machine
  (about 3.5 to 11 GB per process). What was checked: every `import` in the 184 files resolves either to a
  vendored file or to a module of Mathlib, DifferentialGeometry or ShortTime (109 external imports;
  Mathlib and DifferentialGeometry were checked against the local `.lake/packages` copies, ShortTime against
  `../MM` at 7b15b01).
- The original repository compiled these files with clean axioms ([propext, Classical.choice,
  Quot.sound]). Logs are in `s1spec/writeup/lean_bridge/` (`raw_*.log`, `full_build_Curve_All_3130aef.log`).
  The full-build log is from the original repo, not from this vendored copy.
- To build: `lake exe cache get && lake build CWFlow.Curve.All` (use `LEAN_NUM_THREADS=1..2` if memory is low).

## Main results (CWFlow/Curve)
T1 `cdf_exists_maximal` (maximal CDF flow; if T_max < ∞ then ∫κ² ds is unbounded); T2 `cdf_extends_of_bounded_energy`;
T3 `aNormSq_eq_curv_sq`; T4 `IsFlow.lift_isCDF`; T5 `IsCDF.curv_sq_le`; `exists_isLagImmersion_of_periodic`,
`cdf_maximal_shift`, `isFlow_glue_endpoint`; n=1 port `cdf_exists`, `cdf_unique`, `cdf_unique_up_to_reparam`,
`cdf_symmetric`, `flowVel_eq_neg_kss` (N1.lean). See `s1spec/writeup/lean_bridge/REVIEW_W7.md`.

## Closure by directory
| Directory | Files | Lines |
|---|---:|---:|
| CWFlow/Bridge/ | 11 | 2332 |
| CWFlow/Curve/ | 16 | 3042 |
| CWFlow/Extension/ | 76 | 23517 |
| CWFlow/Flow/ | 30 | 5247 |
| CWFlow/Interp/ | 19 | 5424 |
| CWFlow/Moser/ | 7 | 1300 |
| CWFlow/Scalar/ | 15 | 2281 |
| CWFlow/Weinstein/ | 10 | 1386 |
| **total** | **184** | **44529** |

## Full file list (lines, path)
```
   404  CWFlow/Bridge/ChartJet.lean
   291  CWFlow/Bridge/Contract.lean
    52  CWFlow/Bridge/Data.lean
   168  CWFlow/Bridge/Defs.lean
   104  CWFlow/Bridge/Existence.lean
   163  CWFlow/Bridge/Fields.lean
   172  CWFlow/Bridge/GraphChart.lean
   375  CWFlow/Bridge/Realizable.lean
    91  CWFlow/Bridge/ScalarSolution.lean
   131  CWFlow/Bridge/ScalarTensor.lean
   381  CWFlow/Bridge/Smooth.lean
    19  CWFlow/Curve/All.lean
   502  CWFlow/Curve/Arclength.lean
   222  CWFlow/Curve/CircleChart.lean
   120  CWFlow/Curve/Criterion.lean
    41  CWFlow/Curve/Defs.lean
   384  CWFlow/Curve/Energy.lean
    41  CWFlow/Curve/Lift.lean
    88  CWFlow/Curve/Lift3.lean
   104  CWFlow/Curve/LiftA.lean
    54  CWFlow/Curve/LiftC.lean
    76  CWFlow/Curve/LiftV.lean
   171  CWFlow/Curve/MaxShift.lean
   341  CWFlow/Curve/N1.lean
   387  CWFlow/Curve/N1Basic.lean
   347  CWFlow/Curve/Periodic.lean
   145  CWFlow/Curve/Periodic2Circle.lean
   240  CWFlow/Extension/ChartAbsorb.lean
   249  CWFlow/Extension/ChartBounds.lean
   548  CWFlow/Extension/ChartDeriv.lean
   244  CWFlow/Extension/ChartEnergy.lean
   432  CWFlow/Extension/ChartGraphEvol.lean
   115  CWFlow/Extension/ChartSup.lean
   258  CWFlow/Extension/CovSymbols.lean
   465  CWFlow/Extension/Cover.lean
   158  CWFlow/Extension/EnergyBound.lean
   331  CWFlow/Extension/EnergyCompare.lean
   182  CWFlow/Extension/EnergySymbol.lean
   476  CWFlow/Extension/EvolPointwise.lean
   134  CWFlow/Extension/EvolWeights.lean
   264  CWFlow/Extension/FlowBridge.lean
   168  CWFlow/Extension/FlowClosed.lean
    41  CWFlow/Extension/FlowM.lean
   297  CWFlow/Extension/GNorm.lean
   219  CWFlow/Extension/GaugeTerm.lean
   440  CWFlow/Extension/GlobalIneq.lean
    83  CWFlow/Extension/GlobalSetup.lean
   170  CWFlow/Extension/GraphAt.lean
   364  CWFlow/Extension/GraphBridge.lean
    56  CWFlow/Extension/GraphChart.lean
   291  CWFlow/Extension/GraphCoeff.lean
   234  CWFlow/Extension/GraphLocal.lean
   417  CWFlow/Extension/GraphPiece.lean
   329  CWFlow/Extension/GraphRadius.lean
   348  CWFlow/Extension/GraphRadiusLift.lean
   349  CWFlow/Extension/GraphRadiusLocal.lean
   534  CWFlow/Extension/GraphTower.lean
   498  CWFlow/Extension/IBPBounds.lean
   128  CWFlow/Extension/ImmCoord.lean
   460  CWFlow/Extension/ImmDiv.lean
   261  CWFlow/Extension/ImmPot.lean
   455  CWFlow/Extension/ImmReparam.lean
   270  CWFlow/Extension/ImmSmooth.lean
   543  CWFlow/Extension/IntM.lean
   102  CWFlow/Extension/JetIBP.lean
   321  CWFlow/Extension/JetOfSymbol.lean
   157  CWFlow/Extension/JetSymm.lean
   192  CWFlow/Extension/JetTime.lean
   680  CWFlow/Extension/LeadingLower.lean
   703  CWFlow/Extension/LeadingMain.lean
   172  CWFlow/Extension/LimitChart.lean
   523  CWFlow/Extension/LimitChartAux.lean
   236  CWFlow/Extension/LimitM.lean
   213  CWFlow/Extension/LimitMAux.lean
   180  CWFlow/Extension/LimitScaled.lean
   170  CWFlow/Extension/LimitT.lean
   347  CWFlow/Extension/LocalEvolution.lean
   122  CWFlow/Extension/LocalWeight.lean
   850  CWFlow/Extension/Partition.lean
   455  CWFlow/Extension/Persist.lean
   386  CWFlow/Extension/PieceFlow.lean
   227  CWFlow/Extension/PotReparam.lean
   366  CWFlow/Extension/PullZ.lean
   252  CWFlow/Extension/Rescale.lean
   184  CWFlow/Extension/SliceSmooth.lean
   458  CWFlow/Extension/Slope.lean
   203  CWFlow/Extension/SymbolDeriv.lean
   172  CWFlow/Extension/SymbolLow.lean
   378  CWFlow/Extension/SymbolProd.lean
    60  CWFlow/Extension/SymbolTime.lean
   234  CWFlow/Extension/SymbolTotalD.lean
   493  CWFlow/Extension/Symbols.lean
   183  CWFlow/Extension/TensorOps.lean
   186  CWFlow/Extension/TimeGraph.lean
   527  CWFlow/Extension/TimeGraphFlow.lean
   242  CWFlow/Extension/TimeGraphUnif.lean
   465  CWFlow/Extension/VelBound.lean
   304  CWFlow/Extension/VolChart.lean
   498  CWFlow/Extension/VolEvol.lean
   367  CWFlow/Extension/VolM.lean
   258  CWFlow/Extension/WeightIntegral.lean
   297  CWFlow/Extension/WeightTop.lean
   303  CWFlow/Extension/Weights.lean
   167  CWFlow/Flow/DiffeoChart.lean
    50  CWFlow/Flow/ExtensionTheorem.lean
   240  CWFlow/Flow/FlowExistence.lean
    56  CWFlow/Flow/FlowS.lean
   117  CWFlow/Flow/FlowUnique.lean
   227  CWFlow/Flow/Gauge.lean
   250  CWFlow/Flow/GaugeEvolution.lean
   163  CWFlow/Flow/GaugeGraph.lean
   235  CWFlow/Flow/GaugePotential.lean
   271  CWFlow/Flow/GlobalGauge.lean
   148  CWFlow/Flow/GraphImm.lean
   468  CWFlow/Flow/GraphUnique.lean
   106  CWFlow/Flow/InducedMetric.lean
    94  CWFlow/Flow/JointFamily.lean
    98  CWFlow/Flow/LiftTopology.lean
    53  CWFlow/Flow/LocalInverse.lean
    44  CWFlow/Flow/Main.lean
   117  CWFlow/Flow/MaximalFlow.lean
   323  CWFlow/Flow/NormalFlow.lean
   295  CWFlow/Flow/NormalTransport.lean
   199  CWFlow/Flow/Potential.lean
   122  CWFlow/Flow/Reparam.lean
   315  CWFlow/Flow/Restart.lean
   144  CWFlow/Flow/RestartScalar.lean
   273  CWFlow/Flow/ShortUnique.lean
    89  CWFlow/Flow/TanAlgebra.lean
   160  CWFlow/Flow/TanField.lean
   208  CWFlow/Flow/Tangential.lean
    41  CWFlow/Flow/UniqueMain.lean
   174  CWFlow/Flow/Velocity.lean
   196  CWFlow/Interp/BallHamilton.lean
   347  CWFlow/Interp/BallInterp.lean
   299  CWFlow/Interp/BallPotential.lean
   143  CWFlow/Interp/ClosureExt.lean
   171  CWFlow/Interp/EnergyGronwall.lean
   357  CWFlow/Interp/Hamilton.lean
   625  CWFlow/Interp/HamiltonStep.lean
   224  CWFlow/Interp/Jets.lean
    99  CWFlow/Interp/LogConvex.lean
   432  CWFlow/Interp/OdeJet.lean
   623  CWFlow/Interp/OdeTower.lean
   318  CWFlow/Interp/ParamInverse.lean
   262  CWFlow/Interp/RadialLift.lean
   219  CWFlow/Interp/RayLift.lean
   192  CWFlow/Interp/RealLemmas.lean
    92  CWFlow/Interp/SupBall.lean
   287  CWFlow/Interp/SupBound.lean
   326  CWFlow/Interp/ThreeOrders.lean
   212  CWFlow/Interp/TowerBound.lean
   105  CWFlow/Moser/Chart.lean
   137  CWFlow/Moser/Field.lean
   231  CWFlow/Moser/Flow.lean
   194  CWFlow/Moser/Forms.lean
   165  CWFlow/Moser/Natural.lean
   344  CWFlow/Moser/Theorem.lean
   124  CWFlow/Moser/Transfer.lean
   140  CWFlow/Scalar/Flat.lean
   163  CWFlow/Scalar/Graph.lean
   187  CWFlow/Scalar/GraphTime.lean
   125  CWFlow/Scalar/Invariance.lean
   222  CWFlow/Scalar/JetFunction.lean
   305  CWFlow/Scalar/JetRec.lean
   134  CWFlow/Scalar/LinearChange.lean
    57  CWFlow/Scalar/Multilinear.lean
   195  CWFlow/Scalar/Pointwise.lean
    58  CWFlow/Scalar/PotFamily.lean
   349  CWFlow/Scalar/Principal.lean
   136  CWFlow/Scalar/Reparam.lean
    64  CWFlow/Scalar/Smooth.lean
    93  CWFlow/Scalar/Symm4.lean
    53  CWFlow/Scalar/TaylorComp.lean
    60  CWFlow/Weinstein/Bridge.lean
    85  CWFlow/Weinstein/Cn.lean
   238  CWFlow/Weinstein/Coord.lean
   132  CWFlow/Weinstein/Existence.lean
   238  CWFlow/Weinstein/Global.lean
   235  CWFlow/Weinstein/Glue.lean
    95  CWFlow/Weinstein/Interface.lean
    96  CWFlow/Weinstein/Lift.lean
    78  CWFlow/Weinstein/Overlap.lean
   129  CWFlow/Weinstein/ZeroSection.lean
```
