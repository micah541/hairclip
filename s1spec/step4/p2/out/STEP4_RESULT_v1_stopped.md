# STEP 4 RESULT: stopped at the switch, because design hypothesis T5 (I) fails

Status (Nov 10, 19:50 PT). P1 is finished. The switch (A), the P2 run (B) and T9 (C) were **not** done. The reason is in §1. Nothing was patched.

## 1. Finding: the stored base does not pass through P1's centre
- **What T5 assumes.** T5 §1.2 (I) says: "The base is started from that centre: ‖γ_P1(t_a) − (c̃+λ̃Γ̃)(τ_a)‖_C² ≤ 1e-14". The T7/T8 budget Ȳ ≈ 2.3e-12 is described as "dominated by P1's start value".
- **What was built.** T1's base (`t1/polished/`, the input of T6, T7 and T8) was seeded from the float HLS snapshot `feas2/snaps/p2_k20.npz` (T1_REPORT §0, "Seed"; risk 1). That snapshot is the feas2 run with tail_tol 1e-10, plus symmetrisation and a closure bump. The base is therefore a different trajectory from the true one: its blow-up clock is offset by Δt = 7.0e-6.
- **Measured at the switch.** P1 slab 95 (k = 26.2537786) against base slab 7, at the best clock τ* = 1.0727:

| quantity (rescaled, T5 §1.1) | value | against |
|---|---|---|
| A = ‖χe‖_ρ | **2.76e-9** | |
| B = ‖(χe)″‖_ρ, V₀ = ‖ζ∞e′‖_W∞ | 1.18e-9, 1.33e-9 | |
| **Y(τ0)** (μ_e = 0.036, ε = 1e-5) | **2.77e-9** | T8 requirement 1.60e-11: **×173 over**; 2.15e-11 (final C_γ0): ×129; T5/T7 budget 2e-12: ×1390 |
| sup\|e\| on x ≤ 34.25 / x > 34.25; sup\|δκ̃\| | 7.4e-7 / 1.9e-4; 2.0e-7 | |
| trend: A at k = 21.0 / 23.48 / 26.25 | 4.40e-9 / 3.51e-9 / 2.76e-9 | decays at ≈ 0.53 per unit τ |

- **How the numbers were obtained.** This is a float computation, first order in e, and it is *not* certified. Both curves are tip-normalised (κ(0) = K0, tip at −D_S1). The clock τ* and the p̃, ñ coefficients minimise ‖χe‖_ρ. To first order this is the minimum over all modulations (scale, axial shift, clock), so T5's constraint choice can only give a larger A.
- **Cross-checks.** e built from positions and e built from the ODE e′ = δθ − κ̃w, w′ = κ̃e agree to 4 digits. The three P1 checkpoints give a smooth trend.
- **Consequence for P2.** By (Y), Ȳ(τ) ≥ Y(τ0)e^{−ϱ(τ−τ0)} with ϱ ≤ 0.053. That gives Ȳ ≥ 1.4e-9 at τ ≈ 14.3 (k ≈ 850, when the window material crosses the ramp) and Ȳ ≥ 6.6e-10 at τ1.
  - Items 3 and 4 then fail by about ×50–100.
  - Item 1 loses its margin: its error share is about 4.9e-8 against 4.96e-8.
  - P2 cannot close with this base. The P1 ball (C² = 1.13e-13, ×3.8 over T5's 3e-14) contributes ≲ 1e-11, which is 300× less than the base mismatch.
- **True vs certified decay.** The true mismatch decays at about 0.53 per unit τ, so the true error is harmless by τ ≈ 10. The certified rate is about 0.04–0.05, and that is what (Y) can use.

## 2. (SE): a rotation does not settle it; how to settle it
- **The defect is a tilt, not an off-axis far point.** T6's turning defect is δ = 1.6506e-14. It is constant on slabs 5–187: ∫κ̃ is conserved exactly by the Leg₀ test function, so δ was inherited from the seed. The far point lies on the axis to float accuracy: |Im Γ̃(ℓ̃)| ≲ 5e-13 at slab 5, against δ|Re Γ̃(ℓ̃)| = 4.9e-12, with Re Γ̃(ℓ̃) = −0.76ℓ̃. Hence q̃(ℓ̃) = δ Re Γ̃(ℓ̃), which agrees with T6's bound q̃/ℓ̃ ≤ 1.31e-14 = 0.79δ.
- **Why a rotation fails.** The only rotations that are symmetries of the rescaled flow are those about the origin. R_{−δ} makes the far tangent vertical, but:
  - it moves the far point off the axis by |q̃(ℓ̃)| (up to 4.6e-9);
  - it tilts the tip by δ and moves it off the axis by D_S1δ = 4.1e-14, so tip parity is lost.

  That is two conditions with one angle, so exact (SE) cannot be reached by rotating.
- **Settlement, for the rebuilt base.** The comparison curve must be exactly σ-symmetric:
  - (i) Seed with turning exactly π in arb. P1's sine class has exact turning, and the scheme then keeps δ ≡ 0.
  - (ii) Remove the residual-driven far-point drift h(τ) = Im Γ̃(ℓ̃) with a physically static rigid shift of the far part. It is blended in physical arclength over a fixed macroscopic window (e.g. s ∈ [6, 12]).

  Its extra residual is the physical CDF residual of a static deformation of size |h|λ. That is an explicit term which T6 must certify.
- **On the old base.** The same construction (a blended rotation by −δ about the far point, physical amplitude ≤ 0.76δL/2 ≈ 3e-13) would give exact (SE) there too. This is a proposal: it is not implemented and not certified.

## 3. Ways forward (coordinator's decision)
- **A. Recommended: the design's own prescription.** Re-seed T1 from P1's exact centre: slab 93 (k = 21.0, Δτ ≈ 0.9 before the switch, covering the start-up layer) or slab 95. Then re-polish and re-run T6, T7 and T8. P1 is not rerun.
  - Expected Y(τ0): the P1-ball level, ≈ 5–10e-12, to be measured. That leaves a margin of about ×2–3 against 1.6e-11, because C² is ×3.8 over T5's assumption.
  - Cost: T1 about 0.7 h, T6 0.9 CPU-h, T7 about 4.5 CPU-h here plus nl7 3.3 CPU-h off-site, T8 3 min. The P2 driver and T9 are cheap.
  - φ1 and φ2 (`t7/phi_store.json`) may stay as fixed functions; T9 must use the same ones.
- **B. Keep the base.** This needs a certified decay rate of at least 0.39 on [τ0, 14.3], e.g. T5 fallback (c): a Lyapunov matrix on the first N modes. The measured true rate is 0.53, so the margin is thin. The frozen far mismatch (sup|e| ≈ 1.9e-4 at x > 34) must also be taken into F_j(τ0). This is a new T7 certificate type and is risky.

## 4. Chain status, end to end (certified unless marked)

| link | result | margin / status |
|---|---|---|
| Γ0 (T2, `t2/gamma0_final.json`) | M = 3762 sine modes; b₂ closed by Krawczyk (half-width 8.14e-24, contraction 3.0e-12); embedded, slot half-gap ≥ 0.0248 | certified |
| P1 (T3, 95 slabs, 5.5 CPU-h) | t = 4.597846334999414998 ± 6.9e-19 since Γ0; k_max ∈ [26.253778600561, 26.253778600564]; L = 48.32562154110466 ± 4.0e-15; X0 = 0.84121520310116 ± 1.5e-15 | C² 1.13e-13 vs 3e-14: **×3.8 over**; C⁵ lead 3.52e-10 vs 1e-9: ×2.8 |
| base (T1) | residuals certified (T6, slabs 5–187, margins ≥ 17) | **wrong trajectory, (I) fails (§1)** |
| T7 | min c_I = 0.0528 vs 0.045 (×1.17); λ < 2μ1 by 18.7% | provisional: c2r3 is done and constants.json was re-merged at 19:43, but not re-read here; μ_e, ϱ_c, K_e are placeholders (per T8); nl7 is running |
| T8 | requirement Ȳ ≤ 1.60e-11 (2.15e-11) | open: N3 far box; F_j(τ0) |
| switch / P2 / T9 | not run | blocked |

**Files and commands** (float probes, 1 process each at nice 10, under 1 min each):
```
cd /home/claude/s1spec/step4/p2
python3 -I f_err2.py ../t3/p1/ckpt/slab_00095.json 7   # switch: A, B, V0, Y, clock, sup|e|
python3 -I f_err.py  ../t3/p1/ckpt/slab_00095.json 7   # cross-check of A from positions
python3 -I f_err2.py ../t3/p1/ckpt/slab_00094.json 4   # trend (also 00093 with base slab 1)
python3 -I f_end.py 5,187                              # far-end geometry for (SE)
```

## 5. Draft statement: the step-4 theorem
**Theorem 4.** Let Γ0 be the closed curve X(α) = X0 + (L0/2π)∫₀^α e^{iθ}, with θ = −π/2 + α + Σ_{n≤3762} b_n sin nα, and with the exact data of `t2/gamma0_final.json`:
- b_n exact doubles for n ≠ 2;
- b₂ the unique closing value in its certified interval;
- L0 = 0x1.9e9c2d62f0874p+6 and X0 = (x0, 0).

Γ0 is smooth, embedded and symmetric about the e₁-axis. Let γ(t) be the curve diffusion flow with γ(0) = Γ0. Then there are t₁ > t_a = 4.597846334999… (time since Γ0), λ₁ > 0 and c₁ ∈ ℝ such that:
- (a) γ exists on [0, t₁] as a smooth, axis-symmetric solution, and k_max(γ(t₁)) ∈ [2.03e4, 2.36e4];
- (b) Γ₁ = (γ(t₁) − c₁e₁)/λ₁ satisfies the step-3 hand-off hypotheses at τ₁ (v3.1; functionals as in lemmaS.pdf p. 1, lemmaO_v3.tex Revision and step3_proof.tex):
  - (1) |α| ≤ 1.5e-4 and √Λ ≤ 5e-8;
  - (2) shells on [14, 60] have (E0, E1, E2) ≤ (2e-7, 1e-4, 0.2);
  - (3) windows [y, 1.5y], y ∈ [50, 253], ω-normalised, with v = f − 1.3056αy^1.0246: E1–E3 ≤ 1e-5 and E0 ≤ 2e-7;
  - (4) K4 ≤ 0.06;
  - (5) far-cylinder jets in λ₁ units: |h″| ≤ 8.8e-3, |h‴| ≤ 7e-5, |h⁗| ≤ 6.2e-7, slope ≤ 0.21.

Hence, by step 3, γ develops a finite-time singularity modelled on S1.

Proof chain: T2 → P1 (validated, [0, t_a]) → switch (IFT for λ₀, c₀, τ₀) → Theorem T5 on [τ₀, τ₁] with (R) from T6, (C) from T7 and (F) from T8 → T9. **Open: the switch, P2 and T9, blocked by §1.**
