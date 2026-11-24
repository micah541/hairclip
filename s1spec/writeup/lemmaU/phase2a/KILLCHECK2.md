# Kill check 2, Lemma U — independent review (Nov 13, 2026)

**Verdict: GO WITH CONDITIONS.** The rule passes with room: late-centre M1 is ×2.39 / ×2.28 / ×2.10 at KS = 1 / 1.33 / 2 once A's export is consumed correctly (U_model2_A1's ×2.57 / ×2.36 is not the honest number), and E₂ᴬ = 0.0696 [cert]. The horizon problem is a mismatch between two statements, not a defect of the method. The one structural item is the outer-cutoff hypothesis H-O2: as A wrote it, no lemma can supply it; it must be redone in L² form before anything else is built on (A).

Labels: **[cert]** arb, project script. **[est]** project float model on certified inputs. **[rev]** my scratch run (`killcheck2_scratch/kc2_model.py` → `kc2_model.out`: U_model2's unmodified source with the data term changed as in "Horizon" below; its self-check passes and row R0 reproduces A1 exactly).

## Margins (weight (4,4,0.3), data56, box (0.15, 0.52), late centre)

| quantity | value | label |
|---|---|---|
| E₂ᴬ repaired, tightened (H2) | 0.0696 (rule ≲ 0.1) | [cert], given H-O placeholders (they contribute 1e-6) |
| (A) core M1, tightened / v3.1 (H2) | ×4.36 / ×2.73 | [cert], same condition |
| Rate margin, Ĝ₂₂ = 338259 | ×29 | [cert] |
| U_model2_A1 as run (D from τ₁, slope ignored), KS 1 / 2 | M1 ×2.57 / ×2.36 | [est], **not valid for late centres** |
| **Honest: window-uniform data C0 + C1·T**, KS 1 / 1.33 / 2 / 4 | **M1 ×2.39 / ×2.28 / ×2.10 / ×1.75**; M0 ×4.2 | [rev] |
| Same, present centres, KS 1 / 2 | M1 ×2.84 / ×2.39; M0 ×2.09 | [rev] |
| H-O placeholders ×10, KS 1 / 2 | M1 ×2.32 / ×2.04; M0 ×3.3 | [rev] |
| E₂ᴬ at the rule limit 0.1, KS 1 / 2 | M1 ×2.26 / ×2.00 | [rev] |
| Fallback box (0.20, 0.60), KS 1 / 2 | M1 ×2.36 / ×2.07 (K not certified there) | [rev] |
| Joint stress: H-O ×10, KS 2, C0 ×1.5, E₂ᴬ = 0.1 | M1 ×1.81 (×1.78 with the fallback box) | [rev] |

The certified effective KS is ≤ 1.33 (P2), so the governing row is ×2.28. Hashes in `A_repair/out/SHA256SUMS` verify; I did not rerun the chain.

## 1. Does U_model2 consume A's export correctly? No, in three places; none changes the verdict.

- **D_slope is ignored and the wrong constant is read.** `U_model2.py:98` takes `D` = (1.40e-7, 7.59e-5), which is the integral from τ₁. For a centre launched at τ_ρ ≫ τ₁ that bound carries D_slope·(τ_ρ − τ₁), unbounded. The correct datum is A's window statement (`A_REPAIR.md` R4, `assembly_rep.json` export): C0 = (2.925e-7, 1.510e-4), C1 = (5.06e-9, 5.85e-9), valid from every τ_a. D₁ doubles; M1 falls 7–11%.
- **D₂ = 71 is never read.** Only D₀ and D₁ enter (z₂, z₃ data); K has no z₄ entry (asserted in my script). Harmless.
- **The SENS rows of the A1 run are not live.**
  - `U_model2.py:307`: with `--inputs`, all three "D_k export" rows are `dict(inp0)`, i.e. the headline three times.
  - `U_model2.py:310`: the E₂ᴬ rows use `inputs_default`, so they are the default run's rows, not variations of A's export.
  - My R2–R5 rows replace them.
- Not certified for the headline: K is certified at the crude box (0.11, 0.52) and g ∈ [3e5, 6e5] (`phase1/ramp/ramp_K_cert.json` note). The headline box is (0.15, 0.52); the fallback has g = 6.29e5. KS = 2 covers it; rerun P2.

## 2. Horizon: fixable by restating Lemma D. Cost: the 7–11% above.

- **Not an artefact.** (A) lives on a fixed rescaled region. Its far part [350, 700] always holds O(0.1) of solution, so forcing, coupling and outer terms are per-unit-τ sources for ever. A horizon-free D cannot exist. `U_ramp.tex:35,148,158` (∫_{τ₁}^∞ f⁺ finite) is false as stated.
- **Not fatal.** Lemma U never needs the total. Its ramp weight for a centre that has passed y = 56 decays like e^{−(τ−τ_ρ)} (μ = 4, c = ¼), and is monotone in τ. With X(T) ≤ C0 + C1·T on every window, summation by parts gives

  ∬ Ψ εz² ≤ S_X·C0 + C1·Σ_X, with S_X = 1.52 and Σ_X = 3.0 [rev].

  The slope term is 1e-4 of the constant at level 1, and 3% with H-O ×10.
- **What the theorem needs.** Physical time is finite, and at fixed physical radius the far field barely moves (Lemma O, time ratio 1.6e-7). τ → ∞ matters only for sup-in-τ bounds, and those are damped: E^A takes rate/γ′, Lemma U takes the decaying weight. No bound has to be integrable to τ = ∞.
- **Check needed (1 h, mechanical):** the verifier's J_m(T) ≤ J_m + s_m·T must hold on every window, not only from τ₁. v3.1 already uses it that way (`lemmaS.tex:256`, layer-cake). Confirm in `unconditional()`.

## 3. H-O1 / H-O2

- **H-O1 (box ½, 1 on [350, 700]):** fine. Source: B3U box where r < 190λ₁, D5a where r ≥ 190λ₁.
- **H-O2 pointwise (|z₂|, |z₃| ≤ 1) has no source.** f is the slope, so z₃ is a fourth derivative of the curve: Lemma O does not bound it (`step3_proof_v32.tex:265`), and Lemma U carries levels 0–2 only, which give sup z₀, z₁. `outer_rep.py:97` uses z₂ pointwise cubed at y = 350, and `:89` bounds z₃ pointwise at level 2.
- **The plan already had the right form:** `launch_lemma.tex:119` (D3′) takes z₀,₁,₂ from B3U energies in L², z₃ by Young. A's R5 departed from it. Redo R5 that way, as on the inner ramp.
- **Sensitivity is negligible.** The weight at 350 is 5e-12. E₂ᴬ adds: 1e-6 at placeholders, 8e-4 at M̄₂ = 10, 8.5e-3 at M̄₃ = 1000, 7.5e-3 at all ×10 (E₂ᴬ = 0.077). E₀ is the sensitive level (×14 at all ×10); M0 stays ×3.3. T3's L² estimate was +4e-5 on E₂ᴬ [est].
- **Circularity:** a bootstrap loop, not a logical circle. (A) uses U at hypothesis level (2B̂); U uses the derived E^A; the loop gain through the outer zone is under 1e-3. It needs the open–closed argument written out with H-O as a consequence of B3U + D5a.

## 4. New ingredients

- **φ = exp(K_b e^{−s}):** sound. It is bounded by e^{K_b/28}, co-moving decay beats the b-rate by ×71 or more. Γ_m and K need a rerun with φ (≤ 1e-4 [est]).
- **data56, global box (0.15, 0.52):** sound; rate ×29 [cert]. Zone III at hand-off is about 0.09 and D5a moves it by √(1.6e-7) terms, so (0.15, 0.52) is plausible (×1.6) but not certified. The fallback (0.20, 0.60) costs 1% of M1 [rev] and a P2 rerun.
- **Cutoff × jet in Lemma U's own ramp:** the headline route has no integration by parts, so A's ×8 mechanism (IBP-derived jets) does not occur; P2 enumerated the direct classes. Two cautions:
  - the Lemma Y variant (×2.67) does create them, unbounded. Do not use Lemma Y numbers;
  - rerun A's class enumerator (`ramp_rep.py`) on [28, 56] as a cross-check of P2. S33 / foot is 2e-4 of B₂, so ×100 room.
- The earlier reviewer estimate "≲ 1.5e-3" for (A)'s cutoff × jet was wrong by ×8 (mine; A's 0.0118 is certified).

## 5. v3.1 / v3.2 status

Closed step 3 is **conditional on more than (F)** until v3.3 exists:
1. Lemma S(A)'s proof as written is incomplete (top pair, cutoff × jet, U_j, I₃). The repaired statement holds at the v3.1 (H2) with E₂ = 0.2686, core M1 ×2.73, FLUX ×4.9 / ×5.3 [cert]. It has different constants (E, γ′, E_B, K0, K2, I₃) and is not yet in a versioned write-up.
2. **Hypothesis (O), new:** the outer cutoff of (A) on [350, 700]. v3.1 never treated it. Its inputs sit in the far field, which is the (F) region.
3. (F), as before.

`main.tex` should say "conditional on (F) and (O); Lemma S(A) holds with corrected constants (erratum in preparation)". It should not say step 3 is closed except for (F). No numerical threat: every repaired margin is ≥ ×2.

## 6. Recommendation: GO to the bulk, in this order

| # | task | who | hours |
|---|---|---|---|
| 1 | Fix the three statements: Lemma D as C0 + C1·T with summation by parts; R5 in L² form (D3′); bootstrap order with H-O derived. Decide box (0.15, 0.52) | Fable | 2–3 |
| 2 | (A): R5 redone and certified; export C0 / C1 as the interface; H-J window check; interface conventions at y = 28 audited | Opus | 6–8 |
| 3 | Lemma U proof text (Prop. U, T with φ, C, D; no Lemma Y); compile all tex off this machine | Opus | 12–15 |
| 4 | Arb: Σ sweeps, S_X, launch constants, Φ_j, c₂₂₃; K and Ĝ at the widened box with φ; class-enumerator cross-check | Opus, Sonnet for reruns | 10–12 |
| 5 | Zone III box from D5a in the S₁ chart; Lemma R ramp (the typed 0.0161); cone family, fixed L | Opus | 8–10 |
| 6 | Step-4 hand-off reruns at the tightened (H2), manifest, kits | Sonnet | 3–5 |
| 7 | v3.3 assembly, v3.1 erratum, main.tex / STEP4_RESULT | Sonnet, Opus for the erratum | 6–8 |
| 8 | Independent review | Fable | 3–4 |

**Total: 50–62 agent-hours, including the 3–5 of step-4 reruns** (was 45–55 + 3–5).

**Stop and report if:** the certified late M1 after task 4 is below ×1.5; D5a cannot give zone III ≤ (0.20, 0.60); R5 in L² form adds more than 0.01 to E₂ᴬ; or H-J fails on windows.

**Needs Micah's decision:** the wording change in `main.tex` (section 5), and going ahead with the bulk.
