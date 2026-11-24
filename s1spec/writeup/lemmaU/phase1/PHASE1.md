# Lemma U, Phase 1 (approved by Micah, Nov 13 2026: "Start Phase 1 only")

Source of tasks: `writeup/step3_gap_v32/gonogo/GONOGO.md` (risks 1–3; staging) and `HANDOFF.md` there.
Weight: (mu, nu, a) = (4, 4, 0.3), a in [0, 0.3] kept free; inner ramp [28, 56]; (H2) = (1e-8, 1e-6, 1e-3), with (H3) included.

| dir | task | owner |
|---|---|---|
| ghat/ | P1: certified (arb) Gamma_m, kappa~, kappa2~ over the geometry; certified Ghat_jj with general-weight L4 classes; fixed point in g; also on the radius-dependent box | agent A |
| ramp/ | P2: honest general-weight ramp constants K_jk on [28, 56] in arb (|coef| tables), L4 ramp foot terms with C_eta^Xi | agent B |
| launch/ | P3: launch lemma (paper) + T3 fix/rerun | agent C |
| (root) | P4: T1 rerun with certified inputs -> checkpoint (continue if M1 margin >= ~x1.3); Fable checkpoint review | coordinator |

Rules: CLAUDE.md conventions. Never edit writeup/lemmaS*, step3_proof*, step4/. New files only. Every number labelled [cert]/[est]/[paper].

## Status
- **P1 done** (ghat/P1_REPORT.md). All [cert]:
  - (4,4,a), a ∈ [0, 0.3]: Γ_m = 4^m, κ̃ = 8, κ̃₂ = 64.
  - Crude box: Ĝ₂₂ = 239750.34, rate margin ≥ 2.5637, Gw ≤ 1.6395. The crude box governs y = 28 for both families.
  - Radius-dependent box: Gw ≤ 1.397 [cert arithmetic], but it needs a localised Ĝ derivation that does not exist yet.
  - Umax must be the support's value (0.09): rate margin 2.45.
  - Open: the ∂_τΘ transport term of the moving weight; no sign shown yet.

## FINDING in closed v3.1 Lemma S(A) (needs Micah's decision; nothing edited)
- The certified run (lemmaS_F10p_M0.52.json args, capture_F2/assembly_arb) has `'s33': 'assumed'` on the core: sup|z₂| ≤ M₂* = 3 (weights_exact.py:24, l.172/246 default).
- lemmaS.tex l.19/l.349 states the core {2,3,3} is "closed by the L⁴ identity and IBP", and l.91 calls M₂* = 3 superseded. The ramp uses ibp2 correctly.
- P1 [cert arithmetic]: the IBP route gives core Ĝ₂₂ 6.57e4 → 7.77e4 (margin against 1.54e5: 2.34 → ~1.98). γ′₂ ≈ 0.143 → 0.124 [est], so E^A moves slightly (incl. the C1 rerun's 0.048).
- Either justify sup|z₂| ≤ 3 on the core, or rerun (A) with s33=ibp in a scratch copy (gonogo/lemmaS_A_rerun pattern).
- **P3 done** (launch/LAUNCH.md, launch_lemma.tex: NOT compiled, no LaTeX on this machine).
  - Launch rule: present centres at τ₁, from (H2)+(H3); late centres at ρ/λ = 56, from E^A plus the neighbour term.
  - ν = 2 is ruled out (neighbour factor 2.10 [est]). Box by radius (zone III from Lemma O D5a).
  - Order: verifier → D5a → (A) with cutoff [350, 700] → U-crude → U-fine → D5b/H6.
  - Launch-limited M1 for (4,4,0.3): ×2.13 / 2.05 / 1.91 [est] (KS = 1/2/4); conservative variant ×1.98 / 1.81 / 1.57.
  - T1 with (H3): ×1.666 (KS=1), ×1.406 (KS=4).
  - T3 fixed: Young route adds 3.2e-5 to 4.1e-5 to E₂ᴬ [est].
  - "To certify" list: LAUNCH tex §9.
- **(A) with core s33 = IBP**: Micah approved a rerun in scratch (`gonogo/lemmaS_A_ibp/`, run_chain_ibp.sh).
- **(A) with core s33 = IBP: DONE [cert]** (scratch copies; v3.1 untouched).
  - `gonogo/lemmaS_A_ibp/`, tightened (H2):
    - E = (1.2110055e-7, 8.1741747e-5, 0.048044859), against 0.048044818 with the assumed sup;
    - γ′₂ = 0.12355 (was 0.14306); core Ĝ₂₂ = 77720.9, margin 1.9771 (was 2.3378);
    - every other margin is identical: S0v 5.38, M0 6.58, M1 4.86, ramp, windows, |f| ≤ ½, K2/K3.
  - `gonogo/lemmaS_A_ibp_v31H2/`, v3.1 (H2):
    - E₂ = 0.24704515 (v3.1: 0.24704507); margins identical to v3.1 except Ĝ₂₂ (1.98);
    - quotes that become wrong: γ′₂ 0.143 → 0.1236; Ĝ₂₂ core 6.58e4 → 7.77e4; and in the edge-free (B) lines, E_B (y_E = 50) 0.250 → 0.25014, γ′_B 0.2394 → 0.2376, I3 const 1.84e-21 → 1.8420e-21.
  - **Conclusion:** v3.1 Lemma S(A) holds with the {2,3,3} class by IBP as the tex claims. Its certified script used the superseded assumed sup, and some quoted constants need an erratum, to go in v3.3 (new files) after Micah reviews.
- **P2 done** (ramp/P2_REPORT.md).
  - Certified K table. Effective KS ≤ 1.33 (from the new K[1][3]), against a flip at 14.8.
  - C_η^Ξ ≤ 13.349; footw = (1.94, 2068.7); S33 coefficient 99.2, including the (2,2,2) cutoff term.
  - **Open:** the level-2 cutoff cubic term (2,2,3), weight 7.33, is unbounded.
  - **Possible (A) v3.1 omission:** the cutoff × jet terms on [14, 28] (flagged, not established).
  - **Convention clash:** P2's conventions give Ĝ₂₂ ≈ 2.77e5 (Gw ≈ 1.82, rate ×2.22), against P1's 2.40e5. The checkpoint uses the worst.
- **P4 (Sonnet)**: checkpoint/t1_checkpoint.py, running.
- **P4 checkpoint run done** (checkpoint/CHECKPOINT_NUMBERS.md). Worst-Ĝ row: late centres M1 ×2.02/1.93 (launch), ×1.86/1.68 (conservative); the "present conservative" row is ×1.34, which Fable rules an impossible configuration.
- **Fable checkpoint: CONTINUE WITH CONDITIONS** (CHECKPOINT.md). Fable's scratch scripts copied to checkpoint/fable_scratch/.
  - **REAL GAP in closed v3.1 Lemma S(A), confirmed by the coordinator in code:**
    - the ramp pair (j+1, j+1) on [14, 28] is zeroed (lemmaS_F2.py:81 `out[j+1] = 0`, after `_expand` returns the square); lemmaS.tex:209 claims "only J_k, k ≤ j";
    - via J₃ it would give E₂ +3.0; the proposed fix (IBP + Young against the ramp dissipation, 0.32 a_min spare) gives +0.005 to 0.010 [est]; v3.1 core M1 ×2.82 → ~×2.75 [est];
    - also: cutoff × jet terms unbounded on (A)'s ramp (≤ 1.5e-3 on E₂ [est]); U_j (j ≥ 2) dropped in weights_exact.py:81–83 (core Ĝ₂₂ 77721 → 79271);
    - ramp data X_{j+1} assume 2a_min retained dissipation that (A) does not export.
  - Recommends Phase 2a (10–14 h): repair (A) in new files, then a second kill check. Total remaining 55–70 h.
  - **STOPPED here for Micah's decision.**
