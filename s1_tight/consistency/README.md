# Consistency of the tight S1 enclosures with the step-3 hulls (outA5.txt / ibp_boxes)

**Conclusion.** Wherever they overlap, every tight S1 enclosure lies inside the interval hulls of the wide CAPD run `s1spec/step3/outA5.txt`, read in the conventions of `step3nl/agent_mod/ibp/ibp_boxes.py` (the boxes used by the IBP modulation lemma and step 3). So every downstream object refers to the same certified set of shrinkers. What `check_outA5.py` tests (ball arithmetic throughout, all numbers read back exactly), and what it found:
- **Scope.** The test covers all 1,025 stage-C rows (`MC`, σ ∈ [σ_tip, 8]) and the 501 stage-B rows (`MB`) that meet σ ∈ [8, 13].
- **The tight side.** For each outA5 row, the tight enclosure is restricted to the row's τ-range and taken as the union, over the tight rows meeting it, of the exact range of the Lagrange form f(t_c) + f′(t_c)h + f″(hull)h²/2.
- **Coordinate mapping.** On stage C, σ = 8 − c·t with c = τ_tip/T0, using each file's own τ_tip interval. On stage B, σ = 32 − t.
- **Containment tests.** This range is tested for containment, component by component:
  - prof_tip.txt and prof_tight.txt in MC: p, q, k, m, cos θ in Ct2 and −sin θ in St2;
  - ext813_prof.txt in MB: k, m, p and q − σ in the `gen_hs2.chunks` envelopes, q in σ + dq, and (cos K, sin K) in (C, S), where K = ∫_σ^32 k = θ(13) − θ(σ) + ∫₁₃^32 k and the last term is bounded by the hstail envelope (≤ 9.9e-8);
  - σ_tip in the `st` interval.
- **Result.** All rows pass for all three files (`check_outA5.json`, `check_outA5.log`).
  - Smallest relative margins: St2 1.3e-9 at the σ = 8 row, where the outA5 box is pinned by θ(8) = 0 to within 2e-15 and the tight θ is pinned the same way; q 4.6e-7 (MC); q 1.8e-9 (MB, where the q box is σ + dq with σ spanning the row).
  - All other margins are ≥ 1.9e-4.
- **Why containment is expected.** It is consistent with how the runs were seeded. outA5's stage C starts from the step-1 wide set (ε ∈ [1.5345849609e-4, 1.5346679687e-4] ⊃ the tight bracket) with st ∈ [−2.249127, −2.249019], which contains the tight σ_tip.
- **Rerun.** `python3 -I -B check_outA5.py` takes about 12 s.
