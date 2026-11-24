# PROOFS_REPORT: regularity lemma U-1 and the level-2 export G1 (fix round 10)

Nov 20, 2026. Directory `D` = `/home/claude/s1spec/writeup/lemmaU/bulk/fix10_proofs/`. Nothing outside `D` was written. The proof tree was only read.
Texts are not compiled (no LaTeX in WSL); brace and `\begin`/`\end` balance and `\ref` targets were checked by script.

Labels: **[hand]** written argument, not machine-checked; **[sym]** exact sympy, script in `D`; **[cert]** arb, script in `D`; **[cited]** literature or Lean, statement read from source; **[read]** taken from a project file by reading it; **[inh]** inherited from another step-3/step-4 text with that text's own label.

## Verdict
- **G1: proved, conditional only on items Lemma S(A) already rests on at k = 0, 1, plus U-1.**
  - Code and proof agree. No level-2 source escapes the window form, so the ROUTE_DECISION kill check did not fire.
  - C0_2 ≤ 116.3782 and C1_2 ≤ 0.2606004 [cert]. These are the values flux lemma v2 used.
- **U-1: proved [hand], given three inputs that are not proved here.**
  1. (H6), i.e. γ(t_1) is smooth. This is step 4's [A] parabolic-regularity statement.
  2. The Lean theorems. They are read from source but not compiled here, and the n = 1 instance has not been elaborated.
  3. DKS's smooth-extension step.
  - In return, U-1 gives more than was asked: Λ and every energy are C¹ in τ, and the identities hold at every τ, not only a.e.
  - It also proves (P2′) in the form P1_CONFIRMATION §4 asked for.
- **Step 3 is not called closed here.** Still open:
  - the items of FLUX2_IMPL_REPORT §4.3–4.6, other than C8 and G1;
  - an independent review of these two proofs;
  - the inherited inputs listed below.

## Files
| file | content |
|---|---|
| `U1_regularity.tex` | Theorem U-1, Proposition (P2′), Lemmas 1–8, gaps section |
| `G1_export_k2.tex` | Theorem (Lemma D, k = 2), Proposition (level 2 with r_2), Lemma (D2) proved source by source, status |
| `U1_sym.py` → `U1_sym.json`, `U1_sym.log` | sympy S1–S8, all pass, peak RSS 76 MB |
| `G1_check.py` → `G1_check.json` | arb, 36 checks, all pass |
| `SHA256SUMS` | hashes of the above |

Both scripts write an `args` block (argv, env, versions, sha256 of inputs and of the script). They were run with `PYTHONHASHSEED=0 nice -n 19 python3 -B`, one process at a time.

## U-1: step by step
| step | label | where |
|---|---|---|
| Lean velocity at n = 1 equals −k_ss N (general parametrisation) | [sym] S1 | Lemma 1 |
| reflection equivariance, reparametrisation invariance (either orientation) | [sym] S2, S3 | Lemma 1, Lemma 3 |
| curves are Lagrangian; J isometry, J² = −1; dim V = 2 dim E | [sym] S4, [hand] | §1.2 |
| LWP: `exists_flow_of_immersion` (L1) | [cited, Lean] | commit d72b744, Flow/Main.lean |
| uniqueness: `flow_unique_of_immersion` (L2) | [cited, Lean] | Flow/UniqueMain.lean |
| no dim ≥ 2 hypothesis in the final statements; no `sorry`/`axiom` in CWFlow/ | [read] | text search of the source |
| gluing two flows (joint smoothness at the junction) | [hand] | Lemma 2 |
| maximal smooth solution | [hand] + (L1), (L2) | Lemma 3 (`maximal`) |
| symmetry preserved | [hand] + S2, S3 + (L2) | Lemma `sym` |
| blow-up criterion: T_max < ∞ ⇒ sup ∫k² = ∞ | [cited, DKS step inside the Thm 3.1 proof] + [hand] | Lemma `bu` |
| frame identity ⟨Γ_τ,ν⟩ = −κ_ss + (¼ − a)⟨Γ,ν⟩ + b⟨e₁,ν⟩ | [sym] S5 | §3 |
| graph functions over S1 and over the ray are C^∞ in (t, λ, ξ, x) | [hand] (IFT) | Lemma `graph` |
| (P2′): modulation ODE, constraints conserved, continuation until t → T_max or the box fails | [hand]; matrix invertibility [read] MOD_IBP ν ≤ 0.0466 | Proposition `mod` |
| cutoff profiles: S7, S jets at the ends | [sym] S6 | §1.4 |
| differentiation under the integral; IBP with Lipschitz data | [hand] (+ S7 algebra) | Lemmas `leib`, `ibp` |
| D2v formula; level-0 IBP identity | [sym] S7 | Theorem U-1(b) |
| Agmon on an interval of length 3 | [hand] (+ S8) | Theorem U-1(c) |
| Theorem U-1 (a)–(d) | [hand] | §4 |

## G1: step by step
| step | label | where |
|---|---|---|
| r_export[2] = 1/1000; core retained 0.001001 ≥ r_2; recomputed 2 − B_c − used ≥ r_2 | [cert] B1 | `weights_t2.json` |
| plateau (60,175]: retained ≥ r_2 a_c (absolute), recomputed | [cert] B2 | |
| zones [40,50], [50,60], [175,250], [250,350]: no failing cell; level 2 and all_ok | [cert] B3 (flags read in arb from the certified run) | |
| ramp Young users at level 2: Σδ = 0.3255 ≤ slack 0.32585 | [cert] B4 | |
| no other user of level-2 dissipation on [28,350] (ramp users on [14,28], outer users on O) | [read] + [hand] | Prop. `level 2` |
| Steps 1–7 of Lemma S(A) at j = 2 | [inh] | lemmaS_A_v33.tex §2 |
| class list: handlers ⊂ {pair, phi, ibp223}, routes ⊂ {data, young}, top route IBP | [cert] S1 | |
| D_2 = sum of the parts; D_2^c ≤ D_2 | [cert] S2 | |
| (D2) window form, eight sources | [hand] + [read] (assembly_rep.py line by line) | Lemma (D2) |
| H-J-W on every window; Lemma S (cubic) on any interval | [inh] / [read] | lemmaU_v34.tex |
| integration (Lemma D proof at k = 2) | [hand] | §C |
| C0_2, C1_2 recomputed from the components; coup_2 recomputed | [cert] C1 | |
| setting is 4 N2_ref; N̄² ≤ 4 N2_ref; Lemma U's Ê_k ≤ the assembly's E_k | [cert] C2 | |
| flux2_data used exactly these C0_2, C1_2; the quoted 1.1638e2, 2.6061e-1 dominate | [cert] C3 | |
| A_interface upper bounds dominate | [cert] C4 | |
| use on the shell: ∫∫_{Y(S)} z_4² dy/y ≤ (y_max⁵/28) X_2 | [hand] | Remark |

## Gaps and findings (plainly)
1. **(H6) smoothness at t_1** is not proved here. It is step 4's [A] item: "smooth for t > 0 by parabolic regularity" (step4_proof.tex l.672). The Lean theorems need a smooth initial immersion, so they cannot replace it.
2. **Lean.**
   - The theorems are read from source, but the build was not run here.
   - The n = 1 instance has not been elaborated: E = EuclideanSpace ℝ (Fin 1), M = the unit circle, V = ℂ, J = multiplication by i.
   - The gluing lemma and the closure of `IsFlow` under σ∘·∘ι are written mathematics, not Lean.
   - Suggested check for Micah: a short Lean file instantiating both theorems at n = 1.
3. **DKS (D1)**, "bounded ∫κ² ⇒ smooth extension to [0, T]", is a step inside the proof of Thm 3.1, not its statement.
4. **(GS) at τ_1** is implicit in step 4's certification of Λ(τ_1) but not listed among the hand-off items. (GS) means the tip arc is a normal graph over S1 on |x| ≤ 60.
   - The hand-off list should state it, like the (H4) range item.
   - continuation_v33 should add (GS) to the bootstrap next to (G). The open/closed proof is given in U-1.
5. **Read, not re-derived:**
   - the graph equation (`graph_eq.py`, `graph_struct.pkl`);
   - the identification of the modulation matrix with G + E, and the Neumann bound ν ≤ 0.0466 (MOD_IBP.md);
   - the C⁴ matching of ψ̂3, ĥ2 (`c4fun.py`).
6. **Text items for G4. These are not errors in the export.**
   - lemmaS_A_v33.tex still says "r_2 = 0", "Level 2 exports nothing", B_c^{(2)} = 169/100 and "retained (1.25, 1.25, 0)". The v36 values are B_c^{(2)} = 1688999/10⁶ and retained 0.001001.
   - lemmaU_v34 Lemma D is stated for k ∈ {0,1} with r_k = 1.249. It holds verbatim for k = 2 with r_2 = 1/1000.
7. **V36_REPORT §6 rounds two design-setting quotes down.** It quotes C0_2 = 114.782 and C1_2 = 0.0698142; the certified values are 114.78222… and 0.069814238…. Nothing consumes the design setting, and A_interface's `C0_ref`/`C1_ref` are correct upper bounds.
8. **Inherited in G1 (s5):** F_max comes from the verifier's peak √Λ. That Λ ≤ L + 2C′ stays inside it is flux-lemma check C7. It was not redone here.
9. **Process note:** the first version of `U1_sym.py` simplified S1 with generic sympy functions, reached 2 GB RSS and was killed. It was rewritten in jet variables; peak RSS is now 76 MB.
