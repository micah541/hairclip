# Tightened enclosure of the CDF shrinker S1 (step-1 sharpening)

Same theorem and ε-normalisation (Z₀(8) = ε) as the step-1 proof; only the conical-end state at σ = 8 is sharper.

| file | role | rigorous? |
|---|---|---|
| `laplace.py` | exact decaying solution K_j(σ) = (1/2πi)∫ t^j e^{σt+t⁴} dt of k‴ + (σ/4)k = 0 (arb quadrature + tail bound) | yes |
| `tail2.py` | Ŷ(8) normalised by Ẑ₀(8) = 1; Lyapunov–Perron bounds for the linear correction and the cubic remainder (step-1 constants) | yes |
| `jhat.py`, `jhat_rig.py` | Ĵ = ∫₈^∞ Ŷ′² via a double contour integral + explicit remainder (`jhat8.json`) | yes |
| `gen_input2.py` | outward-rounded CAPD parallelepiped for an ε-interval (ε-column carries dq/dε) | yes |
| `bracket2.py` | bisection on the certified sign of M(ε) (`bracket2.json`) | yes |
| `prof/s1prof.cpp` | CAPD: centre enclosure + hull on 51,246 boxes for the downstream certificates | yes |
| `run_tight.sh` | rebuilds and reruns everything | — |

Key outputs:
- State uncertainty at σ = 8: 4.4e−13 (step 1: 3.4e−8).
- M(1.45e−4) ∈ [0.04339455686, 0.04339455704], M(1.62e−4) ∈ [−0.03933409755, −0.03933409733].
- ε* ∈ [1.5346198892965e−4, 1.5346198899299e−4] (relative width 4e−10).
- τ_tip ∈ [10.249077650575, 10.249077651520], so σ_tip ∈ [−2.2490776515195, −2.2490776505754].
- D ∈ [2.50080627654, 2.50080627668], k(0) ∈ [1.53982063762, 1.53982063822].
- Geometry checks (double point, embedding, mirror, tail): all pass.

Downstream (in `s1spec/step3nl`):
- `psi3res3.py`: ‖(L−λ̃)ψ̂3‖_ρ ≤ 1.28e−6 on the half line (was 3.8e−3).
- `lyap_psi3b.json`: updated triangular-lemma constants.
