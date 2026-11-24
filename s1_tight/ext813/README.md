# ext813: tight enclosure of the S1 profile on σ ∈ [8, 13]

This is the export for the step-4 hand-off (`s1spec/step4/STEP4_RESULT.md`, gap G1). It was built by route (a): validated outward integration from the certified σ = 8 state.

**Result.** We get a pointwise enclosure of (p, q, k, m, θ) on all of σ ∈ [8, 13]. Its width is at most **5.7e-10** (that is p; k, m and θ are at most 1.1e-10, q at most 5e-14). The target was ≲ 1e-8, and the old tail envelopes there were ~5e-5. All checks pass.

| file | content |
|---|---|
| `ext813_prof.txt` | the enclosure: 25,000 rows in the `prof/prof_tight.txt` format, then `END` |
| `ext813_points.txt` | enclosures at the exact points σ = 8 + j/4, j = 0..20 |
| `s1ext.cpp` | CAPD program: `prof/s1prof.cpp` with the time reversed |
| `ext813_check.py` | checks 0–5, plus a reader (`rows`, `taylor`, `d1`, `d2`) |
| `ext813_check.json`, `ext813_check.log` | check results |
| `run_ext813.sh` | rebuilds and reruns everything |

## Format
Each row has 23 numbers, exactly as in `prof_tight.txt`:
```
ta tb tc   p_lo p_hi q_lo q_hi k_lo k_hi m_lo m_hi th_lo th_hi    <- enclosure at tau = tc
           p_lo p_hi q_lo q_hi k_lo k_hi m_lo m_hi th_lo th_hi    <- hull over tau in [ta, tb]
```
- **τ = 8 − σ.** This is prof_tight's reversed time, so here τ ∈ [−5, 0].
- **Rows.** Row i (i = 0..24999) covers σ ∈ [8 + 5i/N, 8 + 5(i+1)/N], with N = 25000, so each row has length 2e-4.
  - Rows go outward, in increasing σ.
  - Within a row, ta < tb.
  - The rows tile [−5, 0] exactly: ta(i) = tb(i+1), tb(0) = 0 and ta(last) = −5.
- **Numbers.** ta, tb and tc are doubles: s = σ − 8 = 5i/N is computed in double and then negated exactly. The enclosures hold at exactly these values.
  - All numbers are printed with 17 significant digits, so they read back exactly.
- **`ext813_points.txt`.** Each line is `s  p_lo p_hi q_lo q_hi k_lo k_hi m_lo m_hi th_lo th_hi`, at s = σ − 8 = j/4.

## Coordinates and normalisation (as in prof_tight / s1prof)
- **State.** (p, q, k, m, θ), with ′ = d/dσ:
  - k′ = m, m′ = p/4, p′ = −kq, q′ = 1 + kp, θ′ = k.
  - Hence m = k′, p = 4k″ and q = σ + 4kk′ + 4J, with J = ∫_σ^∞ k′².
  - This is the full profile ODE, not a linearisation: J is carried by q.
- **θ(8) = 0.** That is, θ(σ) = ∫₈^σ k. The θ of prof_tight also starts at 0 at τ = 0.
- **ε.** Normalised by Z₀(8) = ε, with ε* in the `bracket2.json` bracket (Prop. 2′ of `s1_tight.tex`).
- **Step-3 coordinate.** x = σ − σ_tip = τ_tip − τ, with τ_tip ∈ [10.249077650575472, 10.249077651519443].
  - So this file covers x ∈ [τ_tip, τ_tip + 5] ≈ [10.249, 15.249].
  - Existing readers of prof_tight work unchanged. For example, `step3nl/psi3res3.py` stage C uses Xc = TT − tc and h ∈ [tc − tb, tc − ta].

## Pointwise enclosure and derivative bounds
For σ in a row, put h = σ − σ_c, with σ_c = 8 − tc and |h| ≤ 1e-4. Then
```
f(σ) ∈ f(c) + f′(c)·h + f″(hull)·h²/2
```
- c is the centre box, and f′ is the ODE evaluated on c.
- f″ is evaluated on the hull: p″ = −(mq + k(1 + kp)), q″ = mp − k²q, k″ = p/4, m″ = −kq/4, θ″ = m.
- Every higher derivative is a polynomial in the state (apply the ODE again). So the hull bounds the derivatives of any order needed by a quadrature rule.
- `ext813_check.taylor(row, h)` implements this form, and `d1` and `d2` give f′ and f″.

## Achieved width
Widths of the centre and point enclosures (`ext813_points.txt`):

| σ | p | q | k | m | θ |
|---|---|---|---|---|---|
| 8 | 1.6e-11 | 7.1e-15 | 1.8e-12 | 3.1e-12 | 0 |
| 9 | 3.2e-11 | 2.8e-14 | 5.2e-12 | 6.1e-12 | 3.3e-12 |
| 10 | 7.1e-11 | 3.6e-14 | 9.0e-12 | 1.1e-11 | 8.3e-12 |
| 11 | 1.3e-10 | 3.9e-14 | 1.8e-11 | 2.6e-11 | 2.3e-11 |
| 12 | 3.4e-10 | 4.3e-14 | 4.0e-11 | 5.8e-11 | 3.8e-11 |
| 13 | 5.3e-10 | 4.4e-14 | 7.3e-11 | 1.0e-10 | 4.3e-11 |

- **Maximum over all rows (centre):** p 5.62e-10, m 1.02e-10, k 7.32e-11, θ 7.32e-11, q 4.6e-14.
- **Pointwise:** the Taylor form above at any σ ∈ [8, 13] gives the same numbers to 3 digits. Its h-terms add at most 1e-13 (that is p; the others add ≤ 1.4e-14).
- **Hull widths** are the variation across a row: p 2.4e-7, k 4.0e-8, m 5.0e-8, θ 3.1e-8, and q 2.0e-4 (since q ≈ σ). They enter only through the f″·h²/2 term.

**Growth of the non-decaying modes.** From σ = 8 to 13 the centre widths grow by:
- ×34 for p, ×40 for k and ×33 for m;
- the WKB factor is e^{(E(13)−E(8))/2} = 31.2, which with the factors μ^{−1}, 1, μ gives 26.6, 31.2 and 36.7 for k, m and p;
- the ε direction, which is the decaying one, shrinks.

**Beyond σ = 13.** With the `Send` argument of `s1ext`, the maximum width is 1.4e-9 at σ = 14 and 3.3e-9 at σ = 15.

## How it was certified
1. **Seed.** The starting set is `../in_br2.txt`. This is the certified parallelepiped x + C·r0 of Prop. 2′, and it contains the conical-end state at σ = 8 for every ε in the tight bracket.
   - It is the same file, with the same CAPD 1.8.1 build, that produced `prof/prof_tight.txt`.
   - Step 0 of `run_ext813.sh` reproduces `prof_tight.txt` byte for byte.
   - On why this seed is valid, see the provenance note below.
2. **Integration.** CAPD integrates the profile ODE forward in σ from 8 to 13 with a C⁰ Lohner method (`C0Rect2Set`, order 25, automatic steps). The vector field is the s1prof field with the time reversed.
   - The enclosures come from CAPD's rigorous `SolutionCurve`, evaluated at each row centre and over each row.
   - The arithmetic is interval arithmetic with outward rounding throughout (filib, `-frounding-math`). No float step enters an enclosure.
3. **Assumptions.** Nothing beyond the σ = 8 state is assumed. On σ ≥ 8 the tail satisfies the same ODE.

## Consistency checks (`ext813_check.py`; results in `ext813_check.json` and `.log`)
0. **Seed provenance.** `../in_br2.txt` is rebuilt byte for byte by the certified generator, with a ∫K₁² ball that contains the jhat_rig enclosure (see the note below).
1. **Rows.** The rows tile τ ∈ [−5, 0] exactly; tc lies in [ta, tb]; and each centre box lies inside its hull.
2. **Widths.** As reported above.
3. **σ = 8 against prof_tight.** I took four enclosures of the state at σ = 8:
   - the hull of `in_br2`;
   - ext813 at s = 0;
   - prof_tight row 0, Taylor-transported to σ = 8;
   - ext813 row 0, Taylor-transported to σ = 8.

   They have a common point in every component. The common widths are 1.6e-11 (p), 1.8e-12 (k), 3.1e-12 (m), 6.5e-15 (q) and 0 (θ). The hulls of the two row-0 boxes meet.
4. **Tail envelopes.** Every row's range of k, m, p and q − σ lies inside `hstail.profile_bounds(σ_a)`. This also implies the `gen_hs2` chunk bounds.
   - The worst ratios |value|/bound are 0.909 (k), 0.929 (m), 0.930 (p) and 0.40 (q − σ).
   - So the envelopes are only 7–9% above S1, and this check has teeth.
5. **Independent of CAPD.** The tail2/laplace representation gives (k, k′, k″) = εŶ(σ) + P(σ)Z_nl(σ), where:
   - Ŷ = c_n(K₀, K₁, K₂)(σ) is exact (`laplace.py`), and ε lies in the bracket;
   - |Z_nl,0| ≤ rn0/w and |Z_nl,1| ≤ rn1/(wg), with rn0 = 8.6e-13 and rn1 = 4.4e-13 (`tail2.py`, in the step-1 weighted norms).

   Results at σ = 8 + j/4:
   - This representation meets the CAPD enclosure for k, m and p at every point.
   - For σ ≥ 8.5 it lies strictly inside the CAPD enclosure. Below 8.5 the bound rn0 + 2rn1 is the wider of the two.
   - At σ = 13 it is 2e-14 wide for p, against 5e-10 for CAPD, and the CAPD midpoints are within 4e-11 of it.

   If k, m or p are ever needed tighter, or beyond σ ≈ 15, this representation is itself an enclosure of width ≲ 3.5e-12·μ^i/w(σ). For q and θ one would also need ∫K₁² and ∫K₀.

## Note on the seed `in_br2.txt` (provenance)
- **Timeline.** `in_br2.txt` was written on Nov 6 at 18:14. That is earlier than the final `jhat8.json` (18:23) and `gen_input2.py` (18:25).
- **Reconstruction.** Check 0 rebuilds it byte for byte with the current generator and `bracket2.json`. The input for ∫₈^∞K₁² is taken as double(1.4404409694e-9) ± 1e-18, i.e. the 11-digit printed value with the generator's 1e-18 radius floor.
- **Validity.** That ball contains the jhat_rig enclosure [1.4404409694e-9 ± 6.45e-20]. So `in_br2.txt`, and with it `prof_tight.txt` and ext813, rest on a valid enclosure.
- **Consequence for `run_tight.sh`.** The current `jhat8.json` midpoint gives a seed that differs in 3 numbers, all q-entries (δq ≈ 1.7e-18). So `../run_tight.sh`, as it stands, would regenerate a slightly different `in_br2.txt`, and a `prof_tight.txt` that differs in the last bits. Both versions are valid.
- **Choice here.** ext813 deliberately uses the stored seed, so that it continues prof_tight exactly.

## Rerun
```
cd /home/claude/s1_tight/ext813 && bash run_ext813.sh      # about 17 s; one process at a time, nice 10
```
If `capd-config` is not on PATH, set `CAPDCFG=/path/to/capd-config`. The default is `/home/claude/ext/CAPD/build/bin/capd-config`.

The direct call is:
```
./s1ext ../in_br2.txt Send N [order [pointsfile NP]]
```
The shipped files use Send = 5, N = 25000, order = 25 and NP = 20.
