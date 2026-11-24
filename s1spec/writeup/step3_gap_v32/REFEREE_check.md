# Referee check: the far-field localisation gap (`gap_v32.tex`)

**Verdict: the gap is real.** The certified Lemma S(B) window inequalities leave out every fourth-order term at the window edges. No argument in the sources (text, scripts, `rigor/`, verifier) dismisses them. The main trapping argument depends on the crude far-field box that these inequalities produce. The note states the scope essentially correctly, but has two small inaccuracies (§4).

Paths below are relative to `/home/claude/s1spec/`. "v3.1" means `writeup/lemmaS/` and `writeup/step3_proof/`. "v3.2" means `writeup/lemmaS_v32/`, `writeup/step3_proof_v32/` and `writeup/lemmaO_v32/`.

---

## 1. Are the window-edge terms included anywhere? No.

**What the identities are.** They hold only up to total derivatives, and the dropped total derivatives are exactly the edge terms.
- `writeup/lemmaS/lemmaS.tex:101-103` (Prop. 1): each identity holds "$+\partial_s(\dots)$".
- `lemmaS_F1.py:19-28`: `div_form` keeps the quadratic form and checks only that the difference is a total derivative (`is_total_derivative`, l.28; this is the Euler-operator test in `lemmaS_prop1_sym.py:24-30`).
- `lemmaS_B3.py:7-15` does the same at level 3.
- On a window $[s_a(\tau),s_b(\tau)]$ these total derivatives become edge values of the form $\varepsilon\hat\omega\,(z_jz_{j+3},\ z_{j+1}z_{j+2},\dots)|_{s_a}^{s_b}$. With a cutoff they become ramp commutators instead.

**Text.**
- `lemmaS.tex:174` (Prop. B): "Same computation as F1 with no ramp inflow". Its only extra costs are the $b$-flux factor and $e^{Cy_E^{-4}}$, with $C=\max\hat G_{jj}$, which is an interior constant.
- `lemmaS.tex:239` ((B) level 3): "no inflow".
- `step3_proof.tex:233`: "the $a$-flux vanishes and the $b$-flux costs …". This covers transport only. The fourth-order edge flux is not mentioned anywhere.

**Scripts.** The (B) energy is entry × growth + forcing + coupling, with no edge or ramp term.
- `lemmaS_Bbox.py:35`: `EB = Es*exp(Cexp)*exp(8*abar/yE) + 4*(r0*p**j)**2*IB/gp**2 + 2*coup/gp`. Only the interior `YC`/$\hat G$ is used (l.17-20). The inner-ramp `ramp_k` that `F.constants` returns is discarded.
- `lemmaS_F4B.py:44-53`: the same formula at l.50, for the windows entering through $y=48,50$, i.e. the crude bounds for windows launched later.
- `lemmaS_B.py:20-21`: the same formula. Its docstring at l.1 says "no ramp inflow".
- `lemmaS_B3.py:30-35`: $\hat G_{33}$ and $\gamma_3'$ only.

**`rigor/`.** `assembly_arb.py:12-13` ports Bbox and the (B) exit of F4B "line by line". The arb formulas at l.149 and l.188 are identical, with no edge term. `RIGOR_PASS.md` does not mention edges.

**Verifier.** `step3nl/verifier/constants.py:61,65` takes only `S0V`, `S1V`, `FLUX0`, `FLUX2` from Lemma S. Its own edge is the $E_I$ cutoff on $x\in[57,60]$, at a fixed rescaled position inside the (A) core (CFLUX, `lemmaS.tex:140-148`; `step3_proof.tex:86`). Nothing in it concerns (B).

**Outer end of (A).** Also unquantified.
- `lemmaS_F1.py:2` says "outer cutoff far away".
- `assemble` (l.127-148) contains only: hand-off, inner-ramp inflow, inner-ramp commutators, forcing, coupling, and `extra_outer`.
- The forcing is integrated up to `Rout=350` (l.146), so the shell energy is implicitly supported up to $y\approx350$. No cutoff term exists there.
- `lemmaS_F2.py` has only two coefficient boxes, core and ramp (l.2, l.115-116). There is no box for $y>60$.
- `step3_proof.tex:272,290` speaks of "the outer ramp $y>60$" of (A) but gives no computation for it.
- **Correction:** `extra_outer` is *not* an outer-end term and it is *not* zero. It carries the inner-ramp foot and S33 data from `F2.main` (`lemmaS_F2.py:187-214`). In the final `lemmaS_F10p_M0.52.json` it is $3.26\cdot10^{-4}$ (foot) plus $4.75\cdot10^{-4}$ (S33) at level 2. The note's "`extra_outer`$\,=0$" misreads the variable name. Its conclusion, that the outer end is unquantified, still holds.

## 2. Can the omission be dismissed by an argument already in the sources? No.

**(a) Physically fixed windows, so no material flux.** Fails. This argument cancels only the transport boundary term $-c\hat\omega z_j^2$ against the motion of the edges (`step3_proof.tex:233`; transport identity at `lemmaS.tex:103`). The fourth-order operator is diffusive, so its flux through a moving edge is not zero.
- Minor: the edges move at $(\frac14-a)$, not at the full characteristic speed. So a $b$-part of the transport flux remains, and it is also a pointwise edge term. v3.1 folds it into an energy factor $e^{8\bar a/y_E}$ (v3.2 uses $e^{29q\bar a/y_E}$) without a trace bound. It is of size $\sim10^{-8}$ and can be removed by letting the edges follow the exact characteristics, so it is not material.

**(b) Sign or cancellation.** Fails. Fourth-order edge terms such as $z_jz_{j+3}-z_{j+1}z_{j+2}$ have no sign. The sources never analyse the dropped total derivatives; they are discarded by construction (§1).

**(c) Smallness of $\varepsilon$ and of the time ratio.** Gives smallness of the coefficient only, not a bound.
- In the frozen normalisation $\int\varepsilon\,d\tau=q\,y_E^{-4}$, and v3.2's certified ramp coefficients are tiny: `lemmaS_v32/ramp_in_50_75.json` gives $c_{22}=2.4\cdot10^{-4}$ and $c_{23}=3.35\cdot10^{-5}$.
- These coefficients multiply ramp norms $N_k^2$ ($k\le j+1$) for which nothing in step 3 gives a bound in $60\lambda(\tau)\lesssim r<150\lambda_1$. In particular there is no time-averaged $\int_{\rm ramp}z_3^2$ there.
- So the gap is a *missing estimate*, not a quantitative failure. The numerical room is large: the derived $|f|\le0.093$ is a factor 5.38 below $\frac12$.

**(d) Covering or overlap (summing windows, or using neighbours).** Not in the sources.
- Summing the windows into a partition removes interior edges but leaves the two ends of the region. Done with power weights, this is the note's Route P, which is incomplete (note §4).
- A natural variant that the note does not discuss: bound the ramp norms of one window by the energies (B3) and the level-1 dissipation of overlapping windows. The sources contain the template, `lemmaS.tex:183` (F4: $\int z_3^2\le D_1/(a_{\min}\varepsilon\omega)$), but only at a fixed $y=57.75$ inside (A).
- For physically fixed windows this variant is circular: a level-2 edge needs the $z_3$ that comes from level-1 dissipation, and a level-1 edge needs the $z_2$ that comes from level-2 energies.
- My float estimate (not certified): $q\bar N_3^2\lesssim\hat E_1y_e^4/(a_{\min}\hat\omega)\approx1$–$2.5\cdot10^2$ for hand-off-level $\hat E_1=10^{-5}$ at $y_e=50$. That is about 40–80× above the $N_3^2\lesssim3$ that the strip needs (`step3_proof_v32.tex:394`).
- For windows entering at $y=48$, using $E_1^A$ and the $[50,75]$ coefficient as a proxy, the result is about 0.05–0.1, against a crude $E_2\approx0.25$. That is not negligible. So the variant does not dismiss the gap, and making it work would be a new argument.

**(e) Estimates elsewhere in the sources.** None covers these edges.
- **CFLUX/F4** (`lemmaS.tex:140-148,183`): edge at a fixed rescaled $x\in[57,60]$. It works because (A) supplies $z_{0..2}$ there and its level-1 dissipation supplies $z_3$. A physically fixed edge leaves the (A) core. In the frozen weight, (A)'s dissipation then loses a factor $\omega_e/\omega_A=e^{(\tau-\tau_e)/4}$.
- **Lemma O / B4 at $[150,175]$:** the outer edge in v3.1. `lemmaS.tex:154` asserts "only through $R_1^{-4}$" without any computation, and Lemma O bounds only $h_{xx}$ (sup) and $h_{xxx}$ ($L^2$), not the $z_3\sim h_{xxxx}$ that a level-2 ramp needs. v3.2 adds Lemma R for this, for the outer edge only (`lemmaR.tex:94`). I did not re-verify Lemma R.
- **Old uniformly-local Part (B)** (`lemmaS.tex:135-138`, cosh weights, "no cutoffs"):
  - It is marked superseded (l.91).
  - It is a sketch: it reuses power-weight $\hat G$ for cosh weights.
  - It does not treat the ends of the graph region.
  - The note's "proposed next step" is precisely a rigorous version of it.
- **`step3nl/BOOTSTRAP.md:81-93` ("O1 annulus" short-time stability):** a plan only. Lemma O's certified $t^*$ is far too small for cylinders that fit at $r\approx50$–$75\lambda_1$ (ratio $\gtrsim2.6\cdot10^{-6}$, `lemmaR.tex:94`), and the same holds at $r\approx48\lambda(\tau_e)$.

**About the note's "$k\le j+1$".** This requirement depends on the method.
- It holds for the F1-style Cauchy–Schwarz tables that v3.2 certifies (`ramp_in.py`/`ramp_out.py`, where $c_{j,j+1}\ne0$).
- The linear commutators can be reduced to $k\le j$ by applying Young against the window's own dissipation (`lemmaS_F2.py:2-4`). `step3_proof.tex:209` even asserts "order $\le j$ only" for the (A) ramp.
- The quasilinear $\{2,3,3\}$ foot term on a ramp still needs a cutoff-free $z_3$, or a ramp sup of $z_2$ (`lemmaS.tex:321,330-332`).
- So level 2 needs $z_3$-type data in any method in the sources. Level 1 then needs $z_2$, which in turn rests on level 2. The note's scope "levels 1 and 2" therefore stands.

## 3. Does the main trapping argument depend on the crude far-field box? Yes.

- **Verifier inputs.** W1 (trapping) comes from the verifier (D1, D6). Its inputs $S_0^v$, $S_1^v$, $K_0$, $K_2$ come from (A) (`step3_proof.tex:270,282`; `constants.py:61,65`).
- **(A) uses B3 beyond $y=60$.** D3 takes B3 in the form "$|V|\le\frac12$ on the outer ramp $y>60$" (`step3_proof.tex:272,290`). This is the coefficient box on (A)'s support beyond the core, which extends to $y\approx350$ (`F1:146`), together with (A)'s unquantified outer end.
- **B3 must be re-derived strictly.** The continuity argument needs it (`step3_proof.tex:262,315`). The re-derivation is D4 (`l.273,305`), and D4 is exactly the edge-free computation: `F4B:44-53` for the windows launched later, `Bbox:35` for the hand-off windows (also `step3_proof.tex:233`).
- **Continuation.** Prop. global uses "(B) and Lemma O beyond" to control curvature for (P1) (`step3_proof.tex:321`).
- **Conclusions.** W3 and W4 rest on (B) (`step3_proof.tex:104-105`).
- **H6 in v3.1.** The dependency table (`l.273-274`) has D4 take B4 on $[150,175]$, and D5 take D4 on the strip. The note's claim that the trapping argument is "independent of H6" therefore relies on v3.2's restructuring (D5a/D5b, Lemma R), which I did not verify.

## 4. Verdict and scope

**The gap is real.** Not proved:
- the window energies of Lemma S(B), at levels 1 and 2, for every window edge in $60\lambda(\tau)\lesssim r<150\lambda_1$;
- hence the crude far-field bounds ($E_B$, $|f|\le\frac12$, curvature $\le S_1^f/r$; B3, W3–W4);
- hence the coefficient box and the outer end of Lemma S(A) as D3 uses them;
- hence the closure of the bootstrap and the continuation step;
- also the hand-off-level strip bounds, and through them H6, D4 and D5.

**Not affected (agreed with the note):**
- the verifier, given its inputs;
- the core and inner-ramp computations of (A), given a crude box and an outer-cutoff bound;
- Lemma O and Lemma R as local lemmas (the latter not re-checked).

**Corrections to the note:**
1. "`extra_outer`$\,=0$" is wrong. That variable is the nonzero inner-ramp foot/S33 term; the right evidence is the list of missing outer-end terms in §1.
2. "$k\le j+1$" is specific to the method. With Young against dissipation, the linear part needs only $k\le j$, but the obstruction at level 2 remains through the quasilinear foot term.
3. Labels: in both versions, $|f|\le\frac12$ and curvature $\le S_1^f/r$ are in W4; W3 holds $E_B$ and $|f|\le0.093$.

None of these changes the conclusion. Severity: this is a missing estimate with large numerical room, not a failed inequality. A repair needs a time-averaged ramp bound on $z_3$ (and on $z_2$ for all times) in $r<150\lambda_1$, or uniformly local weights. Neither is in the sources.

**Side remark (outside the claim).** v3.1 evaluated the (B) forcing in entry coordinates (`Bbox:31`, `F4B:45`), and $R\propto y^{1.02}$ is not suppressed. v3.2 has already found and fixed this with the frozen profile (`step3_proof_v32.tex:235-262`, Prop. B′). I did not re-verify the fix.
