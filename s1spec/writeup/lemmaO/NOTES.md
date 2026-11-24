# Lemma O (CDF graph form, local short-time stability) - notes

Files: `lemmaO.tex/.pdf` (statement, proof, appendix of inequalities), `lemmaO_consts.py` (C_j), `sanity.py` (symbolic + PDE check).
Run: `python3 lemmaO_consts.py K2 K3 M t1`;  `python3 sanity.py`.

**Result.** Graph equation verified symbolically (sign: leading term -h_xxxx). Proof = two cutoff energy estimates, E1 = int z^2 g_x^2 (Prop 3) and
E3 = int z^4 g_xxx^2 (Prop 4), plus a weak (mean) bound and 1D interpolation. Since (A) is assumed on the whole [0,t1], the energy
inequalities have *constant* right-hand sides (dissipation absorbs the top-order terms, K2/K3 bound the rest): E1<=L1 t, E3<=C* t, no Gronwall,
no continuity argument, **no smallness of t1** (delta0 = infinity formally; only C0 has a sqrt(t1) additive term). The forcing G(h0) is paired
with the dissipation exactly as suggested (needs h0^(5) in L2 only).

**Numbers (not useful at the requested times).**
| (K2,K3,M), t1 | C1 sqrt t1 | C2 sqrt t1 | sqrt(C3 t1) |
|---|---|---|---|
| (0.1,1,0.1), 3.2e-4 | 0.83 | 16 | 15 |
| (0.3,1,0.3), 3.2e-4 | 1.6 | 26 | 27 |
| (1,1,1), 3.2e-4 | 5.9 | 85 | 86 |
| (10,10,10), 1e-8 | 2.0 | 47 | 53 |

Deviations <= 0.1 need t <= 1.2e-8 / 4.4e-9 / 4.4e-10 (case (i), K2=M=0.1/0.3/1) and 3.5e-14 (case (ii)), i.e. t1 smaller by 3e4-7e5 than requested.
**Bottleneck:** `kappa = K3+2M` (UNWEIGHTED L2 norm of g_xxx on the cutoff transition zone) times cutoff curvature `q = 12Z^2+4Z2` in Prop 4
(kappa/10 -> C3 /75; floor ~1e4 from c3*K2^2, c4*K2^4 terms); for C1, C2: Lambda1, dominated by (cutoff gradient)*(unweighted |g_xx| <= K2+M)*s1 and by kappa again.
Unproved remedy: nested cutoffs feeding the small energy of the outer level into kappa (only helps for t1 <~ 1e-7, and E1 must be chained too).

**Sanity check.** Periodic [-4,4), Fourier + BDF, 4 data sets (smooth, k<=8, large curvature, H^5-marginal spectrum k^-6). All observed deviations are
below the bounds by 4-10 orders (observed g_xx ~ t, bound ~ t^{1/2}); energies E1, E3 with the actual cutoff also far below L1 t, C* t.

**Not rigorous / caveats.** Floating-point evaluation (no interval arithmetic); the max of the degree-8 polynomial in Prop 4 is by root-finding
(+1e-9 relative inflation); cutoffs are C^2 quintic smoothsteps, not C^infty (not needed). The hypothesis (A) is assumed, not derived.

## v2 (sharpening pass): `lemmaO_v2.tex/pdf`, `lemmaO_consts_v2.py`, `sanity_v2.py`
Cylinder [-X,X], additive slack t*(K2,K3,M,X) (slack 0.5M; L2 slack sqrt2*0.5M). Verdict: target t=4e-5 / 1e-5 NOT reached; X=2.4, K=1.5M, M=0.03-0.3:
t2* ~ 0.7-1e-7 (binding, pointwise h_xx), t3* ~ 2-5e-7, t1* ~ 3-5e-6 -> misses by ~410-570 (4e-5) / ~100-140 (1e-5). X=5 misses by 21-61 / 5-15.
Pairing the forcing with the energy does not help (deleting all forcing terms moves t* by <=28%); nested cutoffs make it worse (n=1 best).
Dominant: unweighted ||g_xxx|| (kappa=K3+sqrt(2X)M) x cutoff curvature; then L2->pointwise conversion for g_xx.
