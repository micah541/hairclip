# Step 2: linear stability of the CDF shrinker S1, computer-assisted proof

This bundle extends the step-1 bundle (`s1_cap/`, existence of S1). The statement and proof are in `writeup/step2_proof.pdf`.

**Theorem.** Work in L²(e^{−¼(1+x²)^{2/3}} dx), where x is arclength from the tip. Every eigenvalue of the linearized rescaled CDF operator L at S1 with Re λ > −1/8 is one of the four symmetry modes:
- even class (reflection-symmetric): λ = 1 (time shift) and λ = ¼ (translation along the axis);
- odd class: λ = ¼ (translation across the axis) and λ = 0 (rotation).

Each is a simple zero of the Evans function of its class.

## Layout

| path | content |
|---|---|
| `lincheck.py` | symbolic check of the linearized operator L and of its four symmetry eigenfunctions |
| `evans.py`, `evans_norm.py`, `winding.py`, `locate.py`, `realscan.py`, `numrange.py` | non-rigorous numerics: spectrum, large-region windings, numerical range |
| `apr/sdp.py` | finds the numerical-range certificates (Clarabel SDP); not part of the trust base |
| `apr/verify.py` | **rigorous** check of a certificate in Arb on the whole real line: CAPD profile boxes, tail-theorem bounds, and a u = x^{-1/3} chart at infinity |
| `apr/cert_b0.25_th{0.0,0.48}_tau0.03.json` | the two certificates (Λ₀ = 1.790819 at θ = 0, Λ₁ = 1.768282 at θ = 0.48) |
| `lin/lintail2.py` | **rigorous** Lyapunov–Perron theorem for the admissible plane at σ₀ = 32 (frozen block basis, Krawczyk root enclosures, log-norm bounds) |
| `lin/check_alg.py` | numerical check of the block-basis identities |
| `cap/evans_seg.tmpl.cpp`, `cap/gen_vf.py` | CAPD program (template plus generator of the vector fields) for the Evans function on a λ-segment |
| `cap/seg2.py` | builds the CAPD input for one segment (tail theorem → initial Lohner set) |
| `cap/contour.py` | adaptive driver along the upper half of the contour; writes `contour_result.json` |
| `cap/winding_cert.py` | certified winding numbers from the segment enclosures (argument cones; positive per-segment scales allowed) |
| `cap/interior.py` | Lyapunov–Perron bounds on a cover of the interior (analyticity of the Evans functions) |
| `run_step2.sh` | reruns every rigorous step |

## Requirements

CAPD 6.x (`CAPDCFG=/path/to/capd-config`), python-flint ≥ 0.9, numpy, and the step-1 bundle (`S1=/path/to/s1_cap`). cvxpy with Clarabel is needed only to regenerate the certificates.

The scripts contain absolute paths from the machine where they were run (`/home/claude/...`). Adjust them, or recreate that layout, before rerunning.

## Trust base

- CAPD's interval Taylor solver;
- Arb ball arithmetic;
- the step-1 theorem (profile enclosures and tail bounds);
- the hand arguments written out in the PDF: the quadratic-form identity and null forms, the Lyapunov–Perron estimates, and the argument principle with convex cones.

The order-of-zero = algebraic-multiplicity statement is quoted from the Evans-function literature. It is not reproved here.
