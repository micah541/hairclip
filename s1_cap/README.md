# Computer-assisted proof of the CDF self-shrinker S1 (Conjecture 7.2 of the Ogden–Warren notes)

`s1_proof.pdf` (source `s1_proof.tex`) states and proves the theorem. This folder reproduces every rigorous step.

## Run
    # build CAPD 6.1.0 once:  git clone https://github.com/CAPDGroup/CAPD && cd CAPD && mkdir build && cd build && cmake .. && make -j
    pip install python-flint numpy
    CAPDCFG=/path/to/CAPD/build/bin/capd-config ./run_all.sh      # prints ALL STEPS PASSED  (log: run_all.log)

## Files
| file | role | rigorous? |
|---|---|---|
| `tail.py` | Proposition 2: Lyapunov–Perron constants on [8,∞), contraction check, γ⁽¹⁾ by validated quadrature (Arb, 256 bits) | yes |
| `gen_input.py` | outward-rounded parallelepiped x + C·r0 containing the conical-end states at σ = 8 for an ε-interval | yes |
| `s1proof.cpp` | CAPD: backward integration to σ = −1/2, then Poincaré map to q = 0; prints the enclosure of M(ε) = m | yes |
| `bracket_eps.py` | bisection on the certified sign of M (gives the ε* bracket) | yes (each sign certified) |
| `s1curve.cpp` | CAPD trajectory enclosures (7961 boxes) and tip enclosure on the bracket set | yes |
| `geometry.py` | double point / embedding / mirror / tail checks (lemmas L1–L3 in the PDF) | yes |
| `sanity_tail.py` | numerical cross-check of Proposition 2 against an accurate float integration | no (sanity only) |

## Key outputs
- M(1.45e-4) ∈ [0.0433869, 0.0434043] > 0,  M(1.62e-4) ∈ [−0.0393414, −0.0393249] < 0; Poincaré map validated on the whole interval.
- ε* ∈ [1.53458496e-4, 1.53466797e-4];  D ∈ [2.500798, 2.500813];  k(0) ∈ [1.539789, 1.539858].
- One double point, at distance [0.82859, 0.83170] from 0, crossing angle [52.29°, 52.62°]; rays at ±[106.645°, 106.682°]; turning [196.645°, 196.682°].

Trust base: CAPD's interval ODE solver and Poincaré map (FILIB intervals), Arb ball arithmetic and quadrature, and the hand estimates of Proposition 2.
