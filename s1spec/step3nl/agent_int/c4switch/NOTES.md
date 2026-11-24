# c4switch: inner-region switch psi3^,h2^ -> psi3~,h2~ (C^4). Scripts: c4fun.py (fast corrections), int_s1.py, int_s2.py, int_s3.py. Results: c4switch_S1.json, c4switch_S2_1.json, c4switch_S3_*.json, c4switch.json.
S1: d = new - old = (deg-9 Hermite corr) - (deg-7 corr); zero for x>=28. Norms are box-sup based, rigorous.
S2: stages C, Bt (both functions via symbolic Taylor form), B, T (original analytic tail, extension unchanged). psi3: also reported with the old ray+chunk stage for sigma in [8,20] (1.0505e-6).
S3: n3_forcing with C^4 jets.
Caveats: new h2 residual is larger (2.77e-6 vs 1.57e-6): jet-4 corrections are ~1e-6 at interfaces and the residual contains u'''' ; worst at x=0 (h2~'''' corrected at the tip).
Q2 is evaluated on psi3~. Box widths: S1 400/piece, S2 width_scale 1, S3 original grids. Tail stage T taken verbatim from originals (valid since extension identical).
