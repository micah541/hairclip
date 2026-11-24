import numpy as np, gram_float as GF
orig = GF.prof
def prof_cut(s, freeze_tail=True):
    P = orig(s)
    if s > 8.0: P = np.array([0.0, P[1], P[2], P[3]])
    return P
sol, up = GF.run(); G4a, G3a, tipa = GF.gram(sol, up)
GF.prof = prof_cut
sol, up = GF.run(); G4b, G3b, tipb = GF.gram(sol, up)
GF.prof = orig
np.set_printoptions(precision=3, linewidth=150)
print("G4 change from p-forcing on sigma>8:\n", G4a - G4b)
C, S, ws, za, zb, Ia, Ib, yy, II = up(sol.y[:, -1])
print("y_p(tip) =", yy[0])
