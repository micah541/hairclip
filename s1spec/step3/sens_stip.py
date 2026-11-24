import numpy as np, gram_float as GF
from scipy.integrate import solve_ivp
import pickle
T = pickle.load(open('tail_mid.pkl','rb'))['Tm']
base = GF.STIP
res = {}
for d in [-5.4e-5, 0.0, 5.4e-5]:
    GF.STIP = base + d
    sol, up = GF.run(T=T, sig_end=base + d)      # end at the shifted tip
    G4, G3, tip = GF.gram(sol, up)
    za, zb = tip['za'], tip['zb']
    W13 = za[1]*zb[3]-zb[1]*za[3]
    Gs = (G3+G3.T)/2
    res[d] = (W13, Gs, tip['ws'])
    print(f"d={d:+.1e} W13={W13:.6f}  G3 diag={np.diag(Gs)}  G3[0,2]={Gs[0,2]:+.6e} ws0={tip['ws'][0]:+.6e}  mineig={np.linalg.eigvalsh(Gs)[-1]:+.4e}")
GF.STIP = base
print("rel dW13:", (res[5.4e-5][0]-res[-5.4e-5][0])/2/res[0.0][0])
print("dG3 (half-diff):\n", (res[5.4e-5][1]-res[-5.4e-5][1])/2)
