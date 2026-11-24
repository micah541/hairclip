import numpy as np, json
from scipy.integrate import solve_ivp
import gram_float as GF
from gram_float import coef, MLT, SIG0, STIP
wa, wb = GF.lt_frame(); al, be = json.load(open('wstar_coef.json'))
def f(tau, Y):
    c = coef(SIG0 - tau, -GF.OMEGA); M = MLT(c)
    return np.concatenate([-(M @ Y[0:4]) - c['P1']*Y[0:4], -(M @ Y[4:8]) - c['P1']*Y[4:8]])
sol = solve_ivp(f, [0, SIG0 - STIP], np.concatenate([wa, wb]), method='DOP853', rtol=1e-12, atol=1e-16, dense_output=True)
for s in [32, 28, 24, 20, 16, 12, 8, 5, 2, 0, STIP]:
    Y = sol.sol(SIG0 - s); A, B = Y[0:4], Y[4:8]; W = al*A + be*B
    X = s - STIP; rh = np.exp(-(1+X*X)**(2/3)/4)
    print(f"sig={s:7.3f}  |wa|={np.abs(A).max():.3e} |wb|={np.abs(B).max():.3e}  w*={W[0]:+.3e},{W[1]:+.3e}  cancel={max(abs(al)*np.abs(A).max(), abs(be)*np.abs(B).max())/np.abs(W).max():.1f}  sqrt(rho)*w*0={np.sqrt(rh)*W[0]:+.3e}")
