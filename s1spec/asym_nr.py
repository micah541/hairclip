import numpy as np
def thmax(beta):
    A = 2*beta/3
    r = np.linspace(1e-4, 3, 300001)
    Re = -A**4 + 6*A**2*r**2 - r**4 - beta/6
    Im = r*(4*A*r**2 - 4*A**3 - 0.25)
    m = Im < 0
    t = np.min((-Re[m])/(-Im[m]))
    # also r with Im>0 never hurts for th in (0,pi/2) if Re<0 there
    assert np.all(Re[~m] < 0) if np.any(~m) else True
    return np.arctan(t)
for b in [0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45]:
    print(f"beta={b}: asymptotic th_max = {thmax(b):.4f}")
