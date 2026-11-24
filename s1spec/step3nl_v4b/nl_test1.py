import numpy as np, time
from nlflow import Flow, run
fl = Flow(X=34.0, h=0.05)
x = fl.x
for amp in [1e-4]:
    u0 = amp * np.exp(-x**2 / 4)
    t0 = time.time()
    H, v = run(fl, u0, dt=2e-3, T=12.0, every=1.0)
    r = np.diff(np.log(H[:, 1])) / np.diff(H[:, 0])
    print("decay rates per unit tau:", np.round(r, 4), " time", round(time.time() - t0, 1))
