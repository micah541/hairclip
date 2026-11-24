import sys, os, importlib, numpy as np
sys.path.insert(0, '/home/claude/s1spec/step3nl')
os.environ['CUBFINE'] = '1'; os.environ['ALIN'] = 'spec'
for EPS in ('1e-5', '3e-5', '1e-4'):
    for MU, EPSI in (('0.08', '3e-5'), ('0.08', '1e-4'), ('0.03', '1e-4')):
        os.environ['MU'] = MU; os.environ['EPSI'] = EPSI
        sys.argv = ['x', EPS, '1e-3', '1e-6', '1e-6']
        import inner_closure_float2 as M; importlib.reload(M)
        cap = None
        for Ly in np.logspace(-13, -8, 21):
            if M.worst(0.0, Ly)[0] > 0: cap = Ly; break
        r0 = M.worst(0.0, 1e-14)[0] / 1e-14; ra = M.worst(3e-3, 1e-12)[0] / 1e-12
        print(f"eps={EPS} mu={MU} epsI={EPSI}: rate(a=0,1e-14) {r0:+.4f}  rate(a=3e-3,1e-12) {ra:+.4f}  cap Lyap {cap}", flush=True)
