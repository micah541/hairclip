import sys, copy, numpy as np, json
sys.dont_write_bytecode=True
from cross import *
def nax(S): return len(axis_diag(S)['cross'])
S = load_state(f'{OUT}/state_t0.npz')
grid = list(np.arange(0.005, 0.6, 0.005)) + list(np.arange(0.6, 4.6, 0.02)) + [4.59, 4.598, 4.5985]
events = []; prev = (0.0, copy.deepcopy(S), nax(S)); rec = [(0.0, prev[2], S.kmax())]
for tg in grid:
    S.run_to(tend=float(tg)); n = nax(S); rec.append((S.t, n, S.kmax()))
    if n != prev[2]:
        # bisect between prev state and S
        ta, Sa, na = prev[0], copy.deepcopy(prev[1]), prev[2]; tb = S.t
        while tb - ta > 4e-4:
            tm = 0.5 * (ta + tb); Sm = copy.deepcopy(Sa); Sm.run_to(tend=tm); nm = nax(Sm)
            if nm == na: ta, Sa = Sm.t, Sm
            else: tb = Sm.t
        events.append(dict(t_lo=ta, t_hi=tb, n_before=na, n_after=n, kmax=S.kmax()))
        print('EVENT', events[-1], flush=True)
    prev = (S.t, copy.deepcopy(S), n)
json.dump(dict(events=events, rec=rec), open('logs/crossscan.json', 'w'))
print('final', rec[-1])
