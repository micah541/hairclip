import sys, numpy as np
from run_f2 import *
L = float(sys.argv[1])
S = load_state(f'{OUT}/base_{int(L)}.npz'); S.run_to(kstop=1.003 * L); save(S, f'{OUT}/baseX_{int(L)}.npz'); print(L, S.kmax())
