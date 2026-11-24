"""G0_INDEP: certified B&B bound of the rate at one (alpha, Ly) cell for a given G0FIX set (env).  usage: cell.py a_hi L_lo L_hi floors(0/1)"""
import sys, time, json, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), 'vp'))
import bb, constants as K
a, Ll, Lh, fl = float(sys.argv[1]), float(sys.argv[2]), float(sys.argv[3]), sys.argv[4] == '1'
t = time.time(); r = bb.bound_cell(0.0, a, Ll, Lh, tol=5e-4, rtol=0.005, floors=fl)
print(json.dumps(dict(G0FIX=''.join(sorted(K.G0FIX)), RAMPMULT=K.RAMPMULT, a_hi=a, L_lo=Ll, L_hi=Lh, floors=fl, U=r['U'], lb=r['lb'], best=r['best'], secs=round(time.time() - t, 1))))
