import sys, os, bb
a=float(sys.argv[1]); L=float(sys.argv[2])
r=bb.bound_cell(0.0,a,L*(1-1e-9),L*(1+1e-9),tol=5e-4,rtol=0.005,floors=False)
print(f"QM={os.environ.get('QM',0)} QB={os.environ.get('QB',0)} S1V={os.environ.get('S1V',0)} GN={os.environ.get('GN',0)} alpha<={a:g} Ly={L:g}: U={r['U']:+.5f} lb={r['lb']:+.5f}", flush=True)
