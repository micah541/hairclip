import sys, time, json, bb
# usage: cell.py a_hi Llo Lhi floors
a=float(sys.argv[1]); Llo=float(sys.argv[2]); Lhi=float(sys.argv[3]); fl = sys.argv[4]=="1"
t=time.time(); r=bb.bound_cell(0.0,a,Llo,Lhi,tol=5e-4,rtol=0.005,floors=fl,cem_iters=25)
print(json.dumps(dict(U=r['U'],lb=r['lb'],q=r['qstar'],best=r['best'],secs=round(time.time()-t,1))))
