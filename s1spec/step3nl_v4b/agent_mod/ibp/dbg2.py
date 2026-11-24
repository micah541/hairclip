import sys, time
sys.path.insert(0,'/home/claude/s1spec/step3nl/agent_int/v3'); sys.path.insert(0,'/home/claude/s1spec/step3nl'); sys.path.insert(0,'/home/claude/s1spec/step3')
t=time.time()
import c4fun as C
print('import',time.time()-t)
from flint import arb
for xv in (0.3, 5.1, 12.2, 20.5, 27.9, 31.0, 34.0):
    t=time.time(); X=arb(xv,0.005); j=C.PSI.jets_any(X,xv,5); h=C.H2.jets_any(X,xv,5)
    print(xv, time.time()-t, [float(abs(z).upper()) for z in j], [float(abs(z).upper()) for z in h][:3])
