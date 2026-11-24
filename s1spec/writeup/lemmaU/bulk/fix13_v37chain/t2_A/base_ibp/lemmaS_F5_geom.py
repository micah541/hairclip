import numpy as np, json
al=5e-4; p=1.0246; yO=200.; hw=0.75*yO; L=50.*1.0  # cylinder y in [50,350], L=50 lambda1
res={}
for c in (1.4,1.3056):
    f=lambda y: c*al*y**p            # U(y)
    fp=lambda y: c*al*p*y**(p-1)
    y=np.linspace(50,350,3001)
    fc=fp(yO)
    fc=f(yO); slope=(f(y)-fc)/(1+f(y)*fc)
    # curvature in cylinder units (lengths scaled by L)
    curv=np.max(fp(y)/(1+f(y)**2)**1.5)*L
    # S1-ray chart slope at R1=175
    res[c]=dict(U200=f(200),U350=f(350),U50=f(50),Uprime200=fp(yO),max_rel_slope=float(np.max(abs(slope))),
                curvature=float(curv),slope_limit=0.2437,curv_limit=0.12,
                time_ratio=(3/150)**4,tstar=3.70e-7,margin=3.70e-7/(3/150)**4,slope_R1=float(f(175)))
for c,r in res.items(): print(c,{k:round(v,5) if v>1e-3 else v for k,v in r.items()})
json.dump({str(k):v for k,v in res.items()},open('lemmaS_F5_geom.json','w'),indent=1)
