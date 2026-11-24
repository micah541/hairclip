import json, math
p=1.0246; al=1.5e-4; U1=lambda y:1.3056*al*y**p
C1, C2, G3, R4 = 0.00747719, 0.0395462, 0.0239123, 0.366035
s0, m2, n3, M = 0.08, 0.03, 0.004, 0.12
def T6(y0, split, tag):
    co=json.load(open(f'ramp_out_{y0:g}_379.5.json'))['c']; N=[0.0]*4
    for (a_,b_) in ((float(y0),split),(split,379.5)):
        L=50.; k2=m2+C2; k3=n3+G3; k4=math.sqrt(2)*M+math.sqrt(R4); pm=math.tan(math.atan(s0)+C1); Fb=math.tan(math.atan(0.0928)+C1)
        sp_,sf=1+pm*pm,1+Fb*Fb; ln=math.sqrt((b_-a_)/L); ell=math.log(b_/a_); Ub=U1(b_); Ks=k3+3*pm*k2**2*ln; D1=sf**1.5*(b_/L)*k2
        N2f=(b_/L)**1.5*sf**2*(sp_**0.25*Ks+3*Fb*k2**2*ln); N3f=(b_/L)**2.5*sf**2.5*sp_**0.25*(k4+10*(pm+Fb)*k2*Ks+(6+15*(pm*pm+Fb*Fb))*k2**3*ln)
        n=[(Fb+Ub)*math.sqrt(ell),(D1+p*Ub)*math.sqrt(ell),N2f+(D1+p*p*Ub)*math.sqrt(ell),N3f+3*N2f+(D1+p**3*Ub)*math.sqrt(ell)]
        N=[x+y*y for x,y in zip(N,n)]
    conv=379.5/28; D=[conv*sum(co[str(j)][str(k)]*N[k] for k in range(j+2)) for j in range(3)]
    print(tag, y0, split, 'N2=',['%.3g'%x for x in N],'Delta=',['%.4g'%x for x in D])
T6(250,330.,'orig'); T6(250,340.,'split340'); T6(252,340.,''); T6(253,340.,''); T6(253,335.,''); T6(253,345.,'')
