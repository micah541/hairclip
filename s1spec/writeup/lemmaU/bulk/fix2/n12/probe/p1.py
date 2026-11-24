import sys
sys.path.insert(0,'/home/claude/s1spec/step4/p2/t9h')
from flint import arb, acb, ctx
ctx.prec=128
import s1chart
S=s1chart.S1()
print('Theta',S.Theta.str(15),'tip',s1chart.SIG_TIP.str(12))
I=acb(0,1)
for sig in (8,10,14,16,20,28,40,60):
    x=arb(sig)-s1chart.SIG_TIP
    c=S.far_chart(x) if sig>8 else S.tip_chart(x)
    z=c['G']*(-I*S.Theta).exp()
    print(sig,'y1-sig',(z.real-sig).str(6),'h1',z.imag.str(6),'th-Th',(c['th']-S.Theta).str(4))
