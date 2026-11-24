s=open('constants.py').read()
old='SRC = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))   # .../step3nl  (read only)'
assert old in s
s=s.replace(old,'SRC = "/home/claude/s1spec/step3nl_v33"   # G0_INDEP copy: inputs read from the proof tree (read only)')
s+="""
# ======================= G0_INDEP review patch (NOT part of the proof tree) =======================
# G0FIX = set of letters switching on the terms the reviewer found missing/under-covered in rhs_iv.py (see ../G0_INDEP.md sec. 1-2).
G0FIX = set(os.environ.get("G0FIX", ""))
if 'c' in G0FIX:   # typed constants that were rounded DOWN from their certified values -> outward decimals
    AF = 6.95551e-2; CF = (0.44215, 0.5259, 0.42353, 0.42502); DS = (6.841e-3, 1.54184e-2, 2.36482e-2, 4.5627e-2)
    SB, SA, IA = 3.59552e-5, 1.33721e-4, 9.59971e-12
    TQ = dict(TQ, d010=1.6712e-3, d110=1.0865, d300=3.3141); FNW = (1.093e-7, 6.5e-9, 0.1823)
G0 = dict(xa1=0.0600, xb1=0.42, xa0=0.0394, xb0=0.549,          # top-order cross terms, level 0 (V1V2, V0V2)
          cb0=0.2596, cb1=3.592, b0=BETA0, b1=BETA1, bq=QB, bF=BETAF,   # beta coupling: (cb0 + cb1|al|) |beta| V0 ; cb1 = 2*1.796
          cpsi=1.88e-7,                                           # (||zeta psi1_x|| + ||zeta psi2_x||) |(a,b)|
          d200r=0.01331, d110r=0.2968, d300r=3.024,               # ramp [10.25,16.25] nonlinear monomial sups (v2_T3, alpha<=1e-3, eta 1e-2)
          cp1=6.5e-4, cp2=0.42, e1=5.7e-6, f21=2.701e-8, f22=1.156e-7,
          r2lin=6.0, r3a=1.1, r3b=0.63, r3c=0.25)
"""
open('constants.py','w').write(s)
