"""Check psi2 = Ct C + St S = -sin(theta), and the sign of cos(theta), against the float theta profile (as mod_form2_float.py)."""
import sys, numpy as np
sys.path.insert(0, '/home/claude/s1spec'); import evans as E
sys.argv = [sys.argv[0]]
import mod_gram as MG
h = 0.002; x = np.arange(0, 34.3, h); sig = x + E.STIP
P = np.array([E.PROF(s) if s <= 22 else [0, s, 0, 0] for s in sig]).T; p, q, k, m = P
th = -np.pi / 2 + np.concatenate([[0], np.cumsum((k[1:] + k[:-1]) / 2) * h])
np.savez('prof_float.npz', x=x, p=p, q=q, k=k, m=m, th=th)
worst = [0, 0, 0]; n = 0
for b in MG.boxes():
    xc = float(b['X'].mid())
    if xc > 21.9: break_ = True; continue
    i = int(round(xc / h)); 
    e1 = abs(float(b['psi2'].mid()) + np.sin(th[i])); e2 = abs(float(b['cth'].mid()) - np.cos(th[i])) if False else None
    c_pos = abs(float((b['cth']).mid()) - np.cos(th[i]))     # CSIGN None: cth = St C - Ct S
    c_neg = abs(float((b['cth']).mid()) + np.cos(th[i]))
    ep = abs(float(b['p'].mid()) - p[i])
    worst = [max(worst[0], e1), max(worst[1], c_pos), max(worst[2], c_neg)]
    if n < 3 or n % 800 == 0: print(b['tag'], b['ta'], xc, float(b['psi2'].mid()), -np.sin(th[i]), float(b['cth'].mid()), np.cos(th[i]), float(b['p'].mid()), p[i])
    n += 1
print("max |psi2 - (-sin th)|", worst[0], " max |cth - cos th| (cth = St C - Ct S):", worst[1], " vs +cos:", worst[2])
