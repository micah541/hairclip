s = open('gen_hs2.py').read()
old_start = s.index("def ltr(c, w):")
old_end = s.index("def pluecker(c, Wn):")
new = '''def ltr(c, w):
    # scaled frame  wt = (w, w', w'', w''')/rho  of L^T w = lamt w :  d wt/dtau = M_LT wt - P1 wt
    w0, w1, w2, w3 = w
    k, m, q, P1 = c['k'], c['m'], c['q'], c['P1']
    return [f"-{w1}-{P1}*{w0}", f"-{w2}-{P1}*{w1}", f"-{w3}-{P1}*{w2}",
            f"-(-{k}^2*{w2}+({k}*{m}+{q}/4)*{w1}+(1/2-lamt)*{w0})-{P1}*{w3}"]

'''
s = s[:old_start] + new + s[old_end:]
s = s.replace('''        return [f"{c['p']}*{c['R']}", f"C*{c['R']}", f"S*{c['R']}", "wa0", "wb0"]''',
              '''        return [f"{c['p']}*{c['R']}", f"C*{c['R']}", f"S*{c['R']}", f"wa0*{c['R']}", f"wb0*{c['R']}"]''')
s = s.replace('''        return [c['p'], "C", "S", f"wa0*{c['Ri']}", f"wb0*{c['Ri']}"]''', '''        return [c['p'], "C", "S", "wa0", "wb0"]''')
s = s.replace("    v = np.array([1, z, z * z, z ** 3]); sc = 2.0 ** -36", "    v = np.array([1, z, z * z, z ** 3]); sc = 1.0")
open('gen_hs2.py', 'w').write(s)
