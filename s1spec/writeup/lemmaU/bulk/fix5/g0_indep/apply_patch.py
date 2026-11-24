# G0_INDEP: applies the review patch to the COPY vp/rhs_iv.py and vp/unc.py (constants.py patched by hand: SRC path + G0 block at the end).
r=open('rhs_iv.py').read()
r=r.replace("        self.ramp = {lev:", "        self.g0 = {n: _c(v) for n, v in k.G0.items()}; self.CABP = _c(_A(k.CAB) * _A(k.G0['cpsi']))\n        self.ramp = {lev:",1)
old="    H = R0p + c.EPS * R2p + c.MU * (I0p + c.EPSI * I2p)\n"
new="""    # ======================= G0_INDEP review patch =======================
    G0F = K.G0FIX
    if G0F:
        g = c.g0; three = cc(3.0)
        Sr2 = mx(mid[2], regS["12<x<=16.25"][2])                 # sup|w_xx| on the inner ramp [10.25,16.25] (rho-Sobolev, regions 8..12, 12..16.25)
        betaH = (g['b0'] + g['b1'] * al) * A + g['bq'] * (_q * _q)   # x-degree>=1 part of the |beta| bound
        betaL = g['bF'] * fn                                     # x-degree 0 part (strip constant NOT included: see report)
        cb = g['cb0'] + g['cb1'] * al
        if 'x' in G0F:    # level 0: cross terms of the twice-integrated top-order term (1-J^-4) v_xxx
            I0p = I0p + (g['xa1'] + g['xb1'] * supv1) * (V1 * V2) + (g['xa0'] + g['xb0'] * supv1) * (V0 * V2)
        if 'b' in G0F:    # beta coupling, levels 0 and 2
            I0p = I0p + cb * betaH * V0; L = L + c.MU * (cb * betaL * V0)
            I2p = I2p + cb * betaH * D2v; L = L + c.MUEPSI * (cb * betaL * D2v)
        if 'p' in G0F:    # -a psi1_x - b psi2_x
            I0p = I0p + c.CABP * Amx * V0; I2p = I2p + c.CABP * Amx * D2v
            if Amc is not None:
                AmcF = cl['am_phi'] * FnW(al, c)        # joint (alpha,x)-degree 2 with V0: homogeneous part
                L = L + c.MU * (c.CABP * AmcF * V0) + c.MUEPSI * (c.CABP * AmcF * D2v)
                if FL: L = L + c.MU * (c.CABP * cl['am_strip'] * V0) + c.MUEPSI * (c.CABP * cl['am_strip'] * D2v)   # strip constant of |(a,b)|: degree-1 FLOOR term (f1-type), NOT in closure.tail's floor formula (size 1e-31, see report)
        if 'm' in G0F:    # level 0: quadratic/cubic monomials of D on the inner ramp
            I0p = I0p + cl['ar0'] * (g['d200r'] * Sr2 * B + g['d110r'] * Sr2 * W3 + g['d300r'] * (Sr2 * Sr2) * B) * Vq
        if 'q' in G0F:    # level 2, far region: missing cross terms / factors
            I2p = I2p + (cl['d001'] * V4 * (cl['t261'] * V3 + cl['t17'] * V2) + (g['cp1'] + g['cp2'] * supv1) * V3 * D2v
                         + (cl['two'] * cl['d200'] * supv1 * V2 + cl['two'] * cl['d300'] * (supv1 * supv1) * V2) * D2v
                         + (cl['two'] * cl['d100'] * V2 + cl['two'] * cl['A1'] * V1 + cl['two'] * cl['A1'] * V0) * D2v + g['e1'] * (V1 + V0) * V2)
        if 'f' in G0F:    # level 2 forcing: alpha- and alpha^2-residual derivative norms (v3 F1b)
            L = L + c.MUEPSI * ((g['f21'] * al + g['f22'] * (al * al)) * D2v)
        if 's' in G0F and K.RAMP == "table" and 'a' in ADD:    # level 2 on the inner ramp (linear deviation + monomials), [est] model
            I2p = I2p + cl['ar0'] * g['r2lin'] * ramp * D2v
            I2p = I2p + (g['d110r'] * Sr2 * (V3 + cl['ar0'] * (g['r3a'] * W3 + g['r3b'] * B + g['r3c'] * W1)) + cl['ar0'] * g['d110r'] * supv2 * W3
                         + cl['ar0'] * (cl['two'] * g['d200r'] * Sr2 * W3 + three * g['d300r'] * (Sr2 * Sr2) * W3)) * D2v
    H = R0p + c.EPS * R2p + c.MU * (I0p + c.EPSI * I2p)
"""
assert old in r; r=r.replace(old,new,1); open('rhs_iv.py','w').write(r)
u=open('unc.py').read()
old="    if K.FIX and K.CLIN > 0:"
new="""    if 'b' in K.G0FIX:   # G0_INDEP: beta-forcing parts  mu cb bF Fn V0  and  mu epsI cb bF Fn (V4 + 2.61 V3 + 1.7 V2)
        cbF = (A_(K.G0['cb0']) + A_(K.G0['cb1']) * al) * A_(K.G0['bF']) * fn
        l['V0'] = l['V0'] + mu * cbF; l['V2'] = l['V2'] + mu * eI * cbF * (A_('1.7') + A_('1.305') * (k2 + C1W / k2)); l['V4'] = l['V4'] + mu * eI * cbF * (1 + A_('1.305') / k2)
    if 'f' in K.G0FIX:
        f2 = A_(K.G0['f21']) * al + A_(K.G0['f22']) * al * al
        l['V2'] = l['V2'] + mu * eI * f2 * (A_('1.7') + A_('1.305') * (k2 + C1W / k2)); l['V4'] = l['V4'] + mu * eI * f2 * (1 + A_('1.305') / k2)
    if K.FIX and K.CLIN > 0:"""
assert old in u; u=u.replace(old,new,1); open('unc.py','w').write(u)
