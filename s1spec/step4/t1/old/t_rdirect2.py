import numpy as np
from base2 import Base; from arbops import SlabPolish; import resid
from polish_chain import polish_chain
from rep import DS1
B = Base('C11'); SPs = polish_chain(B, 10, 3, verbose=False); SP = SPs[-1]
for th_ in (0.0, 0.37):
    sm = resid.samples(SP, th_); o, f = resid.norms(sm)
    X, H = sm['X'], sm['H']
    th = -np.pi / 2 + resid.cumint(sm['k0'], H); G = -DS1 + resid.cumint(np.exp(1j * th), H)
    Gt = resid.cumint(1j * sm['tt'] * np.exp(1j * th), H); nu = 1j * np.exp(1j * th)
    rd = (Gt * np.conj(nu)).real + sm['k2'] - (0.25 - sm['a']) * (G * np.conj(nu)).real - sm['b'] * nu.real
    for xm in (5, 20, 100):
        m = X <= xm
        print(f'theta={th_}: x<={xm}: max|r_direct| {np.abs(rd[m]).max():.3e}  max|r_kernel| {np.abs(f["r"][m]).max():.3e}  max|diff| {np.abs(rd[m]-f["r"][m]).max():.3e}   max|tt| {np.abs(sm["tt"][m]).max():.2e}')
