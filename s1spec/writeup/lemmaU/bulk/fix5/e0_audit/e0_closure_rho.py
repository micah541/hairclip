"""As e0_closure.py, plus the tail threshold rho_min of closure.tail changed IN MEMORY (source text of closure.py read, the default
'rho_min=0.002' replaced, executed as __main__; the file on disk is untouched).  Purpose: with the sharp E0 constants the first anchor with
rate <= -0.002 comes two steps earlier (larger a*), which makes the J's worse; a later anchor uses the better rate.
usage (cwd = output dir): python3 -B e0_closure_rho.py NAME PI3 MODA RHO_MIN ALPHA1"""
import sys, os
V = '/home/claude/s1spec/step3nl_v4/verifier'
name, pi3, moda, rmin, a1 = sys.argv[1], float(sys.argv[2]), float(sys.argv[3]), sys.argv[4], sys.argv[5]
os.environ['MODA'] = repr(moda); os.environ['TAG'] = 'e0audit_' + name
sys.path.insert(0, V)
import constants as K
K.PI3 = pi3; K.MOD_ALPHA = moda; K.ALIN = dict(K.ALIN, sq=59.33, BB=6.22)
src = open(V + '/closure.py').read()
old = "def tail(a, L, orc, log, rho_min=0.002):"; assert src.count(old) == 1
src = src.replace(old, f"def tail(a, L, orc, log, rho_min={float(rmin)!r}):")
print(f"E0 audit closure {name}: PI3={K.PI3} MOD_ALPHA={K.MOD_ALPHA} ALIN={K.ALIN} f572=5.72 rho_min={rmin}", flush=True)
sys.argv = [V + '/closure.py', a1]
g = {'__name__': '__main__', '__file__': V + '/closure.py'}
exec(compile(src, V + '/closure.py', 'exec'), g)
