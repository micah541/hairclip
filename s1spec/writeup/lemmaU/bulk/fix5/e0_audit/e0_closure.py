"""E0 audit: full closure march of the fin13v4 verifier (code unmodified, read only) with E0 constants patched in memory.
usage (cwd = output dir): python3 -B e0_closure.py NAME PI3 MODA ALPHA1     (f572 stays 5.72; ALIN sq/BB set to the rounded-up 59.33/6.22)
Env as verifier/run_v4.sh (MODA env is overridden by the argument)."""
import sys, os, runpy
V = '/home/claude/s1spec/step3nl_v4/verifier'
name, pi3, moda, a1 = sys.argv[1], float(sys.argv[2]), float(sys.argv[3]), sys.argv[4]
os.environ['MODA'] = repr(moda); os.environ['TAG'] = 'e0audit_' + name
sys.path.insert(0, V)
import constants as K
K.PI3 = pi3; K.MOD_ALPHA = moda; K.ALIN = dict(K.ALIN, sq=59.33, BB=6.22)
print(f"E0 audit closure {name}: PI3={K.PI3} MOD_ALPHA={K.MOD_ALPHA} ALIN={K.ALIN} f572=5.72 (unchanged)", flush=True)
sys.argv = [V + '/closure.py', a1]
runpy.run_path(V + '/closure.py', run_name='__main__')
