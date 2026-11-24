"""Monitored run for the Lyapunov constants (angle restarted at sigma = 8), then the constants."""
import gen_hs3 as G3, subprocess, sys
G3.build_A2('inA5.txt')
out = G3.run('inA5.txt', env={"REORG": "1e300", "MONALL": "0.01"})
open('outA5.txt', 'w').write(out)
assert "RESULT OK" in out
subprocess.run([sys.executable, "lyap_consts2.py"], check=True)
