"""Regenerate lemmaS_F10p_M0.52.json (args read by rigor/assembly_arb.py l.88) with the tightened (H2) and the core {2,3,3} by IBP (--s33 ibp).
Same argv as rigor/capture_F2.py except --tag F10p.  Run from this directory."""
import sys, os
sys.dont_write_bytecode = True; HERE = os.path.dirname(os.path.abspath(__file__)); os.chdir(HERE); sys.path.insert(0, HERE)
import lemmaS_F2 as F2
sys.argv = ['lemmaS_F2.py', '--J', '6.0273e-12', '8.4673e-12', '1.8708e-12', '--handoff', '2e-7', '1e-4', '0.2', '--sqLam', '5.3075e-8', '--alpha1', '1.5e-4',
            '--core', '0.11', '0.52', '--ramp', '0.05', '0.9', '--mhr', '0', '--slopes', '8.6369e-19', '1.2134e-18', '2.6808e-19', '--rho', '2.693e-3', '--M1s', '0.52',
            '--J4', '6.2362e-8', '--sJ4', '8.9363e-15', '--rampl4', '--rs33', 'ibp2', '--s33', 'ibp', '--l4', '--abar', '2.19103e-7', '--beta', '3.33811e-9', '--tag', 'F10p']
F2.main()
