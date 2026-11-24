"""Regenerate lemmaS_F10p_M0.52.json (args read by rigor/assembly_arb.py l.88) with the tightened (H2) and the core {2,3,3} by IBP (--s33 ibp).
Same argv as rigor/capture_F2.py except --tag F10p.  Run from this directory."""
import sys, os
sys.dont_write_bytecode = True; HERE = os.path.dirname(os.path.abspath(__file__)); os.chdir(HERE); sys.path.insert(0, HERE)
import lemmaS_F2 as F2
sys.path.insert(0, '/home/claude/s1spec/writeup/lemmaU/bulk/fix13_v37chain'); import v36cfg as CFG   # fix6_v36 (v36): verifier source + switches V36_VER/V36_R3/V36_R2
sys.argv = CFG.f2_argv('F10p')   # v36: verifier inputs read from the raw JSON, rounded outward (v36cfg.verifier_inputs); fin12f4 = the v35 typed argv (asserted)
F2.main()
