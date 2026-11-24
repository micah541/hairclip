# probe/ : feasibility experiment only. NOT verifier v4.
`rhs_probe.py` is a copy of `step3nl_v33/verifier/rhs_iv.py` with the terms of EI_SPEC §5 added behind env `PROBE=1`
(`PBOX` = 1e-3 | 2.5e-4 selects the [est] constant set `PROBE_CONSTS`; `PFAC` multiplies all new terms). `PROBE=0` reproduces fin12f4 bit for bit.
All new constants are estimates; the branch-and-bound is the rigorous fin12f4 one. Results are therefore [est]. Known deviations from the spec:
Amc still uses FnW (as fin12f4); T12, Hardy and Q2B terms are omitted (tiny); CF is raised by 2e-3 instead of RQ and most of R5;
the Sobolev form of S1z is the cruder one; `unc.py` was not updated for the new L terms.
It was run in the session scratchpad next to unmodified copies of bb.py, closure.py, constants.py, ival.py, unc.py, tail_cert.py, mod_sups.py.
Files: `rhs_probe.diff` (against the v33 file), `cell.py`, `run_probe.sh`, logs and rate tables of the two marches (x1, x3).
