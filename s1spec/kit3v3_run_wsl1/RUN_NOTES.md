# kit_step3 v3 first run (WSL, Nov 15 2026): 63 PASS, 1 FAIL, 0 WARN, 1 SKIP
- zip e40b444a..., fresh unzip, `nice -n 19 python3 run_all.py; python3 check.py`, about 40 min, Python 3.14.4, flint 0.9.0.
- FAIL `U: ramp_U.json`: only difference is the recorded provenance hash `u1_ghat_sha` (5a06d24c... reference vs 4b86cc91... here): ramp_U hashes the u1_ghat.json it read, whose args block differs in the kit. All numbers identical. Kit fix: check.py should treat *_sha provenance fields of args-bearing inputs as provenance (or ramp_U should read the byte-identical pinned seed).
- SKIP: zone3 intermediate (no reference kept), as documented.
