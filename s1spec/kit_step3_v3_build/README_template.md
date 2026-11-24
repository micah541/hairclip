# Reproduction kit v3 (kit_step3 v3): step-3 certificates of the v3.3 chain

What it reproduces (all python-flint arb with outward rounding, Python 3, no compiler, no network):

| part | what it reproduces | expected time (estimate, 1 core) |
|---|---|---|
| `smoke` | versions, sha256 of every shipped input, syntax, imports | seconds |
| `verifier` | verifier tag **fin12f4** (`FLUX0=1.64e-20`, `FLUX2=8.0e-20`, slopes and F0 x4 against fin11): `tail_cert.py`, `closure.py 1.5e-4` and `1.25e-4`, `mod_sups.py` (sup abar = 2.19103e-07, sup beta = 3.33811e-09) | about 8 min (closure 215 s + 195 s) |
| `ibp` | modulation lemma (unchanged from kit v2) | about 10 min |
| `lemmaSA` | repaired Lemma S(A) (`t2_A_v35`): base chain `run_F10p_ibp -> capture_F2 -> weights_exact -> assembly_arb`, then `hjw -> ic x4 -> a3_signs -> weights x2 -> ramp x2 -> outer_L2 -> assembly x7 -> export_t2` giving `A_interface.json` | about 3.5 min (base chain 2.5 min, t2 chain 40 s) |
| `lemmaU` (alias `zone3`) | Lemma U certificates and zone III, one circular chain: `ramp_U`, `u1_ghat`, `u2_K`, `u34_sweeps`, assemble pass 1, `zone3`, assemble pass 2, `zone3`, assemble pass 3 | about 5 min for u1..u34 + assemble; ramp_U and zone3 not timed [TODO: measure] |
| `lemmaO` | Lemma O (unchanged from kit v2) | about 2 min |
| `text` | `check_keys_v34.py` (keys, sha256 pins, 184 quoted numbers, side conditions of `lemmaU_v34.tex` against the pinned outputs) and `check_balance.py` | seconds |

Everything is run one process at a time. Please do not run other heavy jobs at the same time. `Ctrl-C` is safe; rerun the part.
The new set (verifier, lemmaSA, lemmaU) is about 20 min on the reference machine; with ibp and lemmaO about 35 min. Measure on the first run (`kit_out/timing.json`).

## Install

Linux / WSL:
```
python3 -m venv kitenv && . kitenv/bin/activate
python3 -m pip install -r requirements.txt
```
Windows PowerShell (Python 3.13 from python.org or `py` launcher):
```
py -3.13 -m venv kitenv
.\kitenv\Scripts\Activate.ps1
py -m pip install --upgrade pip
py -m pip install -r requirements.txt
```
Pins that matter: `python-flint==0.9.0` (version of every reference run; other versions may change the last digits, and `t9_lemmaO_v33` differs by about 1e-9 across flint versions) and `sympy==1.14.0` (the kit loads sympy pickles `prop1_forms*.pkl`, `tables_exact.pkl`, `graph_struct.pkl`, `hs_uforms.pkl`).
Unzip to a short path (for example `C:\kit_step3_v3`), not into a OneDrive-synced folder.

## Run

```
python run_all.py smoke            # first: verifies the shipped inputs (sha256), versions, syntax
python run_all.py verifier lemmaSA lemmaU text      # the v3.3 parts (about 20 min); any subset, by name
python run_all.py                  # smoke + every part, in the order verifier ibp lemmaSA lemmaU lemmaO text
python check.py                    # PASS / FAIL / WARN / SKIP table; also kit_out/check_report.txt
python check.py --help             # the comparison rules
```
Use `python3` on Linux. `run_all.py` sets `PYTHONHASHSEED=0`, `PYTHONUTF8=1`, `PYTHONDONTWRITEBYTECODE=1`, single-threaded BLAS and the verifier environment itself; do not set them by hand. Logs and outputs: `kit_out/logs/*.log`, `kit_out/<part>/*`, `kit_out/timing.json`, `kit_out/env.json`.
Parts are independent: `lemmaU` uses the `A_interface.json` of your `lemmaSA` run if present, otherwise the shipped reference copy (it prints which; check.py gives a WARN in that case). `lemmaSA` and `lemmaU` use the shipped (pinned) verifier tables, not your verifier run.

## Order inside `lemmaU` (the chain is circular)

`u1_ghat.py` reads `ramp_U.json` (zone III ramp, `foot_out`), and `ramp_U.py` reads the primary `Y_next` of `u1_ghat.json`; `assemble_U.py` reads a zone3 interface, and `zone3.py` reads the Bhat of an assembled `U_cert`. The kit starts from two pinned seeds (`seeds/u1_ghat.json`, `seeds/zone3_interface.v34.json`, copied into `out/` at the start of the part), runs the steps in the order of `V33_INDEX.md` 4.6 (passes 1, 2, 3) and checks the fixed points (`check.py`: recomputed Y_next = seed; pass 3 = pass 2 apart from args). Details and what could not be confirmed: `BUILD_REPORT.md` (in the builder directory, not shipped in the zip) lists them as `[TODO]`.

## Reading the check table

`check.py` compares with `reference/` (pinned outputs of the proof tree, run with python-flint 0.9.0, hash seed 0):
V1 every JSON of the new parts byte for byte after dropping the provenance `args` block (and reading path-valued strings by file name): PASS if identical; WARN if only numeric leaves differ by <= 1e-9 relative; FAIL otherwise. V2 the 184 numbers quoted in `lemmaU_v34.tex` are re-evaluated exactly on your outputs. V3 the shipped reference files carry the sha256 pinned in `V33_INDEX.md`. Verifier tables, ibp and lemmaO follow the old kit rules R1-R3 (`--help`). Any FAIL is a finding: send `kit_out/check_report.txt` and `kit_out/logs/`.
Headline pins: `A_interface.json` 78fe2f8e..., `U_cert_v35.json` 960809fd..., `zone3_interface.json` 45e8d4dc..., `ramp_U.json` 55317ce5..., step-4 manifest 8bc47e83... (21 stages, shipped as a reference only), `rate_table_fin12f4_a0.00015.json` b20e8715....
The file hashes of the outputs differ from the pinned ones because `args` records paths and script hashes (the shipped scripts are the proof-tree scripts with kit-relative paths); the content comparison above is what counts.

## Layout

`KIT_ROOT` (marker file; every patched script finds the kit root by walking up to it), `run_all.py`, `check.py`, `requirements.txt`, `inputs/` (sha256 manifests), `seeds/`, `reference/`, `step3nl/` (verifier and ibp code), `step3nl_v33/` (pinned fin12f4 table), `writeup/lemmaO`, `writeup/lemmaO_v32`, `writeup/lemmaS/rigor` (pinned tables), `writeup/lemmaU/{phase1,phase2a,bulk}` (lemmaSA tree `bulk/t2_A_v35`, `bulk/t4_cert`, `bulk/t5_zone3`, `bulk/text`), `step4/p2/out_v33` (pinned step-4 inputs `b0_dprof.json`, `n1_chart.json`), `step3`, `lin`, `s1_cap`, `s1_tight`, `cap` (support files of the verifier / ibp, as kit v2).
Scripts that had hard-coded absolute paths of the proof tree are patched (marked `[kit patch]`; paths are computed from `KIT_ROOT`). Nothing in the kit reads outside its own folder; no network, no admin rights.

Not in this kit: the step-4 stages (kit B), the old v3.1 `lemmaS` part (tag fin11), `u5_enum`.
