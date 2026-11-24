# BUILD_REPORT: kit_step3 v4

Built by `kit_step3_v4_build/build.py` (2026-11-20 14:14) from the proof tree `<proof-tree>` = /home/claude/s1spec. Nothing certified was run by the build.

**Zip `kit_step3_v4.zip`: 16296685 bytes, sha256 `0f899eb2de74cf51e25aafb9ac0c401343279a9a34f1dbf0dc7a2286809ce853`.** Tree `kit_step3_v4/`: 548 files, 46283949 bytes.

Step-4 manifest the objects are pinned by: `step4/p2/out_v38/FINAL_manifest_v38.json` sha256 `be909bbf42543ab017bb9c81ba58227c80fc4f010cc1b50a66a515b848f90418`.

## Source files asserted at build time

| file | sha256 | pinned by v38 manifest |
|---|---|---|
| `step3nl_v4c/runs/P/rate_table_fin13v4c_a0.00015.json` | `0f2abffba2a98a406ff81b070dfdf8274b501d73722db9d9361494a5449649bd` | yes |
| `step3nl_v4c/runs/P/rate_table_fin13v4c_a0.000125.json` | `6aaeff628dfaf5b2f056030a7feb5deafb2d878ccf89a0e8b8cefaade14a031e` | yes |
| `step3nl_v4c/runs/P/modsups_fin13v4c/final_alpha_points_5e-8.json` | `88d568f336dc6f251ea3a664d710d425307cee9ba514067febeece5cf92b60d9` | yes |
| `step3nl_v4c/v4bcert/out/strip_cubic_v4c.json` | `15ccc9cb20f95ee6ad5447886846aaeadc1914c28a8ab8b759485ce13d2707cc` | yes |
| `step3nl_v4c/v4bcert/out/e0_A2_budget_v4c.json` | `ef7647fa84ab4d985c66fb07dfbad374d446df56457f3d8895616215eba5d95f` | yes |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/out/A_interface.json` | `c48791b2eb78ff7bd7d87944651d726f41029f91cb629a945e97e6fafa967b42` | yes |
| `writeup/lemmaU/bulk/fix13_v37chain/t4_cert/out/U_cert_v36_pass3.json` | `8a3af159936c74c07db3be9b493ac736922b5af1fc2c7858be6fe1072bdb89f0` | yes |
| `writeup/lemmaU/bulk/fix13_v37chain/t5_zone3/out/zone3_interface.json` | `9fcafb6bcbe8da9770218ba319e70f89ee92b66f3ac1b58011abc45b3d1ea600` | yes |
| `writeup/lemmaU/bulk/fix13_v37chain/n1/out/n1_chart_fin13v4c.json` | `459fa058df41fd3190c75fc91116fe46a8272a8a45ad25b1aac5938cb73f193c` | yes |

## Notes

- support trees step3, step3nl (minus its v3.3 verifier), lin, s1_cap, cap, s1_tight, writeup/{lemmaO, lemmaS, lemmaO_v32, lemmaU/phase1, lemmaU/phase2a, lemmaU/bulk/t1_statements}, reference/{ibp, lemmaO} are the patched copies of kit_step3 v3 (unchanged); added step3nl/graph_eq.pkl, s1_tight/prof/prof_tight.txt, s1_tight/ext813/ext813_prof.txt.
- (A) chain v37: t2_A copied without out/, prov/, ref_phase2a/, probe/ and without the 5 base-chain outputs (regenerated); 315 files.
- flux pinned input added (inputs.json pinned): step3nl_v4c/v4cert/c2_t3d.py
- flux pinned input added (inputs.json pinned): step3nl_v4c/v4cert/out/c2_t3d.json
- flux pinned input added (inputs.json pinned): step3nl_v4c/verifier/run_v4c.sh
- flux pinned input added (inputs.json pinned): writeup/lemmaU/bulk/fix7_flux2/FLUX2_DESIGN.md
- flux pinned input added (inputs.json pinned): writeup/lemmaU/bulk/fix7_flux2/flux2_lemma.tex
- grep proof: no .py file of the kit contains "/home/claude" (build.py refuses otherwise). 76 data files contain it inside provenance fields (args blocks, recorded paths); they are parsed as data only.

## [TODO] / limits (stated plainly)

- [TODO] Windows: not run on Windows yet; the WSL test run is recorded in s1spec/kit3v4_run_wsl1/RUN_NOTES.md.
- [limit] The verifier's certification chain (v4cert, v4bcert stages other than strip_cubic_v4c and e0_A2_budget_v4c) is NOT rerun; its outputs are shipped as pinned inputs (the rate tables pin them by sha256).
- [limit] The interface certificates C1_ke and GS_tau1 need step-4 data (t6c jets, base slab 180); they are in kit_step4b v4, not here.

## Patches (38 entries)

| kit file | kind | old | new |
|---|---|---|---|
| `step3nl_v4c/agent_int/v3/c4fun.py` | generic | `2 absolute path literal(s)` | `(_S1 + '/...')` |
| `step3nl_v4c/v4bcert/strip_cubic_v4c.py` | generic | `3 absolute path literal(s)` | `(_S1 + '/...')` |
| `step3nl_v4c/v4bcert/v4bcommon.py` | generic | `12 absolute path literal(s)` | `(_S1 + '/...')` |
| `writeup/lemmaU/bulk/fix13_v37chain/v36cfg.py` | generic | `12 absolute path literal(s)` | `(_S1 + '/...')` |
| `writeup/lemmaU/bulk/fix13_v37chain/n1/n1_chart.py` | generic | `11 absolute path literal(s)` | `(_S1 + '/...')` |
| `writeup/lemmaU/bulk/fix13_v37chain/r3/r3_forcing.py` | generic | `2 absolute path literal(s)` | `(_S1 + '/...')` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_ibp/lemmaS_symbolic.py` | generic | `1 absolute path literal(s)` | `(_S1 + '/...')` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_ibp/run_F10p_ibp.py` | generic | `1 absolute path literal(s)` | `(_S1 + '/...')` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_ibp/rigor/assembly_arb.py` | generic | `2 absolute path literal(s)` | `(_S1 + '/...')` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_ibp/rigor/capture_F2.py` | generic | `1 absolute path literal(s)` | `(_S1 + '/...')` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_v31H2/lemmaS_symbolic.py` | generic | `1 absolute path literal(s)` | `(_S1 + '/...')` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_v31H2/rigor/assembly_arb.py` | generic | `1 absolute path literal(s)` | `(_S1 + '/...')` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/rigor/arlib.py` | generic | `2 absolute path literal(s)` | `(_S1 + '/...')` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/rigor/assembly_rep.py` | generic | `3 absolute path literal(s)` | `(_S1 + '/...')` |
| `writeup/lemmaU/bulk/fix13_v37chain/t4_cert/t4lib.py` | generic | `2 absolute path literal(s)` | `(_S1 + '/...')` |
| `writeup/lemmaU/bulk/fix13_v37chain/t4_cert/u34_sweeps.py` | generic | `1 absolute path literal(s)` | `(_S1 + '/...')` |
| `writeup/lemmaU/bulk/fix13_v37chain/t5_zone3/ramp_U.py` | generic | `4 absolute path literal(s)` | `(_S1 + '/...')` |
| `writeup/lemmaU/bulk/fix13_v37chain/t5_zone3/z1_box.py` | generic | `1 absolute path literal(s)` | `(_S1 + '/...')` |
| `writeup/lemmaU/bulk/fix13_v37chain/t5_zone3/z3lib.py` | generic | `1 absolute path literal(s)` | `(_S1 + '/...')` |
| `writeup/lemmaU/bulk/fix13_v37chain/t5_zone3/zone3.py` | generic | `3 absolute path literal(s)` | `(_S1 + '/...')` |
| `step4/p2/t9h/basegeo.py` | generic | `1 absolute path literal(s)` | `(_S1 + '/...')` |
| `step4/p2/t9h/clam.py` | generic | `3 absolute path literal(s)` | `(_S1 + '/...')` |
| `step4/p2/t9h/h1.py` | generic | `3 absolute path literal(s)` | `(_S1 + '/...')` |
| `step4/p2/t9h/h1_summary.py` | generic | `1 absolute path literal(s)` | `(_S1 + '/...')` |
| `step4/p2/t9h/h1dt.py` | generic | `2 absolute path literal(s)` | `(_S1 + '/...')` |
| `step4/p2/t9h/ift.py` | generic | `3 absolute path literal(s)` | `(_S1 + '/...')` |
| `step4/p2/t9h/s1chart.py` | generic | `2 absolute path literal(s)` | `(_S1 + '/...')` |
| `step4/p2/t9h/s1chart2.py` | generic | `4 absolute path literal(s)` | `(_S1 + '/...')` |
| `step4/p2/t9h/win.py` | generic | `3 absolute path literal(s)` | `(_S1 + '/...')` |
| `step4/p2/t9h/win_slab.py` | generic | `1 absolute path literal(s)` | `(_S1 + '/...')` |
| `writeup/lemmaU/bulk/fix13_v37chain/flux/code/chart4_shell.py` | generic | `7 absolute path literal(s)` | `(_S1 + '/...')` |
| `writeup/lemmaU/bulk/fix13_v37chain/flux/code/f2common.py` | generic | `4 absolute path literal(s)` | `(_S1 + '/...')` |
| `writeup/lemmaU/bulk/fix13_v37chain/flux/code/flux2_coef.py` | generic | `3 absolute path literal(s)` | `(_S1 + '/...')` |
| `writeup/lemmaU/bulk/fix13_v37chain/flux/code/inputs_manifest.py` | generic | `7 absolute path literal(s)` | `(_S1 + '/...')` |
| `step3nl_v4c/v4cert/c2_t3d.py` | generic | `1 absolute path literal(s)` | `(_S1 + '/...')` |
| `step4/p2/v38/G1_check_v38.py` | generic | `3 absolute path literal(s)` | `(_S1 + '/...')` |
| `step4/p2/v38/U1_sym_v38.py` | generic | `1 absolute path literal(s)` | `(_S1 + '/...')` |
| `step4/p2/v38/GS_boot_v4c_v38.py` | generic | `7 absolute path literal(s)` | `(_S1 + '/...')` |

## Data files containing /home/claude (provenance only)

- `reference/cert/e0_A2_budget_v4c.json`
- `reference/cert/strip_cubic_v4c.json`
- `reference/chain/n1/n1_chart_fin13v4c.json`
- `reference/chain/r3/r3_forcing.json`
- `reference/chain/r3/r3_symcheck.json`
- `reference/chain/t2_A/A_interface.json`
- `reference/chain/t2_A/a3_signs.json`
- `reference/chain/t2_A/asm_design_N20.json`
- `reference/chain/t2_A/asm_design_N2ref.json`
- `reference/chain/t2_A/asm_design_N2ref4.json`
- `reference/chain/t2_A/asm_design_pm3ibp.json`
- `reference/chain/t2_A/asm_hi_N2ref.json`
- `reference/chain/t2_A/asm_lo_N2ref.json`
- `reference/chain/t2_A/asm_probe_N2ref.json`
- `reference/chain/t2_A/hjw_rep.json`
- `reference/chain/t2_A/ic_design.json`
- `reference/chain/t2_A/ic_hi.json`
- `reference/chain/t2_A/ic_lo.json`
- `reference/chain/t2_A/ic_probe.json`
- `reference/chain/t2_A/outer_L2.json`
- `reference/chain/t2_A/ramp_rep.json`
- `reference/chain/t2_A/ramp_rep_pm3ibp.json`
- `reference/chain/t2_A/weights_t2.json`
- `reference/chain/t2_A/weights_t2_fb.json`
- `reference/chain/t4_cert/U_cert_for_zone3_v36.json`
- `reference/chain/t4_cert/U_cert_for_zone3_v36_pass1.json`
- `reference/chain/t4_cert/U_cert_v36_pass1.json`
- `reference/chain/t4_cert/U_cert_v36_pass2.json`
- `reference/chain/t4_cert/U_cert_v36_pass3.json`
- `reference/chain/t4_cert/u1_ghat.json`
- `reference/chain/t4_cert/u2_K.json`
- `reference/chain/t4_cert/u34_primary.json`
- `reference/chain/t5_zone3/ramp_U.json`
- `reference/chain/t5_zone3/zone3_interface.json`
- `reference/chain/t5_zone3/zone3_interface_pass1.json`
- `reference/chain/v36cfg.json`
- `reference/flux/out/chart4_shell.json`
- `reference/flux/out/flux2_cert.json`
- `reference/flux/out/flux2_coef.json`
- `reference/flux/out/flux2_data.json`
- `reference/flux/out/flux2_forms.json`
- `reference/flux/out/flux2_sanity.json`
- `reference/flux/out/gen_texts.json`
- `reference/flux/out/inputs.json`
- `reference/proofs/G1_check.json`
- `reference/proofs/GS_boot_v4c.json`
- `reference/proofs/U1_sym.json`
- `reference/step4/FINAL_manifest_v38.json`
- `reference/verifier/final_alpha_points_5e-8.json`
- `reference/verifier/rate_table_fin13v4c_a0.000125.json`
- `reference/verifier/rate_table_fin13v4c_a0.00015.json`
- `step3nl_v4c/runs/P/modsups_fin13v4c/final_alpha_points_5e-8.json`
- `step3nl_v4c/runs/P/rate_table_fin13v4c_a0.000125.json`
- `step3nl_v4c/runs/P/rate_table_fin13v4c_a0.00015.json`
- `step3nl_v4c/v4bcert/out/e0_A2_budget_v4c.json`
- `step3nl_v4c/v4bcert/out/e0_identities.json`
- `step3nl_v4c/v4bcert/out/e0_lemma_consts.json`
- `step3nl_v4c/v4bcert/out/e0_strip_drift.json`
- `step3nl_v4c/v4bcert/out/e0v4b_constants.json`
- `step3nl_v4c/v4bcert/out/e2_consts.json`
- `step3nl_v4c/v4bcert/out/e2_identity_check.json`
- `step3nl_v4c/v4bcert/out/e2v4b_constants.json`
- `step3nl_v4c/v4bcert/out/nl_regions_a2.5e-4.json`
- `step3nl_v4c/v4bcert/out/r2_strip_ibp.json`
- `step3nl_v4c/v4bcert/out/strip_cubic_v4c.json`
- `step3nl_v4c/v4bcert/out/tc_weighted.json`
- `step3nl_v4c/v4cert/out/c2_t3d.json`
- `step3nl_v4c/v4cert/out/constants_v4.json`
- `step3nl_v4c/verifier/run_v4c.sh`
- `step4/p2/out_v33/b0_dprof.json`
- `writeup/lemmaU/bulk/fix11_review/d_cubic_strip.json`
- `writeup/lemmaU/bulk/fix13_v37chain/t5_zone3/out/ramp_U.json`
- `writeup/lemmaU/bulk/fix3/flux/flux_cert.json`
- `writeup/lemmaU/bulk/fix7_flux2/FLUX2_DESIGN.md`
- `writeup/lemmaU/bulk/t5_zone3/out/zone3_interface.json`
- `writeup/lemmaU/phase2a/U_ramp/ghat_design_arb.json`

## Every file (path, bytes, sha256)

| path | bytes | sha256 |
|---|---|---|
| `KIT_ROOT` | 0 | `e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855` |
| `README.md` | 2125 | `df62dead025604985eed3b4cd56c5213dad01bb2da3a9473262b46ffdec7e45e` |
| `cap/gen_evans.py` | 2235 | `73dd6b5d80e32fe98c4175b3cfa8acba727a25758120443eacca6ebeda86e195` |
| `check.py` | 15738 | `9405af602e2519d84acf209fd476196ab910e3e755c11d58df19fe635efdcc0f` |
| `inputs/MANIFEST.sha256` | 59749 | `ecbd977e7022155d6d7846858aaa40b184fc9eb412df5261ffa90d3c124e2e82` |
| `inputs/REFERENCE.sha256` | 7721 | `e4d96d2374f14ec813c0a1a2a88c09dee81af484d0cbc1ab2d947c9d57a72fd7` |
| `kit_guard.py` | 7645 | `970e4e50d7de034944c49143e5f4fdf95532b0b08778c305a38e7c90dc6cf390` |
| `lin/lintail2.py` | 20794 | `e0a488e368e37bfcecd1c0966e1557416a58009d9bccae3f0c7ff987bcfe31d6` |
| `reference/cert/e0_A2_budget_v4c.json` | 2349 | `ef7647fa84ab4d985c66fb07dfbad374d446df56457f3d8895616215eba5d95f` |
| `reference/cert/strip_cubic_v4c.json` | 7221 | `15ccc9cb20f95ee6ad5447886846aaeadc1914c28a8ab8b759485ce13d2707cc` |
| `reference/chain/base_ibp/lemmaS_F10p_M0.52.json` | 4829 | `fada6fba751135f2fd9f374a98b516b1b1e21700a367083f7a6fa3c8fad7b3a4` |
| `reference/chain/base_ibp/rigor/assembly_arb.json` | 11779 | `6c126bcc86aea50f95f921b3ddfb01a47db9f4a9d33bc136974cc1deb48b8a57` |
| `reference/chain/base_ibp/rigor/cap_F10p.json` | 6619 | `cfd1102e07717d299afa2f33c717ff8580f5a254441c4cca6b6ff817b0ed47da` |
| `reference/chain/base_ibp/rigor/lemmaS_rigcap_M0.52.json` | 4831 | `de81eb78511edc35559f4580e30f31f8a6579d566fa8ba7e0c18591c57dcda95` |
| `reference/chain/base_ibp/rigor/weights_exact.json` | 15438 | `0d17e927088056ff37fbb6785fc1d60388b03fa9e2947bf2c9c43c81adca7281` |
| `reference/chain/n1/n1_chart_fin13v4c.json` | 5652 | `459fa058df41fd3190c75fc91116fe46a8272a8a45ad25b1aac5938cb73f193c` |
| `reference/chain/r3/r3_forcing.json` | 4351 | `2ceae90c803e6291a0807d96d6633669976dc1f126e27de568733442da168948` |
| `reference/chain/r3/r3_symcheck.json` | 771 | `fd3322cef2cad38bf3c91b2f1604c234f70d8f06153a98802d21f82b9c058b3a` |
| `reference/chain/t2_A/A_interface.json` | 118594 | `c48791b2eb78ff7bd7d87944651d726f41029f91cb629a945e97e6fafa967b42` |
| `reference/chain/t2_A/a3_signs.json` | 3364 | `e4d9bccb3e03790fd283885b736c9597ebf70294d4452c14972ad56fcf4375e0` |
| `reference/chain/t2_A/asm_design_N20.json` | 132269 | `1af3f008b3f0f1890cf0cac3e64f7803e94302afcad1e2085b3c36405ec475d3` |
| `reference/chain/t2_A/asm_design_N2ref.json` | 133512 | `447d5dfa02cf218a2be2bd981e3e48f83d9f7ab5c6d4d519e0b99818093b84ae` |
| `reference/chain/t2_A/asm_design_N2ref4.json` | 133518 | `87641c798c693e0a3bde2909185d54b4a49f2d6bd3057fc55e6a56e2fe15e822` |
| `reference/chain/t2_A/asm_design_pm3ibp.json` | 132869 | `10cb6bb9f1b32995683cdaf036f5b55dd13217261c6551d98a32d79d3b7368d3` |
| `reference/chain/t2_A/asm_hi_N2ref.json` | 133481 | `874d2f0df9b32d971443acfc64cecdf7050572cd7ed17102e76ba1875a3c8ca5` |
| `reference/chain/t2_A/asm_lo_N2ref.json` | 133491 | `1066369be92b8cc4cf94b5169bccc1cef71311488a243582aeb6d20407c62956` |
| `reference/chain/t2_A/asm_probe_N2ref.json` | 133381 | `4078ffda615c915811ba12b38715fbe9ea62b61f7f9f7dd8f5953075b08e4f7c` |
| `reference/chain/t2_A/hjw_rep.json` | 4069 | `b8c122a523b782faaf3ecc3217fdbe2afac64b06c8b80234fc5cf3f283969a8d` |
| `reference/chain/t2_A/ic_design.json` | 2361 | `040e73c831f22c4200c80c9f970ac0e5e6147d24302c264bcb149f36b12818f1` |
| `reference/chain/t2_A/ic_hi.json` | 2354 | `545b486ffc890f8cb227b27d69f4155095ae27378e3a4299e23baa09a6a2148e` |
| `reference/chain/t2_A/ic_lo.json` | 2347 | `08f60d718ba07ffc6d1d2e50264531ac44c212614e4aa736b70654b67889b485` |
| `reference/chain/t2_A/ic_probe.json` | 2379 | `585fc8540a69c9283f8a07acb8afe2c0981cd37fb4a0df0b5aae8f35f4b04df3` |
| `reference/chain/t2_A/outer_L2.json` | 20507 | `29f82b3986b7f55595ef04e5d49bfc151bb19216a1f87967d04a68c91c30ac0a` |
| `reference/chain/t2_A/ramp_rep.json` | 26125 | `ffb80e1f8c8ce3af76a5c24fbc5149c74e46695d948e4438e4f110f3eb293cf7` |
| `reference/chain/t2_A/ramp_rep_pm3ibp.json` | 24827 | `51789f32a86192d8481ed4de6a08094b91a95e8ec5cabe4916521dcf4e8f7a39` |
| `reference/chain/t2_A/weights_t2.json` | 48803 | `8a37f374feecf0aad0caab6c29eec7bd027271a71617ac483532710097044bd0` |
| `reference/chain/t2_A/weights_t2_fb.json` | 48862 | `5eca84d4a4a42f48553af7a8e08bf4177e8b73c09488b08fb706d9934c4bf39f` |
| `reference/chain/t4_cert/U_cert_for_zone3_v36.json` | 843 | `a553034d0905826e19e559dc1b7bac75659355fa6d034e2713476ae3df887b33` |
| `reference/chain/t4_cert/U_cert_for_zone3_v36_pass1.json` | 849 | `309ae64d9d4fc6328c58e376b64c118fa09efbf480761208ee871b58b768f936` |
| `reference/chain/t4_cert/U_cert_v36_pass1.json` | 35314 | `b21c3634754de88802285d5f9a5c2a03983be84484dc0c7f504e15d4073a8954` |
| `reference/chain/t4_cert/U_cert_v36_pass2.json` | 35369 | `645b3a80ac2fc0427dbd78faa9ee963366d12d7d5752ea45f0ea488495171517` |
| `reference/chain/t4_cert/U_cert_v36_pass3.json` | 35357 | `8a3af159936c74c07db3be9b493ac736922b5af1fc2c7858be6fe1072bdb89f0` |
| `reference/chain/t4_cert/u1_ghat.json` | 26865 | `6f9790f82927e7cd5160443c863e6cfbe834e7b857b949988c1d1ac18545cd4f` |
| `reference/chain/t4_cert/u2_K.json` | 9300 | `d559a3e48bb1e2360053ae2243213e48f09b46ee719116a8f67b7fe5d5e27d8e` |
| `reference/chain/t4_cert/u34_primary.json` | 32363 | `a97046dc9eba82140aeb4cd4e33ceffcd10057671e241eeb006f9282158bc98f` |
| `reference/chain/t5_zone3/ramp_U.json` | 6631 | `55317ce5a31e8e3c3b3c4b86e85ef9569bb6d5e99190b68b0f21e4cba1a0c0de` |
| `reference/chain/t5_zone3/zone3_interface.json` | 18924 | `9fcafb6bcbe8da9770218ba319e70f89ee92b66f3ac1b58011abc45b3d1ea600` |
| `reference/chain/t5_zone3/zone3_interface_pass1.json` | 18940 | `93d46c40275e832c1f65561faf807fd812b5327d3a1aa0f5d4d7922555ee0773` |
| `reference/chain/v36cfg.json` | 619 | `bf2b51de3836149499292334aee88e03905ac83531954c0aa54b7f15d7982a9b` |
| `reference/flux/out/chart4_shell.json` | 11085 | `88307694104cc8eb23925a1af5f0e7ea0f2150c1eb0228c489d401b6f4b8456f` |
| `reference/flux/out/flux2_cert.json` | 15288 | `597da287203590fd42ecf79f0846f8e574e85a10a4f9ba4408f322e94cf0976a` |
| `reference/flux/out/flux2_coef.json` | 19820 | `0bd3798731288c2c20c78e850b42778dfa0f0dd6d57e8ca9b091e798cc42c21e` |
| `reference/flux/out/flux2_data.json` | 7231 | `f54a7774bf94112cbd02210033dfa69a9a1114ab919a66a25eaf03738505b0af` |
| `reference/flux/out/flux2_forms.json` | 4503 | `921b1fdc0027739dd0a39f907f2d7ef6948363f1df9213da0a8de0f77cad4f8b` |
| `reference/flux/out/flux2_forms.pkl` | 5445 | `83c2a91dd65476f8f25ebe3fa29328465759e9c64b42bd4d7ef950cf872c09af` |
| `reference/flux/out/flux2_sanity.json` | 2167 | `95e633074adba3c0b0e6bc30050cf607a615df457c7ae3262222e11f2bd9d06c` |
| `reference/flux/out/gen_texts.json` | 2615 | `965f8d3751478107f7be956db8bdbb749a99be6b165a02499ec5180968e3c412` |
| `reference/flux/out/inputs.json` | 6013 | `5e76e7fae3b106dc937d2271a15cd238880d63877547d038953297bb605538ef` |
| `reference/flux/texts/BPhi_continuation_v36.tex` | 878 | `9a4b159b8ccd85b23c2bb04d572a3ada4fe4fffa81371b51ab609e0ab97c2359` |
| `reference/flux/texts/absorption_lemma_v2.tex` | 2474 | `ff97ad56c3096e26aab95fc9806c895969b4e520f784ac18a6e1d11ac52fa7d6` |
| `reference/flux/texts/erratum_flux2.tex` | 1420 | `f370b2f35747234941f2d85c22086084f0670ca00dda1c66cd022869809d3c0d` |
| `reference/flux/texts/flux2_lemma_v2.tex` | 22030 | `14ef530260351ab298005d6621fe378943f7822ea04284aed7306ef901012490` |
| `reference/flux/texts/lemmaD_k2_sketch.tex` | 2690 | `823aa2d3fd21c7ebdb6123dc1af9dfe68efc3a279ca7c8f823addd65481fa64d` |
| `reference/flux/texts/numbers_flux2.json` | 3714 | `32256e075adfaf7c00a89dd2963073837a7dd5c19af240a23992b29972af0207` |
| `reference/ibp/final_5e-4.json` | 12642 | `39caabfe877457bebdb9ea366663640c738ae0a54980d4e4b172981cfc2717f9` |
| `reference/ibp/run_0.0005_0.02_n8.json` | 6546 | `71f7bc690ed4f397c8a851ff92586659fd495427c4b206c163f3bf1e600b443f` |
| `reference/ibp/verify_symbolic.json` | 123 | `33ab1e675525a33c62011cf9ac88cdb011aa5671fe9059fa008b6120c658c53b` |
| `reference/lemmaO/crosscheck_v2.txt` | 273 | `464e66e9988dca155b7c349585abe7334356afd080157caf599a11a5ef7a1ba4` |
| `reference/lemmaO/lemmaO_v3_consts.json` | 8594 | `b00318e23070730898bbfe6634c45e6bfff16d5c1d1e1fd10207ffe9a392e6dc` |
| `reference/proofs/G1_check.json` | 5934 | `38217e1bbe9ec3ab8f3986c463ab273878295dd65e81e65b3ea4942b6ad94820` |
| `reference/proofs/GS_boot_v4c.json` | 3180 | `b606fb5d283e7bfda1afa8c291754927297cfe8944f7495fd92c8d4b81cd5402` |
| `reference/proofs/U1_sym.json` | 1757 | `c3e7efb80a694145f87a12991c51fa1fe9706e98700a5ec72c7fbecbea9ffbf5` |
| `reference/step4/FINAL_manifest_v38.json` | 278537 | `be909bbf42543ab017bb9c81ba58227c80fc4f010cc1b50a66a515b848f90418` |
| `reference/verifier/fin13v4c_status.txt` | 205 | `014833fb3cb2deffeae5c77584aabfab74865b09ebaf7f61a6a38ab985556209` |
| `reference/verifier/final_alpha_points_5e-8.json` | 3596 | `88d568f336dc6f251ea3a664d710d425307cee9ba514067febeece5cf92b60d9` |
| `reference/verifier/rate_table_fin13v4c_a0.000125.json` | 16252 | `6aaeff628dfaf5b2f056030a7feb5deafb2d878ccf89a0e8b8cefaade14a031e` |
| `reference/verifier/rate_table_fin13v4c_a0.00015.json` | 17488 | `0f2abffba2a98a406ff81b070dfdf8274b501d73722db9d9361494a5449649bd` |
| `reference/verifier/tail_cert_log.txt` | 889 | `dfa9dc17495b5faf9620fe578b3e4542c5bfc5f7261ef5dbcf4b2c6ea2b41741` |
| `reference/verifier/verify_log_fin13v4c_a0.000125.txt` | 3118 | `21ac2739a53c30314ecd46964116a5fe5462bc06d45ae4f8e7aaeb2ccf4b4080` |
| `reference/verifier/verify_log_fin13v4c_a0.00015.txt` | 3290 | `241786ef0b45089dc600b4954d9a608ef2b7a292562cc383168b6602c95be68c` |
| `requirements.txt` | 130 | `78487c095c62b8793f926c14f2662ce7bfed4e13f7b6332337b39d9bc22e658e` |
| `run_all.py` | 18996 | `f093b0a6eb107cc514a996a59069d15d1843fc383d31981e25fea097830051ea` |
| `s1_cap/gen_input.py` | 4112 | `9f3b0c3ef4c671f5282b3991947898cca2ea0a3a1635d032350f937267c8b091` |
| `s1_cap/tail.py` | 4737 | `43a8d545a75f77d2fe3cd9208bb9d93ece64c4e0c0183122462c95c516d9e541` |
| `s1_tight/bracket2.json` | 199 | `c7793f237994a16aa827a0e14ce89f7472bcf2f8e5aaea73f363c0ef25778429` |
| `s1_tight/ext813/ext813_prof.txt` | 12575417 | `49e01a0933e65aaf8ff22b25d87b34dfe62960691ee128c9546646a2feb908aa` |
| `s1_tight/laplace.py` | 2178 | `831dcbeabce197896be8de28004a5be863474a8a698d606c7fc4fbc423c487f9` |
| `s1_tight/prof/prof_tight.txt` | 24321633 | `c55625ae0b1356a367270cd2d3a208cbb5eed9dae61b5187b2987a708fa07c77` |
| `s1_tight/tail2.py` | 3918 | `16ac67433c3bf10acbeceaa2f14352b9773ac0d547a6db2cf08d0d938b0a46bd` |
| `step3/gen_hs2.py` | 16544 | `60dece69b394ef9b8d9ecb5bf50a05c91d3a1cced6fc5dc759f18474717fa155` |
| `step3/gen_hs3.py` | 14011 | `0d309d81c11781126f2eaef52350401d5b12a0b813266dcd63da13c70ce0c1e3` |
| `step3/gram_post7.py` | 8481 | `fce99790077d73f4bcda95e67314869b0b28c8302560afe2fa08380b9d6ad665` |
| `step3/hs_uforms.pkl` | 4310 | `226db25a9f4a858098a3aeb647920f66ba1f013821e75d21b4f5af65b6258823` |
| `step3/hstail.py` | 8683 | `0354b4114dc29e4c65f6f775f7cf525e0aeb934854de37e8ecbf90790aac9f5e` |
| `step3/lyap_consts2.json` | 778 | `4221c83e658bf7d2cdba61951a0762e4ca9987bd287a7c950900a2a8d9289d87` |
| `step3/outA5.txt` | 3625539 | `e1c0bd6c45b60ff5b73c1b3850f040de38a2eb50d4f6703c9cd9933e313fc33c` |
| `step3nl/agent_int/planB/planB_ramp.json` | 4386 | `c8fc5bce02b14cf494ee1bd0d477de4a6b2ff6bbd3ae62834031c3b341b106eb` |
| `step3nl/agent_int/v3/c4fun.py` | 7963 | `f5d83a283278d979a8a9c08be52253b752ea1375d09ff4ead853399938e8da50` |
| `step3nl/agent_mod/gram8_g1s4.json` | 3909 | `f18ce9b2cb2b63c5470560fb66033162a7c626b95c99765e965ae73ff90baa13` |
| `step3nl/agent_mod/ibp/ibp_boxes.py` | 5937 | `311f2ca21804d6c9358b03e4558ada7721016267b70d61e4db0f2cf7a39cbe7e` |
| `step3nl/agent_mod/ibp/ibp_consts.py` | 5744 | `75c9cc65ccd375d89d517ecb3b870459c8830ce054e2f0575eca10dc793912a5` |
| `step3nl/agent_mod/ibp/ibp_final.py` | 7169 | `aab5921e5e43b7a058c3caa03355f809e0b058712696b2cdb18c6d89d52e74bb` |
| `step3nl/agent_mod/ibp/ibp_gen.py` | 289765 | `8adbe20a282428e2d15634734ede43e3f307e7c009b0ddf25a09a4627d05af5a` |
| `step3nl/agent_mod/ibp/ibp_run.py` | 4798 | `120c8b480ebc3509a7fd5fb82f71d7a87bb7172fc5d23fecf64d73de42f2dd5e` |
| `step3nl/agent_mod/ibp/ibp_sym.py` | 7513 | `5ddc6d5248329baad91674faa09bfdc6175b86701ae974bcfacd16bf5e9306e8` |
| `step3nl/agent_mod/ibp/simplex.py` | 2053 | `9ce5df93f3116d50fe459a20590f4709c4df3eec0881662a3a10713df5ec4f25` |
| `step3nl/agent_mod/lemma/lem_integrals.json` | 2111 | `987afb8843489308d04538949c577c6af4e3310d69dbc1fc1f83fb01b01b0403` |
| `step3nl/agent_mod/lemma/lyap_extra.json` | 697 | `73427d2242cbf467c8846b99c706958f037c076da7f5be4def192ead7a7c6917` |
| `step3nl/agent_mod/mod_forms.py` | 7887 | `6f30fae14b19a0bf00472783193923508a81e835187f7404a02c9eec3d24af68` |
| `step3nl/agent_mod/mod_gram.py` | 7693 | `00b0f2f9f0612c4e43061c9348cdc79b45eb001e4b71fa07c532ff3f861afd24` |
| `step3nl/cubic_E0fine_0.001_0.01_0.01.json` | 5358 | `9966342863778614733e7bd29ea24a7ca80b61e78487cb4efbd7c1929c7703a6` |
| `step3nl/graph_eq.pkl` | 19477 | `68b92b8756add344708dd0f2b2c89d3ca879f76d0c346a0dc1ce8fdc29d5352f` |
| `step3nl/graph_struct.pkl` | 1849 | `646b444445b82eb435b094ba7e283a5e619c6cb3bbfa67f495070f9a18f4195d` |
| `step3nl/h2fun.py` | 7618 | `12fb4ca0619e9f015c72ca6965fb3815a5f01f5152ef283cd62b4713df354702` |
| `step3nl/h2md_14.json` | 7736 | `d810075bb58ad65645a07139eaa929ee870c6164f4bc921db357e451683728db` |
| `step3nl/h2res.json` | 670 | `19525dd139af17f9292105b44d22e9ee458020b8b120eab133fceac01fd51ea6` |
| `step3nl/h2res.py` | 17056 | `bed8d6935cc7c576b130c216c84bdebfc3a0bbb6555943eacd1e2129857ea56d` |
| `step3nl/psi3fun.py` | 7403 | `f94b652ee33a3a54a88d5339ab04e0034eb8c5d49fa78c7efe2d7c958092b019` |
| `step3nl/psi3ip.json` | 422 | `a8cd8f3b1268e7c05dc0d80fa1130980aa70117a05ef87ca5ff40b5c1b106704` |
| `step3nl/psi3md2_14.json` | 7784 | `f0b3348aace69762cb63dee523ac520aacae7246c66aa0571b7e6073340ef2a5` |
| `step3nl/q2_derivs.pkl` | 1762 | `f9dd27d690049910e2a7320e3e9c3e761d063cac2d2eccd40b5efb65a8845dfc` |
| `step3nl_v4c/agent_int/planB/planB_ramp.json` | 4386 | `c8fc5bce02b14cf494ee1bd0d477de4a6b2ff6bbd3ae62834031c3b341b106eb` |
| `step3nl_v4c/agent_int/v3/c4fun.py` | 7986 | `393a6f276788e4a267deb6aa45f660edb17cebe42a648aa7f924656ef5436827` |
| `step3nl_v4c/agent_mod/bloop/bloop.json` | 3244 | `ece9bd1eb8a795eab48ae3b66e183095813a85df67376fa32f4968ed93520a29` |
| `step3nl_v4c/agent_sob/sob_constants_P2.json` | 2976 | `8412627f71b7746b822dc6a1e14b82fec186a96fc62442fb32de591923cca0e6` |
| `step3nl_v4c/cubic_E0fine_0.001_0.01_0.01.json` | 5358 | `9966342863778614733e7bd29ea24a7ca80b61e78487cb4efbd7c1929c7703a6` |
| `step3nl_v4c/graph_struct.pkl` | 1849 | `646b444445b82eb435b094ba7e283a5e619c6cb3bbfa67f495070f9a18f4195d` |
| `step3nl_v4c/runs/P/fin13v4c_1.25e-4.out` | 3118 | `21ac2739a53c30314ecd46964116a5fe5462bc06d45ae4f8e7aaeb2ccf4b4080` |
| `step3nl_v4c/runs/P/fin13v4c_1.5e-4.out` | 3290 | `241786ef0b45089dc600b4954d9a608ef2b7a292562cc383168b6602c95be68c` |
| `step3nl_v4c/runs/P/fin13v4c_status.txt` | 205 | `014833fb3cb2deffeae5c77584aabfab74865b09ebaf7f61a6a38ab985556209` |
| `step3nl_v4c/runs/P/fin13v4c_tail.out` | 889 | `dfa9dc17495b5faf9620fe578b3e4542c5bfc5f7261ef5dbcf4b2c6ea2b41741` |
| `step3nl_v4c/runs/P/modsups_fin13v4c/final_alpha_points_5e-8.json` | 3596 | `88d568f336dc6f251ea3a664d710d425307cee9ba514067febeece5cf92b60d9` |
| `step3nl_v4c/runs/P/modsups_fin13v4c/mod_sups.out` | 130 | `eed1f90c44e2578570d1c88cc26f4a8873b7c5e3c89ac5b0e9025e6143519d30` |
| `step3nl_v4c/runs/P/rate_table_fin13v4c_a0.000125.json` | 16252 | `6aaeff628dfaf5b2f056030a7feb5deafb2d878ccf89a0e8b8cefaade14a031e` |
| `step3nl_v4c/runs/P/rate_table_fin13v4c_a0.00015.json` | 17488 | `0f2abffba2a98a406ff81b070dfdf8274b501d73722db9d9361494a5449649bd` |
| `step3nl_v4c/runs/P/tail_cert_log.txt` | 889 | `dfa9dc17495b5faf9620fe578b3e4542c5bfc5f7261ef5dbcf4b2c6ea2b41741` |
| `step3nl_v4c/runs/P/verify_log_fin13v4c_a0.000125.txt` | 3118 | `21ac2739a53c30314ecd46964116a5fe5462bc06d45ae4f8e7aaeb2ccf4b4080` |
| `step3nl_v4c/runs/P/verify_log_fin13v4c_a0.00015.txt` | 3290 | `241786ef0b45089dc600b4954d9a608ef2b7a292562cc383168b6602c95be68c` |
| `step3nl_v4c/v4bcert/e0_A2_budget_v4c.py` | 4619 | `c91cda33bd605d9482a0b5c8f02c98e9b5e9ca3c2d0caab06cb24e27dc4b55c8` |
| `step3nl_v4c/v4bcert/out/e0_A2_budget_v4c.json` | 2349 | `ef7647fa84ab4d985c66fb07dfbad374d446df56457f3d8895616215eba5d95f` |
| `step3nl_v4c/v4bcert/out/e0_identities.json` | 993 | `d35b5fbe0e3be6043131ec342f0acd24cf648b4525c5f49513ed08d75b5c825b` |
| `step3nl_v4c/v4bcert/out/e0_lemma_consts.json` | 13863 | `34bb7ce730b3702bf0339c778b16eef78a5972759d7d9bb8a0bee593b238a708` |
| `step3nl_v4c/v4bcert/out/e0_strip_drift.json` | 3468 | `4ee153b4350df880bdc0d6240909cfc6fb801f67bca77fdaab2b765a382185a3` |
| `step3nl_v4c/v4bcert/out/e0v4b_constants.json` | 5106 | `c209638bc8ed34655076562d275142d0c84cf09195fcfbe680171cba3323297f` |
| `step3nl_v4c/v4bcert/out/e2_consts.json` | 8797 | `486c203259189e3164c644d54b66043476ea8be80798223062e95252acf9b288` |
| `step3nl_v4c/v4bcert/out/e2_identity_check.json` | 769 | `3ac1369a600f890bd0ee8f77eb41f485ddea9a41967d8f36e3df340d0eeeddae` |
| `step3nl_v4c/v4bcert/out/e2v4b_constants.json` | 5555 | `7d6ada650e1845f890640f96455b136311ffc5eb747d5a4bcb5b9a63f834bd88` |
| `step3nl_v4c/v4bcert/out/nl_regions_a2.5e-4.json` | 4141 | `bdbe8094a6278a936f69922834e889ee3d7648380fe472e355c4c57d9405fcc7` |
| `step3nl_v4c/v4bcert/out/r2_strip_ibp.json` | 4856 | `0daa73c4462b24c29106cab6176b5eeea01ee2a6380c65fd6ba854b78100126a` |
| `step3nl_v4c/v4bcert/out/strip_cubic_v4c.json` | 7221 | `15ccc9cb20f95ee6ad5447886846aaeadc1914c28a8ab8b759485ce13d2707cc` |
| `step3nl_v4c/v4bcert/out/tc_weighted.json` | 1226 | `e93f83724c4b1076dd17fd8770867e9accd4f810b0b1b9e8554ae2c877d13512` |
| `step3nl_v4c/v4bcert/out/tc_weighted.pkl` | 3794 | `da22cb16c2de16736d6ecacedba9bfc84bb0a4922bd2c0aad97e99dbf64f4ecd` |
| `step3nl_v4c/v4bcert/strip_cubic_v4c.py` | 8341 | `fc2ad9f0c86dd03402939c26da2cebf78860458c01edca118913cd95a0c4f531` |
| `step3nl_v4c/v4bcert/v4bcommon.py` | 2184 | `168065beaf5333b24d99a3cd3b85bf4d31a9ddd01301b1d10f85db405cf4e2a6` |
| `step3nl_v4c/v4cert/c2_t3d.py` | 13059 | `8668f46c2adaa692f07b53c0e62370e38749969fd7b8183f983e529695464995` |
| `step3nl_v4c/v4cert/out/c2_t3d.json` | 33734 | `3051e721189ad17e9905a55e9bcd761c7d7cb0c14e58667d51117e8782b470aa` |
| `step3nl_v4c/v4cert/out/constants_v4.json` | 60756 | `590102ee0a6980051e244ffa94e2a90fec96cf1c7ded7c34a495d51b4e4e6fd2` |
| `step3nl_v4c/verifier/bb.py` | 10569 | `76bc1e24d787f53bf97333af6b894ad67456342114aa001404f05de498d47352` |
| `step3nl_v4c/verifier/closure.py` | 24807 | `67583dc7da0335d62f3465d06b6aa739909e5e90ad959bf537a16347ba6b08e5` |
| `step3nl_v4c/verifier/constants.py` | 14644 | `3a5ca5a43aaf4d77e97ac8e0be4c0468ffd78d70aea6a3d0bc6fd3fb2233f9a9` |
| `step3nl_v4c/verifier/ival.py` | 5669 | `74b65bc2505a30cddd5db0eff1ecc9abca355651bf76743a60d293087a3f9bc9` |
| `step3nl_v4c/verifier/mod_sups.py` | 4745 | `61da77f116ee9b4414abc4a7eaddd71557c04f192742b11611859932f277942a` |
| `step3nl_v4c/verifier/rhs_iv.py` | 31665 | `725faa162d7ae9d3ab5941cb05f40edb4dc344a8cb72da3c60db3049d4ed5a53` |
| `step3nl_v4c/verifier/run_v4c.sh` | 1313 | `698ea177156cf74ff978beb3c96dbf88ef0fe6d9588822d828d9fb1ef4c9ceb7` |
| `step3nl_v4c/verifier/tail_cert.py` | 3871 | `b346e6369e25f6bbc804b192a6540c5fd82d5cf5b89dc15d20dcff8b0bea9094` |
| `step3nl_v4c/verifier/unc.py` | 4604 | `b9b9a41a20053304438f6861668f2824d910d68f626a6934224c7d2f5d55db83` |
| `step4/p2/out_v33/b0_dprof.json` | 7827 | `7fdf378c5f59f098263aec8f39bfff7f365e0d25cbbcd154c3f1010bd4a2cd8e` |
| `step4/p2/t9h/basegeo.py` | 3158 | `e5ef6401fc46a2084c91b1804547d11421e93eb5a7ba1a14ebf9d435f99470c6` |
| `step4/p2/t9h/clam.py` | 4089 | `f2c0e662190135a17d7a8f704356aeccfdf6e3e6185871797b9acbc62965dae5` |
| `step4/p2/t9h/h1.py` | 19905 | `b90e20acdce98ae07dc0a0a0774a386eb7dbab18800821e220eef4aa6c560a51` |
| `step4/p2/t9h/h1_points.py` | 3468 | `23e58ee5563812df8a860a46ec01dd757c6a82a4bfcccdf7a447e9847152e0bb` |
| `step4/p2/t9h/h1_summary.py` | 3066 | `92c69d705fbd5981c383818cccdaa6cdf3a5873f0f5b10b001766ba31489f0c3` |
| `step4/p2/t9h/h1dt.py` | 11955 | `457ea8241e3dbfabafa20bc8a19b208df26d862cdd6191f65647cbad3459d39b` |
| `step4/p2/t9h/ift.py` | 8469 | `229b31171cb6521a0d46384b52ae789b2589a2afcdb1ab1bb641872ab9fbfd48` |
| `step4/p2/t9h/s1chart.py` | 6933 | `6f49900bd3177c4626479c04c6101938e0d4887831c105f10ed3e56b2821c0f0` |
| `step4/p2/t9h/s1chart2.py` | 7620 | `b9f27979a0a5f24afb09c78aa85d02f93963b2a916dc71e3a0841e7d37efc33f` |
| `step4/p2/t9h/win.py` | 6733 | `2ad3f9320893d2473e36fd37fc0c98cfbe0600743b1e2ba14c6b2ef5e945704e` |
| `step4/p2/t9h/win_slab.py` | 2043 | `9b28025beb738fb60aa1c193cd0b871080784b69df6cd778f0eae823c37dca96` |
| `step4/p2/v38/G1_check_v38.py` | 11469 | `b04118d86c48dd4c2814c3919cd2fe64bfdd947f65880c5ea359a8620fd85278` |
| `step4/p2/v38/GS_boot_v4c_v38.py` | 7141 | `7d99ce719a8479c4e431483637a90d48fb5861d9cf9c17c088c7888c6488c2f9` |
| `step4/p2/v38/U1_sym_v38.py` | 9674 | `6b483e66362f53f01c07c610d511ac4c296d5a41b4e1a462fa66925478981836` |
| `writeup/lemmaO/lemmaO_consts.py` | 5729 | `8b252a5ef2301698b0e49491fab1d43da18930f7c43396f2e6c087f17025c1ef` |
| `writeup/lemmaO/lemmaO_consts_v2.py` | 3973 | `7eb9af4a6709ccf5eb925c7a82213e7ff5be4ff8eb97bbcb2b9a80e438a7c88a` |
| `writeup/lemmaO/lemmaO_far_margin2.py` | 946 | `897d5c05595ba29757d68b0015fbbc3d1fc9adabf71cd9f615ae2caae817d16a` |
| `writeup/lemmaO/lemmaO_v3_interval.py` | 9212 | `0a1a799cc98ee38e896eea20dce4b24851ad81d251d6f6885b82ed6b31084d0b` |
| `writeup/lemmaO_v32/lemmaO_v32_consts.json` | 9493 | `c64ae6b534f69aac8fd028d5494a39579a9b7b506b714bc52742fe45dc327f9c` |
| `writeup/lemmaO_v32/lemmaR_consts.json` | 1982 | `56b67482b972e8ff9d83565d2d4586c5cef230a6c22d4f4f237713e3329a4820` |
| `writeup/lemmaS/lemmaS_cflux.json` | 1235 | `f08e90fb25d0d8893a60109827334dc24f53fc37ed0cdf3f6e6d5407ca61bd65` |
| `writeup/lemmaS/rigor/assembly_arb.json` | 11745 | `c8944912a74f38637612fbc6275d56238757464e3baf1a494cdddb49b60ec7e8` |
| `writeup/lemmaS/rigor/assembly_arb.out` | 11905 | `25631ad86bd8c338bc82652d8c044bc78f4dcdba557654daf5a0c53fe327ee88` |
| `writeup/lemmaS/rigor/tables_exact.pkl` | 2839 | `861e3ac8b0bab6edec62221b299e3f80d543a8263fc83baa7fe6804c30e0ce17` |
| `writeup/lemmaU/bulk/fix11_review/d_cubic_strip.json` | 3937 | `6a71e36d2c2ba3974885ccbe64f8bc85ee3bbe06c6093a13df98044f7200ca42` |
| `writeup/lemmaU/bulk/fix13_v37chain/V37_LOG.md` | 2057 | `c36c72c9c1417cad7591f686559a109c8b5f966c9803cf7a21cf7adc2ef0a701` |
| `writeup/lemmaU/bulk/fix13_v37chain/flux/code/chart4_shell.py` | 16319 | `81418ebb2e200731e3b219f6065bd3edd002ac64a5518c08c6d8a7fe7933b282` |
| `writeup/lemmaU/bulk/fix13_v37chain/flux/code/f2common.py` | 5256 | `9dc80e4b1c333c9e12b302b523bf86f4c53cb731e023b3cc2a25f8da3d5bf312` |
| `writeup/lemmaU/bulk/fix13_v37chain/flux/code/flux2_assemble.py` | 15561 | `4ea7671661a1c8a332b8bc388232fe9359c49f0d17679e8870174f36e7d3d287` |
| `writeup/lemmaU/bulk/fix13_v37chain/flux/code/flux2_coef.py` | 22171 | `3f4d6940c09e415f4e0c0028e667b5badc5028d8836a948ea718a66c8894d8fa` |
| `writeup/lemmaU/bulk/fix13_v37chain/flux/code/flux2_data.py` | 12068 | `882e7c264673eb914c824a48b130699103b169463964da315afe5bd3dc401ac9` |
| `writeup/lemmaU/bulk/fix13_v37chain/flux/code/flux2_forms.py` | 17587 | `74b94fdac8e0a11a095b980099ab4cee913cf96c6a1e9122539caf4c1fbafcac` |
| `writeup/lemmaU/bulk/fix13_v37chain/flux/code/flux2_sanity.py` | 6616 | `231674c62c2796438f2fe896e26e75614c0dbd49bbfecdedb007addd3530e6b6` |
| `writeup/lemmaU/bulk/fix13_v37chain/flux/code/gen_texts.py` | 15608 | `6f2d8e87ccebc6337aab6c1a7e4a87e0c9d890902387d8dabc5066e16ef6e5d1` |
| `writeup/lemmaU/bulk/fix13_v37chain/flux/code/inputs_manifest.py` | 4889 | `93145590128947500c647029fb8c4b41e336e27ae2237be5d99e257f6f8d35c8` |
| `writeup/lemmaU/bulk/fix13_v37chain/n1/n1_chart.py` | 17044 | `fe27133f73dc32ce2964697f87edb26792ab8f692f1c5fca30bbc95e6592fe4e` |
| `writeup/lemmaU/bulk/fix13_v37chain/r3/cmpjson.py` | 1617 | `961ef418efa852ebf9716efcca9fcc9893e3ae0fa61649d9be7495d00764cad5` |
| `writeup/lemmaU/bulk/fix13_v37chain/r3/compare_r3.py` | 5306 | `b1ea8a103efff14bc4ecfa76cc88a5bce460ebeeb009acfac570329e70841f36` |
| `writeup/lemmaU/bulk/fix13_v37chain/r3/key_table.py` | 6731 | `ae8bb9e7b0eeac80d546a1f5ed32a59efd85bde61a97be33fb67e35cbf346908` |
| `writeup/lemmaU/bulk/fix13_v37chain/r3/r3_forcing.py` | 10952 | `bc76aacb1f0c3e75e6dc2c4754740fbb86f2c8ad275a8d6252fe1200ecc29a2b` |
| `writeup/lemmaU/bulk/fix13_v37chain/r3/r3_symcheck.py` | 2893 | `44218da9383b6310bae7cc526b9ca7907fdf78d6fead14b3a7f49a5c4fc9b509` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_ibp/F2_ibp.log` | 951 | `7c0b9f6a5aec9ca8239b2021709b7d6666c341cc423a3372abc000787b5a0cac` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_ibp/F2_v33.log` | 951 | `7c0b9f6a5aec9ca8239b2021709b7d6666c341cc423a3372abc000787b5a0cac` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_ibp/F2_v36.log` | 951 | `5a96d157bbc035a71156f8d1f5fb440adf29dc5116bd6cf2e76ac7c3a83e59c9` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_ibp/cflux_forms.pkl` | 4517 | `a862bc04ad6273a8d348029ca455af24367ab4e1d74e40f5f6071c43e2ab47a5` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_ibp/lemmaS_B.json` | 666 | `b38f615bc26567d410a53b89aba8e93362cbc2b73b66bc283f240b840daeb2b9` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_ibp/lemmaS_B.py` | 1477 | `38ff9f934df8cbe70b5158e80f9f822b438a5362ed887941fdf06580d360d7a9` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_ibp/lemmaS_B3.json` | 235 | `a67d2ff97d151c76e37c91302011d84440193bdb954f86a5dc747fae0689b075` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_ibp/lemmaS_B3.py` | 2147 | `71637489053d65219af0d5f2996c2f52025a4d7e33f4d64ed227b056053a3cab` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_ibp/lemmaS_BR0_M0.6.json` | 4120 | `6b01aac83f96d434ca202d6d9da4983f265240c6938d21c4f049b4ea69e17287` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_ibp/lemmaS_BR14.4_M0.6.json` | 4121 | `ec2552abd5a0e0eb273126f95612a02139eba2a63cb819dc16f9c3677a5f5685` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_ibp/lemmaS_BR7.2_M0.6.json` | 4120 | `a3fa4d5c874fd006e95204dfc7144e2500aa6333db5c6a794241f25bce95ef44` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_ibp/lemmaS_Bbox.py` | 4031 | `9c0fb29f05da5bf0d9102a90e3665c377c9e7cc7dd0dc6b41b9a92ffd0e54bb3` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_ibp/lemmaS_Bbox_1e-05.json` | 3543 | `b2181cae4ea6d8943b2a03896e692048dc5a56b29d01d711f64946c4c99a0142` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_ibp/lemmaS_Bbox_1e-05_a0.000125_dv0.15.json` | 3556 | `cc60b6d95c33240b688ce943cfaa7d5ba6da8b194e4ee7b2fd73b321028c1cd2` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_ibp/lemmaS_Bbox_1e-05_a0.00015_dv0.15.json` | 3554 | `b68a1153af35561088edf8221d343b07bc8f9ca498919b43d6f42292f9ef397a` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_ibp/lemmaS_Bbox_3e-05.json` | 3530 | `bcaa321b78a3064aa502079cf6a4081d58583db823be117b59e499325975f62b` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_ibp/lemmaS_Bbox_3e-06.json` | 3546 | `83ddab479214606fd42350be17d4cb69758d11dfbfc47fa345174783c3f9dccf` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_ibp/lemmaS_Bbox_3e-06_a0.000125_dv0.1.json` | 3571 | `53c8a4b73fa533f4ce49ca3e3ba4ab86762f69c33cf3cfac625bc88503a7b00b` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_ibp/lemmaS_Bbox_3e-06_a0.000125_dv0.15.json` | 3559 | `1e3972ca0a1a5b21e388baa78dafb3a921b4221554c3cde2b6be5cba09df9b41` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_ibp/lemmaS_Bbox_3e-06_dv05.json` | 3549 | `0259a710cb5dbd8db4cbc7c9d4f6a76c61983a5702fb434fcf409ebdd10cc0f6` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_ibp/lemmaS_F1.json` | 21340 | `6c6a13762160891e6d1c1999d3b44959198521076cb32b9cfa8724b7e459d240` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_ibp/lemmaS_F1.py` | 8232 | `08cc332aeddd4bb6424579e6d6482abb3545408f4002a978731d2fa1cb739937` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_ibp/lemmaS_F10_M0.52.json` | 4828 | `c51da3461225199d54c0da2285529e7762c336d9d400bfbad924506cc4c69091` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_ibp/lemmaS_F10a_M0.52.json` | 4780 | `81916706a1202a06ad078b5128d2b2e8819005301215353f461fcf044d1b5da9` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_ibp/lemmaS_F10b_M0.4.json` | 4792 | `daa9d13765e172447bcadac08dad9f72ae11e98c361548e1d46c97e79b74428f` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_ibp/lemmaS_F10b_M0.45.json` | 4797 | `06fac3d2b8188a952297e601e23c2d6add6b712acb9f48a144e6443d61635c9c` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_ibp/lemmaS_F10b_M0.52.json` | 4796 | `367ee9ed19e78a0a3da78331eaa0c57ccbc8f631ae31ee6b65437993842d6770` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_ibp/lemmaS_F1_scan.json` | 2807 | `2ff77d1f611ce5e5a73463b61dca08e20ae618cff1db8896b917d6318069ad8d` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_ibp/lemmaS_F2.json` | 3960 | `30210a0dba8cb19c100131ec407c5556388d22a8dd596d9932e219d795a5140a` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_ibp/lemmaS_F2.py` | 13823 | `cb736fdbf6a81dee286448971f385c521815e947a87d7eb267c5789e4a1b9082` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_ibp/lemmaS_F4.json` | 1473 | `1d1ee5f105b1020aa99f17f3175d19f59274d666165d46591641639d5da5ac10` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_ibp/lemmaS_F4.py` | 2128 | `dc1737a8ac4611c035539e0a5e9280c50c91e7a98af9290bf3d156fa9c158b7b` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_ibp/lemmaS_F4B.py` | 5297 | `ec23c2ba0c34f47b8abef72e8f38f7f27e8bfa98df46bfc61b086cf5932f6094` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_ibp/lemmaS_F4B_F10_M0.52.json` | 1170 | `362968b1083de0edc2533f771c43025d57f3ac90af3f5d5dd188a428da7b4b99` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_ibp/lemmaS_F4B_F10p_M0.52.json` | 1169 | `7615261d270722285e6dafb55d080038b8c75ebb489eb03c8330eda49d0f875e` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_ibp/lemmaS_F4B_F9_M0.52.json` | 1164 | `da530d0d66555979b1200b7f516605d926103e38ec1fb00ad9f2943ab78d6423` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_ibp/lemmaS_F4B_L4b.json` | 1228 | `f546d1b5fc1e5dfbbba37901108d468103487c3ea8e52c0c4096fd4fab552a07` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_ibp/lemmaS_F4B_N5b_M0.52.json` | 1157 | `d09c89ebec6bf6b91ed341b2534be7fcab73ac7cd858b18af9711c3b39d622c9` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_ibp/lemmaS_F4B_P1_M0.6.json` | 1222 | `d535a0086e0b594ef52348428839bb25721967411c219d6e52b502bb9eeb0d10` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_ibp/lemmaS_F4B_P2_M0.6.json` | 1227 | `99dfa5cfc5fb728d3cdb12a19665adbb91020f0c23a07007dfd4ea7710006511` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_ibp/lemmaS_F4B_R5_M0.6.json` | 1170 | `7ff52e13bb6850c3ebf2780d20613de0e0f1060924de357f846e8f9051f75908` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_ibp/lemmaS_F4B_R8_M0.52.json` | 1165 | `a933f761039a11cc893183f819b42dec091a6c47a725dbe2b3baf04fe57de3fb` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_ibp/lemmaS_F5_geom.json` | 733 | `61fbffa478db8f0b9e0d8cc67a2d0501e827b03f3d9ecabe406cf0026b419427` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_ibp/lemmaS_F5_geom.py` | 898 | `c440953010a8d4b8a66aaf65bd3f17f5a47d01cf7d257dded0549543a7ae6e02` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_ibp/lemmaS_F9_M0.52.json` | 4757 | `a60f0d3edd9f870fa2a41dd37f145acb543dcd6c35c340c3294bfd307edf8bc4` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_ibp/lemmaS_GN.py` | 2241 | `f459d78489fc7358e5d8b60c8481460955a8816ea47bb412290f0460bfb6d58e` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_ibp/lemmaS_H6.py` | 1111 | `71b6ac231cf07e43b9faa6da5cf39865c6b09cc2fef3c629dd4a9b6738e6a415` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_ibp/lemmaS_L4.json` | 3988 | `f0d39a54b7c264f81df91b8444ddbfb0dae9bb210111ff483cafa8aa9125f266` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_ibp/lemmaS_L4b.json` | 3998 | `35d9325783d9c4cdaa84aaa285caade2c91ccec9ecbd0b74f7cc052f65fa6c39` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_ibp/lemmaS_L4terms.py` | 8708 | `65ecceee1b3f1274749cc0c4a5cab6592a913ab4db8e16ff57cf553e4be782a5` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_ibp/lemmaS_N1_M0.4.json` | 4590 | `a9d48e139c9b07d7075f6bf5c924effd67ae043703530281eb2670f98b8cad84` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_ibp/lemmaS_N1_M0.5.json` | 4602 | `23a9fd3f80b5fb2acf296af21d2c92522146554e08c9f3adb6f59c76e7751b02` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_ibp/lemmaS_N1_M0.55.json` | 4601 | `198a7cfe906947d3a43ace94c2cce7c1a5d6b1c7513fbaa1951ea591f1627621` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_ibp/lemmaS_N1_M0.6.json` | 4591 | `82ea3ccc91ef2eb84d55fb09904659720eb5d68836fb8b1ef31064f11bffa7ff` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_ibp/lemmaS_N1_M0.65.json` | 4597 | `d73dfb510d000aed9a83b4e79ec8b15a1bd8276947d03e9a7aa28fc387043da7` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_ibp/lemmaS_N1b_M0.575.json` | 4570 | `9aa6d22c2776defda7176841c777d52aa2096e3a63af7ed799316b244aa09318` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_ibp/lemmaS_N5_M0.3.json` | 4599 | `f90ddde1005847ccd7658d4dfd7b6b15d67b09ea616ff274bfe9a1bfee3c6ad2` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_ibp/lemmaS_N5_M0.4.json` | 4595 | `4b757b5ead94125cafe18ba6bbb67bef3603599e5f1b8bd3a66db81b54957c08` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_ibp/lemmaS_N5_M0.5.json` | 4607 | `32e9a362b6f561a58f9881cf5db4ee8c1a67251d24a58852508558a3619de879` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_ibp/lemmaS_N5_M0.6.json` | 4598 | `ef9e401d5f47fa70b1a3c6149d571e4497c10cd8532058ef05a82ddf3c689f80` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_ibp/lemmaS_N5_M0.65.json` | 4604 | `f6bd656b12e33eee527f0560f085b8948b1cbd1b1406094d3afaa33e8b8eb7da` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_ibp/lemmaS_N5b_M0.52.json` | 4570 | `745bbdaeb6eb388414ee9059d93e7ab26ec895334c519530dae111617289036a` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_ibp/lemmaS_P10000_M0.52.json` | 4594 | `7a730960b771ceba751763b4a45d8a12a050195bc8daaa0f05134e67c6ee0c89` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_ibp/lemmaS_P1000_M0.52.json` | 4595 | `afbd50f30d5e78674a304905e6e61f3699ebef8eec9056b6d8621ae8d142909e` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_ibp/lemmaS_P1_M0.5.json` | 4140 | `ff476c422aed01e019e72794b9022f0ba730ad0327379381b70983c2a01dbbf1` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_ibp/lemmaS_P1_M0.6.json` | 4134 | `f764274adaae6995eca3d94cda27442c85e8eb4f8eb96617d0614ac3613ae0ad` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_ibp/lemmaS_P1_M0.75.json` | 4143 | `29780ffdc2396751bbb8aa0aafb932bbd12353b23401f29622fcbc63a0a94ed4` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_ibp/lemmaS_P2_M0.5.json` | 4135 | `ad10632395e5840eccc02c3eabd0a7f106da755eba0b3fd77261ece3e73557c1` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_ibp/lemmaS_P2_M0.52.json` | 4772 | `62e1247f9cdf029537397e0e046087f1eda5101d5ef6a853028d80f6c10c1ce6` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_ibp/lemmaS_P2_M0.6.json` | 4124 | `57974b9af4b24c69c08fa8274610e82dd835a71a5f633b4feba41b4a4d9ad43e` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_ibp/lemmaS_P2_M0.75.json` | 4133 | `3bf3106572e331bb756a9fb82320c075eaf44d4d15425aca024d2613e7d8866a` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_ibp/lemmaS_P3_M0.52.json` | 4767 | `26313994ab052ccb59923325f8767899493e0b8d86695a524dd9f19133d98e34` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_ibp/lemmaS_Q5o_M0.4.json` | 4135 | `fdbab7d04d714ad03d34f54acff61823c97ba0099cbad3c028eafaeeda00d7bc` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_ibp/lemmaS_Q5o_M0.5.json` | 4144 | `77fe373a3dc9bef8d1bb6d955dc2cbca67f806adf8c1a9e8b64d14e6b91f063f` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_ibp/lemmaS_Q5o_M0.6.json` | 4137 | `dbab7d393fb164cc3a57854f38de7ef6b409cf8aa1ae9fdf873ee977f90d1e5a` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_ibp/lemmaS_Q5p_M0.4.json` | 4134 | `7c153707248d7090ab9eab20aeb90112b49e451951f51c7b70a6febc11cb632a` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_ibp/lemmaS_Q5p_M0.5.json` | 4139 | `c5da7d76f031793d2e6afc15325d82016a6a263bb2621c83748de180af38b27d` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_ibp/lemmaS_Q5p_M0.6.json` | 4135 | `3c8a86d7a558bf7f422ac3fd40f56b726bbe91cbab5bf7bbe99b009f50d6ff70` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_ibp/lemmaS_R4_M0.5.json` | 5375 | `f84b79293059ac071adfcddd41453438a3d473f36886eb82b5d3631b9f22f4dc` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_ibp/lemmaS_R4_M0.55.json` | 5375 | `3bbe6045f6dcd17e91c7dc5ef31487c1a9b5a06c03c5045660ea9044fcbc22cc` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_ibp/lemmaS_R4_M0.6.json` | 5375 | `800faa0701fa5334f8dd8667258033efdb3a1b860c83a4f1e2893e5bdf890ccd` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_ibp/lemmaS_R5_M0.5.json` | 4556 | `36aa19148610d625e43c984a7f55ec831fa51fa75db5f01f4b21c9d8d50aa354` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_ibp/lemmaS_R5_M0.55.json` | 4552 | `6e30958e952fbb500920ac900b1a39777259e60385b703a0f365596a0a264205` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_ibp/lemmaS_R5_M0.6.json` | 4548 | `ee99475f669235f9e8f0f5d98b365e74ed89cf09879d00ea813433fc24694120` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_ibp/lemmaS_R8_M0.52.json` | 4762 | `8ebd62d82a40aad108c5469b9faedb4e24f91edfdd904bf1127e9cb4211cea58` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_ibp/lemmaS_S33_M0.5.json` | 4587 | `1071a9da1f09f62ca5a5677e3f26e34e2b2f9394444528597a44b20d6c7f0b74` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_ibp/lemmaS_S33_M0.55.json` | 4583 | `c921fc4acd01c2a8dae00eb1259b48a93c6fec961e70a8b09f45f3a225547adc` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_ibp/lemmaS_S33_M0.6.json` | 4576 | `12d50f04cfb5036c1a7b0b61587dd31fc29074cf637b00fdb638521c6356c770` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_ibp/lemmaS_S33b_M0.65.json` | 4569 | `6ffffb2a5f7400c3b8e0460a43e993801d529b70b4fe346c7ac3c32f340de98c` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_ibp/lemmaS_S33r300_M0.6.json` | 4566 | `04ac999e0e30f847a9843f47c1350bd494394214064d065617a71a7fc8fb8e62` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_ibp/lemmaS_S33r30_M0.5.json` | 4593 | `8b77a1eeb5f63bbdeb883d664a9add469319dd49e5be295d723cc932b5cf8383` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_ibp/lemmaS_S33r30_M0.55.json` | 4586 | `7fd5c0aab273a371ba503c23cba17a69687f6ac194c335f7d47ebef452ac3e94` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_ibp/lemmaS_S33r30_M0.6.json` | 4580 | `b53e0819f645855ef6236201a830ed9e19d578e30c10b3dec21fad4b6f4b053f` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_ibp/lemmaS_T0.5.json` | 4010 | `bad02deb488a813f31e4942f123d5b5ac249cb80c6e6510e4c44f68e4a50f133` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_ibp/lemmaS_T0.6.json` | 4000 | `39b81a799e4820b4a3480e3499c37b0b871aaa60a9e695bfebf544abef3094a8` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_ibp/lemmaS_T0.75.json` | 4007 | `bff48523000c260b252e6d810415306ccec83855fe5333a60a72f588e2fdeee9` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_ibp/lemmaS_T1.py` | 772 | `b26bcdc7314df1fc494ce534507563df1e1a823cd9c7f4dd9d737b45809dc72f` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_ibp/lemmaS_assemble.py` | 3535 | `04ce6d7664d8ebc8efbc76cf24f22b5d28bda34712d3a0c581cd2677c858d8bb` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_ibp/lemmaS_assembled.json` | 7448 | `d353bfd2792569f0e2bd31862a48bd62423013e3b821140d95e76108803ff04e` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_ibp/lemmaS_cflux.json` | 1235 | `f08e90fb25d0d8893a60109827334dc24f53fc37ed0cdf3f6e6d5407ca61bd65` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_ibp/lemmaS_cflux.py` | 6114 | `8c9250e3c8a5bc816d89b3a30768ee236508cd5bc203b8a103c1f53f52c98470` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_ibp/lemmaS_cflux_final.py` | 1309 | `fddbc40d6791625a301f9bb311a0c57cd0886dfb93f2eef0bb9ba326397e34b5` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_ibp/lemmaS_cflux_num.py` | 4569 | `0743352f2eba936435de0702e2422648ff6d7fbd0a9acdb529a1a43b06a95e24` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_ibp/lemmaS_consts.json` | 152622 | `5c28abe4eb6eb5e9e974a0dc95bbfe53bb24badc82ea7e6fde73293a1db6fb81` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_ibp/lemmaS_diag.json` | 4001 | `22f2daf77449927cef07a4992460c6915c0c9a762d42f85a54bf16c37acc1a59` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_ibp/lemmaS_interval.py` | 5692 | `9203d04b3f2599b21b7b36cc510c5e78794f646f7447c99f4244e007a56c11d8` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_ibp/lemmaS_mh0.json` | 3987 | `684a2b383ffda70f919697bcd56b40a9a925ad6b4eed14f7fe4b4d4ae2c1699d` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_ibp/lemmaS_mh10.json` | 3986 | `a1f4506a5e895d6f873d222586cc08fa7a9c27710083610b29bc6b7aa569db67` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_ibp/lemmaS_mh3.json` | 3991 | `320731384beb86990eb34344c419dbe361b593e5a068025cd0fd26ecf5de2d18` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_ibp/lemmaS_nested.json` | 1564 | `baab4e292b5c2eac6aa3a309eca39fdecf086f9cc0a7db935cc156f017a4cfac` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_ibp/lemmaS_nested.py` | 4147 | `fbcf2a8d6d1f7b444b283e59cdfc8b8e2ff1ef617a3dbd6a99e3c136c31eb767` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_ibp/lemmaS_prop1_bounds.py` | 5996 | `4265ea18a8eecf334c199e9271b59ee1ba89735baefb585d7e662f7cf3fa3c38` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_ibp/lemmaS_prop1_sym.py` | 3708 | `cbbc172a5872d969afad41cd2234d3e5fcc0b170a4a9bc078b159b0c82e5f52d` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_ibp/lemmaS_symbolic.py` | 1303 | `6933ba1a796b3cff3b1512298c18836df159db9097bc06f3e046cf2cf9b3e241` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_ibp/prev_assumed/F2_tight.log` | 951 | `79260c079b40cd458e371a1f609dfb35faa7abe934ac329cc489ba0db048e81c` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_ibp/prev_assumed/rigor_assembly_arb.log` | 12171 | `92e609fa16cfd809ad85c79be9841e2c1545bf44bbe378abe9b617ce5d104c0a` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_ibp/prev_assumed/rigor_assembly_arb.recheck.log` | 12171 | `92e609fa16cfd809ad85c79be9841e2c1545bf44bbe378abe9b617ce5d104c0a` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_ibp/prev_assumed/rigor_assembly_arb.recheck2.log` | 12171 | `92e609fa16cfd809ad85c79be9841e2c1545bf44bbe378abe9b617ce5d104c0a` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_ibp/prev_assumed/rigor_capture_F2.log` | 951 | `79260c079b40cd458e371a1f609dfb35faa7abe934ac329cc489ba0db048e81c` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_ibp/prev_assumed/rigor_capture_F2.recheck.log` | 951 | `79260c079b40cd458e371a1f609dfb35faa7abe934ac329cc489ba0db048e81c` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_ibp/prev_assumed/rigor_ramp_arb.log` | 173 | `ac517008262d69a13f2f20c6312477c0738fb45be5966646d99016219aed3b01` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_ibp/prev_assumed/rigor_weights_exact.log` | 2004 | `6f94ff75cda7f849b5e1c77dc4fcf92108ca448ede362a5c1b5873ababf85ee0` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_ibp/prop1_forms.pkl` | 2890 | `fbeff3a0d60cec287341eccf7051804e2e90238caf71845f432179d600a989e8` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_ibp/prop1_forms3.pkl` | 2599 | `1ba7a2c50a8e579b728ea7d230955618dc11c61a5d107e6395a3134f6dbf4a7c` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_ibp/ref/assembly_arb.json` | 11745 | `c8944912a74f38637612fbc6275d56238757464e3baf1a494cdddb49b60ec7e8` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_ibp/ref/lemmaS_F10p_M0.52.json` | 4832 | `686358eab56c78ccd00ccd8ea0dada226e1573ca9d6a8ae5b502c4cc041c5533` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_ibp/rigor/arbio.py` | 1076 | `dbf20c5e5b1a289bbf33f6c38af242cfa4d3be7ab2fbbf1948e2317eb04a5e0b` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_ibp/rigor/assembly_arb.py` | 24833 | `ef19f0bdfc8e47e9f56c31989e996f7ca4ec0da0652267611badbee142e2ecd1` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_ibp/rigor/capture_F2.py` | 2476 | `a5981657aced82d4f77696eac222f3aed04d166f1647f1226605795022aef615` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_ibp/rigor/ramp_arb.py` | 4473 | `477ffe1bfcbd195492f07c817e6ecb1b0197f277246d5c3568e1994a1aa4ea91` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_ibp/rigor/ramp_tables.json` | 352 | `764c4172903104e1520925bd06d1dd73479d1df89e16b9c60e0f9ecf4c0b168f` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_ibp/rigor/tables_exact.pkl` | 2839 | `861e3ac8b0bab6edec62221b299e3f80d543a8263fc83baa7fe6804c30e0ce17` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_ibp/rigor/weights_exact.py` | 19974 | `6c14948f5dbc16ad802eb5c0bf18bcb6ee1b7a2f01c26cfabe56381b258a7de0` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_ibp/rigor_assembly_arb.log` | 12433 | `6af221ba7d4234237e4b18d4216c128db916bb059b44d6d23b43dc944ea5ef88` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_ibp/rigor_assembly_arb_v33.log` | 12612 | `58c17d99f02e87be71272c5e552a4bd4eb8aede41a584782997540f29be75ba9` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_ibp/rigor_assembly_arb_v36.log` | 12727 | `c65eb93aea890d5b23575ff991ae8631906d8a0c2e215e0dc02e942878f3c412` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_ibp/rigor_capture_F2.log` | 951 | `7c0b9f6a5aec9ca8239b2021709b7d6666c341cc423a3372abc000787b5a0cac` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_ibp/rigor_capture_F2_v33.log` | 951 | `7c0b9f6a5aec9ca8239b2021709b7d6666c341cc423a3372abc000787b5a0cac` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_ibp/rigor_capture_F2_v36.log` | 951 | `5a96d157bbc035a71156f8d1f5fb440adf29dc5116bd6cf2e76ac7c3a83e59c9` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_ibp/rigor_weights_exact.log` | 2004 | `480b45947ea06426bf4f866512a62cf3145054bd964bc01607522277198451f4` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_ibp/rigor_weights_exact_v33.log` | 2004 | `480b45947ea06426bf4f866512a62cf3145054bd964bc01607522277198451f4` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_ibp/rigor_weights_exact_v36.log` | 2004 | `480b45947ea06426bf4f866512a62cf3145054bd964bc01607522277198451f4` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_ibp/run_F1.py` | 2541 | `59fefb073e9fa44fc2a4b8a0f509af6016d76fe9cb7c1a8baef406b24030a692` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_ibp/run_F10p_ibp.py` | 1044 | `d6a0a9073d2837e6e003dd24088d3aa9af6c1195d5552e1bab6d2daf9daddd08` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_ibp/v3.0_backup/lemmaS_Bbox_1e-05_a0.000125_dv0.15.json` | 3565 | `c68deea0230bbd6fd65987ff8f9cedba994513cb3d6078d46f9f42ebc3de0aea` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_ibp/v3.0_backup/lemmaS_Bbox_1e-05_a0.00015_dv0.15.json` | 3550 | `338655cf1c55782a29488303cc741869473edaa7142bb16571c22a234dfb2755` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_ibp/v3.0_backup/lemmaS_F10_M0.52.json` | 4814 | `526a3e91dc861c69f50fff611e03149e597d19e32f77b70904a1cd5aee2ad591` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_ibp/v3.0_backup/lemmaS_F10p_M0.52.json` | 4812 | `1b1d383caeb5ae5330c57ba533336916da337fa70dbd86be389f0d22959936d6` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_ibp/v3.0_backup/lemmaS_F4B_F10_M0.52.json` | 1167 | `dcf8c9a0fc64b001632221de652e45561eb69b9f34a7d3d7d055760d3c861971` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_ibp/v3.0_backup/lemmaS_F4B_F10p_M0.52.json` | 1167 | `c119cb4234e27a3039e29c4e89433bffeb88a164bdfa3bbba1cad2bc07eb9c1f` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_ibp/v3.0_backup/lemmaS_consts.json` | 152397 | `da78760e05006739e8f306830135c545cd05779702d853a180e74cda1e3b2cc1` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_ibp/v3.0_backup/rigor/assembly_arb.json` | 11710 | `8a67300c17ef2ffdda1b87294d67f736829b0351a901590d363f0065ecc6296c` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_ibp/v3.0_backup/rigor/assembly_arb_v3.0.py` | 23780 | `5a84ec24ec40385d8d1424cebc10c2afab27da0cfe14405ec3317dec58d4d8fc` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_ibp/v3.0_backup/rigor/cap_F10p.json` | 6599 | `6ee1ff9b413a539b9c0f8358ab2f2ea8f9b4f5e8eb1ea144bcf92a98e62e8b9d` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_ibp/v3.0_backup/rigor/capture_F2_v3.0.py` | 2410 | `e85875355c54dd04da6ba0723a25229c44b96298d38ba7ef04a6b86ce8c8cb6a` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_ibp/v3.0_backup/rigor/ramp_tables.json` | 352 | `764c4172903104e1520925bd06d1dd73479d1df89e16b9c60e0f9ecf4c0b168f` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_ibp/v3.0_backup/rigor/weights_exact.json` | 15438 | `5a3fb8acbbc30e7c88303179527bdfd033bb6d89c089d745785f100f68ecb3d3` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_v31H2/F2_ibp.log` | 947 | `909b9e8a59bd2397d1d257038d35beae3a8aa188f1872fe3be260e76a10424c3` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_v31H2/cflux_forms.pkl` | 4517 | `a862bc04ad6273a8d348029ca455af24367ab4e1d74e40f5f6071c43e2ab47a5` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_v31H2/lemmaS_B.json` | 666 | `b38f615bc26567d410a53b89aba8e93362cbc2b73b66bc283f240b840daeb2b9` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_v31H2/lemmaS_B.py` | 1477 | `38ff9f934df8cbe70b5158e80f9f822b438a5362ed887941fdf06580d360d7a9` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_v31H2/lemmaS_B3.json` | 235 | `a67d2ff97d151c76e37c91302011d84440193bdb954f86a5dc747fae0689b075` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_v31H2/lemmaS_B3.py` | 2147 | `71637489053d65219af0d5f2996c2f52025a4d7e33f4d64ed227b056053a3cab` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_v31H2/lemmaS_BR0_M0.6.json` | 4120 | `6b01aac83f96d434ca202d6d9da4983f265240c6938d21c4f049b4ea69e17287` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_v31H2/lemmaS_BR14.4_M0.6.json` | 4121 | `ec2552abd5a0e0eb273126f95612a02139eba2a63cb819dc16f9c3677a5f5685` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_v31H2/lemmaS_BR7.2_M0.6.json` | 4120 | `a3fa4d5c874fd006e95204dfc7144e2500aa6333db5c6a794241f25bce95ef44` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_v31H2/lemmaS_Bbox.py` | 4031 | `9c0fb29f05da5bf0d9102a90e3665c377c9e7cc7dd0dc6b41b9a92ffd0e54bb3` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_v31H2/lemmaS_Bbox_1e-05.json` | 3543 | `b2181cae4ea6d8943b2a03896e692048dc5a56b29d01d711f64946c4c99a0142` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_v31H2/lemmaS_Bbox_1e-05_a0.000125_dv0.15.json` | 3556 | `cc60b6d95c33240b688ce943cfaa7d5ba6da8b194e4ee7b2fd73b321028c1cd2` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_v31H2/lemmaS_Bbox_1e-05_a0.00015_dv0.15.json` | 3554 | `b68a1153af35561088edf8221d343b07bc8f9ca498919b43d6f42292f9ef397a` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_v31H2/lemmaS_Bbox_3e-05.json` | 3530 | `bcaa321b78a3064aa502079cf6a4081d58583db823be117b59e499325975f62b` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_v31H2/lemmaS_Bbox_3e-06.json` | 3546 | `83ddab479214606fd42350be17d4cb69758d11dfbfc47fa345174783c3f9dccf` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_v31H2/lemmaS_Bbox_3e-06_a0.000125_dv0.1.json` | 3571 | `53c8a4b73fa533f4ce49ca3e3ba4ab86762f69c33cf3cfac625bc88503a7b00b` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_v31H2/lemmaS_Bbox_3e-06_a0.000125_dv0.15.json` | 3559 | `1e3972ca0a1a5b21e388baa78dafb3a921b4221554c3cde2b6be5cba09df9b41` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_v31H2/lemmaS_Bbox_3e-06_dv05.json` | 3549 | `0259a710cb5dbd8db4cbc7c9d4f6a76c61983a5702fb434fcf409ebdd10cc0f6` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_v31H2/lemmaS_F1.json` | 21340 | `6c6a13762160891e6d1c1999d3b44959198521076cb32b9cfa8724b7e459d240` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_v31H2/lemmaS_F1.py` | 8232 | `08cc332aeddd4bb6424579e6d6482abb3545408f4002a978731d2fa1cb739937` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_v31H2/lemmaS_F10_M0.52.json` | 4828 | `c51da3461225199d54c0da2285529e7762c336d9d400bfbad924506cc4c69091` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_v31H2/lemmaS_F10a_M0.52.json` | 4780 | `81916706a1202a06ad078b5128d2b2e8819005301215353f461fcf044d1b5da9` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_v31H2/lemmaS_F10b_M0.4.json` | 4792 | `daa9d13765e172447bcadac08dad9f72ae11e98c361548e1d46c97e79b74428f` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_v31H2/lemmaS_F10b_M0.45.json` | 4797 | `06fac3d2b8188a952297e601e23c2d6add6b712acb9f48a144e6443d61635c9c` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_v31H2/lemmaS_F10b_M0.52.json` | 4796 | `367ee9ed19e78a0a3da78331eaa0c57ccbc8f631ae31ee6b65437993842d6770` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_v31H2/lemmaS_F10p_M0.52.json` | 4828 | `735b3c8e54aa34fb746a0bf59542aa9e5d727e8ac6599ebb25af27e07d6f61a5` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_v31H2/lemmaS_F1_scan.json` | 2807 | `2ff77d1f611ce5e5a73463b61dca08e20ae618cff1db8896b917d6318069ad8d` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_v31H2/lemmaS_F2.json` | 3960 | `30210a0dba8cb19c100131ec407c5556388d22a8dd596d9932e219d795a5140a` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_v31H2/lemmaS_F2.py` | 13823 | `cb736fdbf6a81dee286448971f385c521815e947a87d7eb267c5789e4a1b9082` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_v31H2/lemmaS_F4.json` | 1473 | `1d1ee5f105b1020aa99f17f3175d19f59274d666165d46591641639d5da5ac10` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_v31H2/lemmaS_F4.py` | 2128 | `dc1737a8ac4611c035539e0a5e9280c50c91e7a98af9290bf3d156fa9c158b7b` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_v31H2/lemmaS_F4B.py` | 5297 | `ec23c2ba0c34f47b8abef72e8f38f7f27e8bfa98df46bfc61b086cf5932f6094` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_v31H2/lemmaS_F4B_F10_M0.52.json` | 1170 | `362968b1083de0edc2533f771c43025d57f3ac90af3f5d5dd188a428da7b4b99` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_v31H2/lemmaS_F4B_F10p_M0.52.json` | 1169 | `7615261d270722285e6dafb55d080038b8c75ebb489eb03c8330eda49d0f875e` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_v31H2/lemmaS_F4B_F9_M0.52.json` | 1164 | `da530d0d66555979b1200b7f516605d926103e38ec1fb00ad9f2943ab78d6423` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_v31H2/lemmaS_F4B_L4b.json` | 1228 | `f546d1b5fc1e5dfbbba37901108d468103487c3ea8e52c0c4096fd4fab552a07` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_v31H2/lemmaS_F4B_N5b_M0.52.json` | 1157 | `d09c89ebec6bf6b91ed341b2534be7fcab73ac7cd858b18af9711c3b39d622c9` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_v31H2/lemmaS_F4B_P1_M0.6.json` | 1222 | `d535a0086e0b594ef52348428839bb25721967411c219d6e52b502bb9eeb0d10` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_v31H2/lemmaS_F4B_P2_M0.6.json` | 1227 | `99dfa5cfc5fb728d3cdb12a19665adbb91020f0c23a07007dfd4ea7710006511` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_v31H2/lemmaS_F4B_R5_M0.6.json` | 1170 | `7ff52e13bb6850c3ebf2780d20613de0e0f1060924de357f846e8f9051f75908` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_v31H2/lemmaS_F4B_R8_M0.52.json` | 1165 | `a933f761039a11cc893183f819b42dec091a6c47a725dbe2b3baf04fe57de3fb` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_v31H2/lemmaS_F5_geom.json` | 733 | `61fbffa478db8f0b9e0d8cc67a2d0501e827b03f3d9ecabe406cf0026b419427` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_v31H2/lemmaS_F5_geom.py` | 898 | `c440953010a8d4b8a66aaf65bd3f17f5a47d01cf7d257dded0549543a7ae6e02` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_v31H2/lemmaS_F9_M0.52.json` | 4757 | `a60f0d3edd9f870fa2a41dd37f145acb543dcd6c35c340c3294bfd307edf8bc4` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_v31H2/lemmaS_GN.py` | 2241 | `f459d78489fc7358e5d8b60c8481460955a8816ea47bb412290f0460bfb6d58e` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_v31H2/lemmaS_H6.py` | 1111 | `71b6ac231cf07e43b9faa6da5cf39865c6b09cc2fef3c629dd4a9b6738e6a415` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_v31H2/lemmaS_L4.json` | 3988 | `f0d39a54b7c264f81df91b8444ddbfb0dae9bb210111ff483cafa8aa9125f266` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_v31H2/lemmaS_L4b.json` | 3998 | `35d9325783d9c4cdaa84aaa285caade2c91ccec9ecbd0b74f7cc052f65fa6c39` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_v31H2/lemmaS_L4terms.py` | 8708 | `65ecceee1b3f1274749cc0c4a5cab6592a913ab4db8e16ff57cf553e4be782a5` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_v31H2/lemmaS_N1_M0.4.json` | 4590 | `a9d48e139c9b07d7075f6bf5c924effd67ae043703530281eb2670f98b8cad84` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_v31H2/lemmaS_N1_M0.5.json` | 4602 | `23a9fd3f80b5fb2acf296af21d2c92522146554e08c9f3adb6f59c76e7751b02` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_v31H2/lemmaS_N1_M0.55.json` | 4601 | `198a7cfe906947d3a43ace94c2cce7c1a5d6b1c7513fbaa1951ea591f1627621` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_v31H2/lemmaS_N1_M0.6.json` | 4591 | `82ea3ccc91ef2eb84d55fb09904659720eb5d68836fb8b1ef31064f11bffa7ff` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_v31H2/lemmaS_N1_M0.65.json` | 4597 | `d73dfb510d000aed9a83b4e79ec8b15a1bd8276947d03e9a7aa28fc387043da7` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_v31H2/lemmaS_N1b_M0.575.json` | 4570 | `9aa6d22c2776defda7176841c777d52aa2096e3a63af7ed799316b244aa09318` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_v31H2/lemmaS_N5_M0.3.json` | 4599 | `f90ddde1005847ccd7658d4dfd7b6b15d67b09ea616ff274bfe9a1bfee3c6ad2` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_v31H2/lemmaS_N5_M0.4.json` | 4595 | `4b757b5ead94125cafe18ba6bbb67bef3603599e5f1b8bd3a66db81b54957c08` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_v31H2/lemmaS_N5_M0.5.json` | 4607 | `32e9a362b6f561a58f9881cf5db4ee8c1a67251d24a58852508558a3619de879` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_v31H2/lemmaS_N5_M0.6.json` | 4598 | `ef9e401d5f47fa70b1a3c6149d571e4497c10cd8532058ef05a82ddf3c689f80` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_v31H2/lemmaS_N5_M0.65.json` | 4604 | `f6bd656b12e33eee527f0560f085b8948b1cbd1b1406094d3afaa33e8b8eb7da` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_v31H2/lemmaS_N5b_M0.52.json` | 4570 | `745bbdaeb6eb388414ee9059d93e7ab26ec895334c519530dae111617289036a` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_v31H2/lemmaS_P10000_M0.52.json` | 4594 | `7a730960b771ceba751763b4a45d8a12a050195bc8daaa0f05134e67c6ee0c89` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_v31H2/lemmaS_P1000_M0.52.json` | 4595 | `afbd50f30d5e78674a304905e6e61f3699ebef8eec9056b6d8621ae8d142909e` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_v31H2/lemmaS_P1_M0.5.json` | 4140 | `ff476c422aed01e019e72794b9022f0ba730ad0327379381b70983c2a01dbbf1` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_v31H2/lemmaS_P1_M0.6.json` | 4134 | `f764274adaae6995eca3d94cda27442c85e8eb4f8eb96617d0614ac3613ae0ad` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_v31H2/lemmaS_P1_M0.75.json` | 4143 | `29780ffdc2396751bbb8aa0aafb932bbd12353b23401f29622fcbc63a0a94ed4` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_v31H2/lemmaS_P2_M0.5.json` | 4135 | `ad10632395e5840eccc02c3eabd0a7f106da755eba0b3fd77261ece3e73557c1` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_v31H2/lemmaS_P2_M0.52.json` | 4772 | `62e1247f9cdf029537397e0e046087f1eda5101d5ef6a853028d80f6c10c1ce6` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_v31H2/lemmaS_P2_M0.6.json` | 4124 | `57974b9af4b24c69c08fa8274610e82dd835a71a5f633b4feba41b4a4d9ad43e` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_v31H2/lemmaS_P2_M0.75.json` | 4133 | `3bf3106572e331bb756a9fb82320c075eaf44d4d15425aca024d2613e7d8866a` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_v31H2/lemmaS_P3_M0.52.json` | 4767 | `26313994ab052ccb59923325f8767899493e0b8d86695a524dd9f19133d98e34` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_v31H2/lemmaS_Q5o_M0.4.json` | 4135 | `fdbab7d04d714ad03d34f54acff61823c97ba0099cbad3c028eafaeeda00d7bc` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_v31H2/lemmaS_Q5o_M0.5.json` | 4144 | `77fe373a3dc9bef8d1bb6d955dc2cbca67f806adf8c1a9e8b64d14e6b91f063f` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_v31H2/lemmaS_Q5o_M0.6.json` | 4137 | `dbab7d393fb164cc3a57854f38de7ef6b409cf8aa1ae9fdf873ee977f90d1e5a` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_v31H2/lemmaS_Q5p_M0.4.json` | 4134 | `7c153707248d7090ab9eab20aeb90112b49e451951f51c7b70a6febc11cb632a` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_v31H2/lemmaS_Q5p_M0.5.json` | 4139 | `c5da7d76f031793d2e6afc15325d82016a6a263bb2621c83748de180af38b27d` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_v31H2/lemmaS_Q5p_M0.6.json` | 4135 | `3c8a86d7a558bf7f422ac3fd40f56b726bbe91cbab5bf7bbe99b009f50d6ff70` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_v31H2/lemmaS_R4_M0.5.json` | 5375 | `f84b79293059ac071adfcddd41453438a3d473f36886eb82b5d3631b9f22f4dc` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_v31H2/lemmaS_R4_M0.55.json` | 5375 | `3bbe6045f6dcd17e91c7dc5ef31487c1a9b5a06c03c5045660ea9044fcbc22cc` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_v31H2/lemmaS_R4_M0.6.json` | 5375 | `800faa0701fa5334f8dd8667258033efdb3a1b860c83a4f1e2893e5bdf890ccd` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_v31H2/lemmaS_R5_M0.5.json` | 4556 | `36aa19148610d625e43c984a7f55ec831fa51fa75db5f01f4b21c9d8d50aa354` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_v31H2/lemmaS_R5_M0.55.json` | 4552 | `6e30958e952fbb500920ac900b1a39777259e60385b703a0f365596a0a264205` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_v31H2/lemmaS_R5_M0.6.json` | 4548 | `ee99475f669235f9e8f0f5d98b365e74ed89cf09879d00ea813433fc24694120` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_v31H2/lemmaS_R8_M0.52.json` | 4762 | `8ebd62d82a40aad108c5469b9faedb4e24f91edfdd904bf1127e9cb4211cea58` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_v31H2/lemmaS_S33_M0.5.json` | 4587 | `1071a9da1f09f62ca5a5677e3f26e34e2b2f9394444528597a44b20d6c7f0b74` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_v31H2/lemmaS_S33_M0.55.json` | 4583 | `c921fc4acd01c2a8dae00eb1259b48a93c6fec961e70a8b09f45f3a225547adc` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_v31H2/lemmaS_S33_M0.6.json` | 4576 | `12d50f04cfb5036c1a7b0b61587dd31fc29074cf637b00fdb638521c6356c770` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_v31H2/lemmaS_S33b_M0.65.json` | 4569 | `6ffffb2a5f7400c3b8e0460a43e993801d529b70b4fe346c7ac3c32f340de98c` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_v31H2/lemmaS_S33r300_M0.6.json` | 4566 | `04ac999e0e30f847a9843f47c1350bd494394214064d065617a71a7fc8fb8e62` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_v31H2/lemmaS_S33r30_M0.5.json` | 4593 | `8b77a1eeb5f63bbdeb883d664a9add469319dd49e5be295d723cc932b5cf8383` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_v31H2/lemmaS_S33r30_M0.55.json` | 4586 | `7fd5c0aab273a371ba503c23cba17a69687f6ac194c335f7d47ebef452ac3e94` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_v31H2/lemmaS_S33r30_M0.6.json` | 4580 | `b53e0819f645855ef6236201a830ed9e19d578e30c10b3dec21fad4b6f4b053f` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_v31H2/lemmaS_T0.5.json` | 4010 | `bad02deb488a813f31e4942f123d5b5ac249cb80c6e6510e4c44f68e4a50f133` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_v31H2/lemmaS_T0.6.json` | 4000 | `39b81a799e4820b4a3480e3499c37b0b871aaa60a9e695bfebf544abef3094a8` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_v31H2/lemmaS_T0.75.json` | 4007 | `bff48523000c260b252e6d810415306ccec83855fe5333a60a72f588e2fdeee9` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_v31H2/lemmaS_T1.py` | 772 | `b26bcdc7314df1fc494ce534507563df1e1a823cd9c7f4dd9d737b45809dc72f` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_v31H2/lemmaS_assemble.py` | 3535 | `04ce6d7664d8ebc8efbc76cf24f22b5d28bda34712d3a0c581cd2677c858d8bb` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_v31H2/lemmaS_assembled.json` | 7448 | `d353bfd2792569f0e2bd31862a48bd62423013e3b821140d95e76108803ff04e` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_v31H2/lemmaS_cflux.json` | 1235 | `f08e90fb25d0d8893a60109827334dc24f53fc37ed0cdf3f6e6d5407ca61bd65` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_v31H2/lemmaS_cflux.py` | 6114 | `8c9250e3c8a5bc816d89b3a30768ee236508cd5bc203b8a103c1f53f52c98470` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_v31H2/lemmaS_cflux_final.py` | 1309 | `fddbc40d6791625a301f9bb311a0c57cd0886dfb93f2eef0bb9ba326397e34b5` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_v31H2/lemmaS_cflux_num.py` | 4569 | `0743352f2eba936435de0702e2422648ff6d7fbd0a9acdb529a1a43b06a95e24` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_v31H2/lemmaS_consts.json` | 152622 | `5c28abe4eb6eb5e9e974a0dc95bbfe53bb24badc82ea7e6fde73293a1db6fb81` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_v31H2/lemmaS_diag.json` | 4001 | `22f2daf77449927cef07a4992460c6915c0c9a762d42f85a54bf16c37acc1a59` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_v31H2/lemmaS_interval.py` | 5692 | `9203d04b3f2599b21b7b36cc510c5e78794f646f7447c99f4244e007a56c11d8` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_v31H2/lemmaS_mh0.json` | 3987 | `684a2b383ffda70f919697bcd56b40a9a925ad6b4eed14f7fe4b4d4ae2c1699d` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_v31H2/lemmaS_mh10.json` | 3986 | `a1f4506a5e895d6f873d222586cc08fa7a9c27710083610b29bc6b7aa569db67` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_v31H2/lemmaS_mh3.json` | 3991 | `320731384beb86990eb34344c419dbe361b593e5a068025cd0fd26ecf5de2d18` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_v31H2/lemmaS_nested.json` | 1564 | `baab4e292b5c2eac6aa3a309eca39fdecf086f9cc0a7db935cc156f017a4cfac` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_v31H2/lemmaS_nested.py` | 4147 | `fbcf2a8d6d1f7b444b283e59cdfc8b8e2ff1ef617a3dbd6a99e3c136c31eb767` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_v31H2/lemmaS_prop1_bounds.py` | 5996 | `4265ea18a8eecf334c199e9271b59ee1ba89735baefb585d7e662f7cf3fa3c38` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_v31H2/lemmaS_prop1_sym.py` | 3708 | `cbbc172a5872d969afad41cd2234d3e5fcc0b170a4a9bc078b159b0c82e5f52d` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_v31H2/lemmaS_symbolic.py` | 1303 | `6933ba1a796b3cff3b1512298c18836df159db9097bc06f3e046cf2cf9b3e241` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_v31H2/prop1_forms.pkl` | 2890 | `fbeff3a0d60cec287341eccf7051804e2e90238caf71845f432179d600a989e8` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_v31H2/prop1_forms3.pkl` | 2599 | `1ba7a2c50a8e579b728ea7d230955618dc11c61a5d107e6395a3134f6dbf4a7c` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_v31H2/ref/assembly_arb.json` | 11745 | `c8944912a74f38637612fbc6275d56238757464e3baf1a494cdddb49b60ec7e8` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_v31H2/ref/lemmaS_F10p_M0.52.json` | 4832 | `686358eab56c78ccd00ccd8ea0dada226e1573ca9d6a8ae5b502c4cc041c5533` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_v31H2/rigor/arbio.py` | 1076 | `dbf20c5e5b1a289bbf33f6c38af242cfa4d3be7ab2fbbf1948e2317eb04a5e0b` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_v31H2/rigor/assembly_arb.json` | 11744 | `6f60c4f54df92412ef07ed7507d35422f6cd65fc555863812c780f0887928483` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_v31H2/rigor/assembly_arb.py` | 24111 | `91fdcbf07f7bc4305df5c03296819900e4c2460088d4e8658ed3ce0dba60526e` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_v31H2/rigor/cap_F10p.json` | 6615 | `d5f0baa87ffdb16b0e21bd395a4f4f91ffbd7d7fee7e8cff909ad8dfcf77180d` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_v31H2/rigor/capture_F2.py` | 2435 | `99e4e0382500f92860238f489889d37c67d863a4ad403234e75dd6852985fd99` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_v31H2/rigor/lemmaS_rigcap_M0.52.json` | 4830 | `a7487326e07047323e27bed4a3c540d29c17783a7577373c247ebfaea1030785` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_v31H2/rigor/ramp_arb.py` | 4473 | `477ffe1bfcbd195492f07c817e6ecb1b0197f277246d5c3568e1994a1aa4ea91` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_v31H2/rigor/ramp_tables.json` | 352 | `764c4172903104e1520925bd06d1dd73479d1df89e16b9c60e0f9ecf4c0b168f` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_v31H2/rigor/tables_exact.pkl` | 2839 | `861e3ac8b0bab6edec62221b299e3f80d543a8263fc83baa7fe6804c30e0ce17` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_v31H2/rigor/weights_exact.json` | 15438 | `0d17e927088056ff37fbb6785fc1d60388b03fa9e2947bf2c9c43c81adca7281` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_v31H2/rigor/weights_exact.py` | 19974 | `6c14948f5dbc16ad802eb5c0bf18bcb6ee1b7a2f01c26cfabe56381b258a7de0` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_v31H2/rigor_assembly_arb.log` | 12695 | `8065263299dd0ad8b97d948e398d4ee302d6a75efbb7db6497a3d5ca37e28e20` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_v31H2/rigor_capture_F2.log` | 947 | `909b9e8a59bd2397d1d257038d35beae3a8aa188f1872fe3be260e76a10424c3` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_v31H2/rigor_weights_exact.log` | 2004 | `480b45947ea06426bf4f866512a62cf3145054bd964bc01607522277198451f4` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_v31H2/run_F1.py` | 2541 | `59fefb073e9fa44fc2a4b8a0f509af6016d76fe9cb7c1a8baef406b24030a692` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_v31H2/run_F10p_ibp.py` | 900 | `c69ba17561bdb8fd2e28b6c030f59c2a8cd2f8a7f0ce4f56b7616c413549e197` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_v31H2/v3.0_backup/lemmaS_Bbox_1e-05_a0.000125_dv0.15.json` | 3565 | `c68deea0230bbd6fd65987ff8f9cedba994513cb3d6078d46f9f42ebc3de0aea` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_v31H2/v3.0_backup/lemmaS_Bbox_1e-05_a0.00015_dv0.15.json` | 3550 | `338655cf1c55782a29488303cc741869473edaa7142bb16571c22a234dfb2755` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_v31H2/v3.0_backup/lemmaS_F10_M0.52.json` | 4814 | `526a3e91dc861c69f50fff611e03149e597d19e32f77b70904a1cd5aee2ad591` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_v31H2/v3.0_backup/lemmaS_F10p_M0.52.json` | 4812 | `1b1d383caeb5ae5330c57ba533336916da337fa70dbd86be389f0d22959936d6` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_v31H2/v3.0_backup/lemmaS_F4B_F10_M0.52.json` | 1167 | `dcf8c9a0fc64b001632221de652e45561eb69b9f34a7d3d7d055760d3c861971` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_v31H2/v3.0_backup/lemmaS_F4B_F10p_M0.52.json` | 1167 | `c119cb4234e27a3039e29c4e89433bffeb88a164bdfa3bbba1cad2bc07eb9c1f` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_v31H2/v3.0_backup/lemmaS_consts.json` | 152397 | `da78760e05006739e8f306830135c545cd05779702d853a180e74cda1e3b2cc1` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_v31H2/v3.0_backup/rigor/assembly_arb.json` | 11710 | `8a67300c17ef2ffdda1b87294d67f736829b0351a901590d363f0065ecc6296c` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_v31H2/v3.0_backup/rigor/assembly_arb_v3.0.py` | 23780 | `5a84ec24ec40385d8d1424cebc10c2afab27da0cfe14405ec3317dec58d4d8fc` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_v31H2/v3.0_backup/rigor/cap_F10p.json` | 6599 | `6ee1ff9b413a539b9c0f8358ab2f2ea8f9b4f5e8eb1ea144bcf92a98e62e8b9d` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_v31H2/v3.0_backup/rigor/capture_F2_v3.0.py` | 2410 | `e85875355c54dd04da6ba0723a25229c44b96298d38ba7ef04a6b86ce8c8cb6a` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_v31H2/v3.0_backup/rigor/ramp_tables.json` | 352 | `764c4172903104e1520925bd06d1dd73479d1df89e16b9c60e0f9ecf4c0b168f` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/base_v31H2/v3.0_backup/rigor/weights_exact.json` | 15438 | `5a3fb8acbbc30e7c88303179527bdfd033bb6d89c089d745785f100f68ecb3d3` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/rigor/a3_signs.py` | 2673 | `45b80db67f24e5395e2d8c0485da7c9d586f0a3c5cf6e7ddb33f9b8c7e22f85e` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/rigor/arlib.py` | 4233 | `5838a787226c94df363652ca6e303fd74fff507dc5ccc6f86e371520904990b0` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/rigor/assembly_rep.py` | 35352 | `6f8d138113ff146d855c1d0629b908c4b177056797a470689850e7c7b506e324` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/rigor/export_for_U.py` | 2624 | `b53f1162f801049c28cbc9501f6a85605d9dedb030f9b3a71416ec0f826ee4cb` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/rigor/export_t2.py` | 15592 | `601cf98cf1e764af4d24c8d54d9d54957dd7eb2dca7826c114bce9f846ff21c8` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/rigor/hjw_rep.py` | 5215 | `c91056570a2d8f9db139e2fa00d685f5f6ecae400f6bd40c996527bad0593d0c` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/rigor/ic_rep.py` | 2301 | `c70cfb3fbeb8e621d102518acc5e94621f79d2e719d9323304a5567ca3fc3255` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/rigor/outer_L2.py` | 15980 | `1e6a39995ae0ac8effdca9724b9e9f42677653fb7e96b4ae3fb2f4140af22317` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/rigor/outer_rep.py` | 7929 | `925010953d7d009435c3a3d1d89d1303282ed1924ff4c4fa3244d8a899cf95c5` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/rigor/ramp_rep.py` | 16391 | `461caab008f922152b370bec971b0aaf82634a0c43806c06024cdcad067343c2` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/rigor/weights_rep.py` | 12692 | `38fe9a8e91a95c6d92eabe4c6857a42dd78861ed919cff36baed8d669928fd69` |
| `writeup/lemmaU/bulk/fix13_v37chain/t2_A/rigor/weights_t2.py` | 24836 | `ed87327e5f0974b01f768873fb2b206cfcae42956267bea0abee1d7b092b90e2` |
| `writeup/lemmaU/bulk/fix13_v37chain/t4_cert/assemble_U.py` | 22715 | `d14baeea1f0f8bdb703ca58aa5b40cd4c4f288beacf182443d05ff54a488cfce` |
| `writeup/lemmaU/bulk/fix13_v37chain/t4_cert/make_ucert_for_zone3.py` | 996 | `1e0a25c3dc05f0c254e62eb895b9a5a8b7e9ad64ff3753b438cd5aa9f3e020cd` |
| `writeup/lemmaU/bulk/fix13_v37chain/t4_cert/rampphi.py` | 2608 | `c4a8291e7d417620ba6808a3331cc621b8e2ade3373c37bcca2963b268eea3a3` |
| `writeup/lemmaU/bulk/fix13_v37chain/t4_cert/t4lib.py` | 3436 | `546f4433c5038dd5a826ee545f484daa9e2201c8c2851ecfd3eda7c21d673bca` |
| `writeup/lemmaU/bulk/fix13_v37chain/t4_cert/u1_ghat.py` | 11244 | `daca142d8fd2d12626fcb757955f39d7fddc1b189ffa2fa15aeafea705c24f79` |
| `writeup/lemmaU/bulk/fix13_v37chain/t4_cert/u2_K.py` | 13495 | `059a27348804929c6f6dec1972a25f51db1e63f8baadbe36571f356ff8f9492c` |
| `writeup/lemmaU/bulk/fix13_v37chain/t4_cert/u34_sweeps.py` | 11495 | `43501a536cae4c800c3fae699d2094932ce3dc7ea1ac84959131ace7f0995712` |
| `writeup/lemmaU/bulk/fix13_v37chain/t4_cert/u5_enum.py` | 7348 | `8c8508201e7e8c2d07ce5612e1fe6ca5047deced82e8187b97926812b4e12432` |
| `writeup/lemmaU/bulk/fix13_v37chain/t5_zone3/out/ramp_U.json` | 6631 | `55317ce5a31e8e3c3b3c4b86e85ef9569bb6d5e99190b68b0f21e4cba1a0c0de` |
| `writeup/lemmaU/bulk/fix13_v37chain/t5_zone3/ramp_U.py` | 14866 | `0d4566ae3f0dd490fd1964046851125b64bf3c5fd29a6a5e31fb78054bc7a24c` |
| `writeup/lemmaU/bulk/fix13_v37chain/t5_zone3/z1_box.py` | 10607 | `04409b96688c5b592475f5396dd57c266fd3e33740bd92a0b0c94e7e1880111b` |
| `writeup/lemmaU/bulk/fix13_v37chain/t5_zone3/z3lib.py` | 6280 | `193b843d319aef12ddbfcb3ea612b2d2c441c59461a400da774368f5853e5971` |
| `writeup/lemmaU/bulk/fix13_v37chain/t5_zone3/zone3.py` | 23563 | `b78ce0fbd32e4920b8af999131497d2089415e940e1944126ef20392cb736654` |
| `writeup/lemmaU/bulk/fix13_v37chain/v36cfg.py` | 7424 | `1d5e2cd6bfdf5738e054a8f13ec27c91cd2ccb49a40f73253d3855459a826002` |
| `writeup/lemmaU/bulk/fix3/flux/flux_cert.json` | 4796 | `1e4a262ef1b9dfcde084c76381af3596d4fb98d094476cc6fb99e773e9f85a91` |
| `writeup/lemmaU/bulk/fix3/flux/flux_forms.pkl` | 24001 | `56202f4023447377fc701a5a833d2245c61ccf3b8ba9271ba61e57140fef27ba` |
| `writeup/lemmaU/bulk/fix7_flux2/FLUX2_DESIGN.md` | 44737 | `4405dce058cd7409d02f48a0940b348cb70a227204fc422251497238fa253db3` |
| `writeup/lemmaU/bulk/fix7_flux2/flux2_lemma.tex` | 19451 | `06e27c02e026cc85d1599d9ce5a4abd3e1244663e29958bacf624e927aa62ec8` |
| `writeup/lemmaU/bulk/t1_statements/r2_probe/ramp_out_253_379.5_box0.15_0.52.json` | 430 | `b5dd32b18f9f07dec1c8da3541841dfd8e0c57285c57d4f7ae05cc7b0324448d` |
| `writeup/lemmaU/bulk/t5_zone3/out/zone3_interface.json` | 18831 | `45e8d4dcaff2bd75fcea725f80e98625af70a40709fc2cf5d3df78b5becb3ea3` |
| `writeup/lemmaU/phase1/ghat/arbio_local.py` | 1124 | `7e7c5ac37854a38149e327ad5dbee3c73aa6f236dbe47595a6ca1e2a71489ced` |
| `writeup/lemmaU/phase1/ghat/ghat_arb.py` | 19025 | `4222d46a93707e2aebf07945e44b935405be15bf8643b415fabb9970fabf3461` |
| `writeup/lemmaU/phase1/ghat/weight_consts_arb.py` | 9841 | `0083df27fcddf81eb7f609954db84f684de1f67bfefbd60a6e86e74345b1ffb3` |
| `writeup/lemmaU/phase1/ramp/ulib.py` | 16092 | `672ff0eba7884b575f411e080c0602816370d548d6410ca15bffea8446d7f1ac` |
| `writeup/lemmaU/phase2a/U_ramp/ghat_design_arb.json` | 10161 | `0e1fbe40721dbc99d400c23c58fee935067c7803ff8cbd5bf09fd2120c8512c2` |
