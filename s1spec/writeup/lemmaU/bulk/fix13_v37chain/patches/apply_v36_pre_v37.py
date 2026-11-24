"""fix6_v36: the v36 edits, applied on top of the fix4_r3 patches (paths_*, r3_*).  Each replacement must match exactly once.
Run once: python3 -I -B patches/apply_v36.py  (from anywhere).  Then `diff -u X_pre_v36.py X.py` gives the full change vs v35."""
import os, sys
N = '/home/claude/s1spec/writeup/lemmaU/bulk/fix6_v36'
CFGIMP = "sys.path.insert(0, '%s'); import v36cfg as CFG   # fix6_v36 (v36): verifier source + switches V36_VER/V36_R3/V36_R2" % N
EDITS = {
 't2_A/rigor/arlib.py': [
  ("ROOT = '/home/claude/s1spec/writeup/lemmaU/bulk/fix4_r3/t2_A'   # fix4_r3 (R3) copy of fix2/n12/t2_A_v35 (paths only); was the N1 copy (v35) of fix1/b0/t2_A_v34",
   "ROOT = '%s/t2_A'   # fix6_v36 (v36) copy of fix2/n12/t2_A_v35 (+ fix4_r3 R3 patches, gated by V36_R3); was the N1 copy (v35) of fix1/b0/t2_A_v34\n" % N + CFGIMP),
  ("env={k: os.environ.get(k) for k in ('PYTHONHASHSEED', 'PYTHONDONTWRITEBYTECODE', 'OMP_NUM_THREADS')},",
   "env={k: os.environ.get(k) for k in ('PYTHONHASHSEED', 'PYTHONDONTWRITEBYTECODE', 'OMP_NUM_THREADS') + CFG.ENVKEYS},"),
 ],
 't2_A/run_t2.sh': [
  ("set -e; cd /home/claude/s1spec/writeup/lemmaU/bulk/fix4_r3/t2_A",
   "set -e; cd %s/t2_A\n# v36: the chart input follows V36_VER (fin12f4: pinned step-4 v33 n1_chart.json; fin13v4b: n1_chart.py copy rerun on the fin13v4b sups)\n"
   "CHART=$(python3 -I -B -c \"import sys; sys.path.insert(0, '%s'); import v36cfg; print(v36cfg.N1CHART)\")" % (N, N)),
  ("--chart /home/claude/s1spec/step4/p2/out_v33/n1_chart.json\"   # fix B0 + N1", "--chart $CHART\"   # fix B0 + N1 (v36: chart per V36_VER)"),
 ],
 't2_A/rigor/hjw_rep.py': [
  ("RT_DEF = '/home/claude/s1spec/step3nl_v33/verifier/rate_table_fin12f4_a0.00015.json'", "RT_DEF = L.CFG.RT                                   # v36: verifier per V36_VER (was fin12f4)"),
 ],
 't2_A/base_ibp/run_F10p_ibp.py': [
  ("import lemmaS_F2 as F2\n", "import lemmaS_F2 as F2\n" + CFGIMP + "\n"),
  ("""sys.argv = ['lemmaS_F2.py', '--J', '6.0273e-12', '8.4673e-12', '1.8708e-12', '--handoff', '1e-8', '1e-6', '1e-3', '--sqLam', '5.3075e-8', '--alpha1', '1.5e-4',
            '--core', '0.11', '0.52', '--ramp', '0.05', '0.9', '--mhr', '0', '--slopes', '3.4548e-18', '4.8534e-18', '1.0724e-18', '--rho', '2.693e-3', '--M1s', '0.52',
            '--J4', '6.2362e-8', '--sJ4', '3.5746e-14', '--rampl4', '--rs33', 'ibp2', '--s33', 'ibp', '--l4', '--abar', '2.19103e-7', '--beta', '3.33811e-9', '--tag', 'F10p']""",
   "sys.argv = CFG.f2_argv('F10p')   # v36: verifier inputs read from the raw JSON, rounded outward (v36cfg.verifier_inputs); fin12f4 = the v35 typed argv (asserted)"),
 ],
 't2_A/base_ibp/rigor/capture_F2.py': [
  ("import lemmaS_F1 as F, lemmaS_F2 as F2\n", "import lemmaS_F1 as F, lemmaS_F2 as F2\n" + CFGIMP + "\n"),
  ("""sys.argv = ['lemmaS_F2.py', '--J', '6.0273e-12', '8.4673e-12', '1.8708e-12', '--handoff', '1e-8', '1e-6', '1e-3', '--sqLam', '5.3075e-8', '--alpha1', '1.5e-4',
            '--core', '0.11', '0.52', '--ramp', '0.05', '0.9', '--mhr', '0', '--slopes', '3.4548e-18', '4.8534e-18', '1.0724e-18', '--rho', '2.693e-3', '--M1s', '0.52',
            '--J4', '6.2362e-8', '--sJ4', '3.5746e-14', '--rampl4', '--rs33', 'ibp2', '--s33', 'ibp', '--l4', '--abar', '2.19103e-7', '--beta', '3.33811e-9', '--tag', 'rigcap']""",
   "sys.argv = CFG.f2_argv('rigcap')   # v36: as run_F10p_ibp.py"),
 ],
 't2_A/base_ibp/rigor/assembly_arb.py': [
  ("sys.path.insert(0, '/home/claude/s1spec/writeup/lemmaU/bulk/fix4_r3/r3'); import r3_forcing as R3   # fix4_r3 (R3, v36pre): honest forcing in window_loop, f4b (cF, (B) exit)",
   "sys.path.insert(0, '%s/r3'); import r3_forcing as R3   # fix4_r3 (R3, v36pre): honest forcing in window_loop, f4b (cF, (B) exit); v36: gated by V36_R3\n" % N + CFGIMP),
  ("    r0 = R3.r0_honest(BETA, AL, ABAR, ABAR); assert GAM == 1; cmin = D('0.25') - ABAR - ABAR/50",
   "    r0 = R3.r0_honest(BETA, AL, ABAR, ABAR) if CFG.R3ON else D('1.4')*(BETA + D('1.02')*AL*ABAR); assert GAM == 1; cmin = D('0.25') - ABAR - ABAR/50"),
  ("            EB.append(Es[j]*Cexp.exp()*(8*ABAR/yE).exp() + 4*R3.F_honest(j, r0, IB, R3.feps2(j, yE, yb, yE, AL, n=100))**2/gp[j]**2 + 2*coup/gp[j])   # R3",
   "            FB2 = R3.F_honest(j, r0, IB, R3.feps2(j, yE, yb, yE, AL, n=100))**2 if CFG.R3ON else (r0*P**j)**2*IB\n"
   "            EB.append(Es[j]*Cexp.exp()*(8*ABAR/yE).exp() + 4*FB2/gp[j]**2 + 2*coup/gp[j])   # R3 (gated)"),
  ("    cF = R3.F_honest(1, R3.r0_honest(BETA, AL, ABAR, ABAR), Iint, R3.feps2(1, 14, 350, 28, AL))/SQ   # R3 (was 1.4(...)/SQ p sqrt(Iint))",
   "    cF = R3.F_honest(1, R3.r0_honest(BETA, AL, ABAR, ABAR), Iint, R3.feps2(1, 14, 350, 28, AL))/SQ if CFG.R3ON else D('1.4')*(BETA + D('1.02')*AL*ABAR)/SQ*P*sqrt(Iint)   # R3 (gated)"),
  ("            EB.append(E[j]*Cexp.exp()*(4*2*ABAR/yE).exp() + 4*R3.F_honest(j, r0, IB, R3.feps2(j, yE, Rout, yE, AL, n=100))**2/g**2 + 2*coup/g); gpB.append(g)   # R3",
   "            FB2 = R3.F_honest(j, r0, IB, R3.feps2(j, yE, Rout, yE, AL, n=100))**2 if CFG.R3ON else (r0*P**j)**2*IB\n"
   "            EB.append(E[j]*Cexp.exp()*(4*2*ABAR/yE).exp() + 4*FB2/g**2 + 2*coup/g); gpB.append(g)   # R3 (gated)"),
  ("RT = json.load(open('/home/claude/s1spec/step3nl_v33/verifier/rate_table_fin12f4_a0.00015.json'))",
   "RT = json.load(open(CFG.RT))"),
 ],
 't2_A/rigor/assembly_rep.py': [
  ("sys.path.insert(0, '/home/claude/s1spec/writeup/lemmaU/bulk/fix4_r3/r3'); import r3_forcing as R3   # fix4_r3 (R3, v36pre): honest forcing (lambda~ + p/4 != 0, eps-part); see r3_forcing.py",
   "sys.path.insert(0, '%s/r3'); import r3_forcing as R3   # fix4_r3 (R3, v36pre): honest forcing (lambda~ + p/4 != 0, eps-part); see r3_forcing.py; v36: gated by V36_R3\nCFG = L.CFG" % N),
  ("    r0 = R3.r0_honest(BETA, AL, ABAR, ABAR)          # R3: was 1.4(betabar + 1.02 alpha_1 abar) (derived with lambda~ = -p/4); abar = bbar (verifier)\n"
   "    assert AS.GAM == 1; FE2 = [R3.feps2(j, 14, Rout, 28, AL) for j in range(3)]   # R3: eps-part of the forcing, int_{14}^{Rout} omega |d_s^j R_eps|^2 ds",
   "    if CFG.R3ON:\n"
   "        r0 = R3.r0_honest(BETA, AL, ABAR, ABAR)          # R3: was 1.4(betabar + 1.02 alpha_1 abar) (derived with lambda~ = -p/4); abar = bbar (verifier)\n"
   "        assert AS.GAM == 1; FE2 = [R3.feps2(j, 14, Rout, 28, AL) for j in range(3)]   # R3: eps-part of the forcing, int_{14}^{Rout} omega |d_s^j R_eps|^2 ds\n"
   "    else:\n"
   "        r0 = D('1.4')*(BETA + D('1.02')*AL*ABAR); FE2 = None                          # v35 (V36_R3=0, regression only)"),
  ("        F2 = R3.F_honest(j, r0, Iint, FE2[j])**2      # R3: F_j <= r0' p^j sqrt(Iint) + sqrt(FE2_j) (Minkowski); was (r0 p^j)^2 Iint",
   "        F2 = R3.F_honest(j, r0, Iint, FE2[j])**2 if CFG.R3ON else (r0*P**j)**2*Iint      # R3: F_j <= r0' p^j sqrt(Iint) + sqrt(FE2_j) (Minkowski); was (r0 p^j)^2 Iint"),
  ("    return dict(E=E, Eu=Eu, S0=S0, S1=S1, parts=parts, gp=gps, info=info, cmin=cmin, cbar=cbar, r0=r0, Iint=Iint, FE2=FE2)",
   "    return dict(E=E, Eu=Eu, S0=S0, S1=S1, parts=parts, gp=gps, info=info, cmin=cmin, cbar=cbar, r0=r0, Iint=Iint, **({'FE2': FE2} if CFG.R3ON else {}))"),
  ("    for j in (0, 1):\n        gp = R['gp'][j]; inf = R['info'][j]; Dj = inf['D']",
   "    for j in ((0, 1, 2) if REXP[2].upper() > 0 else (0, 1)):     # v36 (G1/F2', r_2): level 2 exported when weights_t2 retains r_2 > 0 on [28, 350]\n"
   "        gp = R['gp'][j]; inf = R['info'][j]; Dj = inf['D']"),
  ("    _V = json.load(open('/home/claude/s1spec/step3nl_v33/verifier/modsups_fin12f4/final_alpha_points_5e-8.json'))['1.5e-4']",
   "    _V = json.load(open(CFG.MODS))['1.5e-4']                                                         # v36: per V36_VER (was fin12f4)"),
  ("    cF = R3.F_honest(1, R['r0'], R['Iint'], R['FE2'][1])/AS.SQ        # R3: F_1/sqrt(Lambda) with the honest F_1 (was 1.4(...) p sqrt(Iint)/SQ)",
   "    cF = R3.F_honest(1, R['r0'], R['Iint'], R['FE2'][1])/AS.SQ if CFG.R3ON else D('1.4')*(AS.BETA + D('1.02')*AS.AL*AS.ABAR)/AS.SQ*AS.P*sqrt(R['Iint'])   # R3 (gated)"),
  ("        R3REC = dict(r0=R['r0'], FE2=R['FE2'], F2=[R['info'][j]['F2'] for j in range(3)], module=R3.__file__)\n"
   "        out[name] = dict(r3=R3REC, verifier_sups=",
   "        R3REC = dict(r0=R['r0'], FE2=R['FE2'], F2=[R['info'][j]['F2'] for j in range(3)], module=R3.__file__) if CFG.R3ON else None\n"
   "        out[name] = dict(**({'r3': R3REC} if CFG.R3ON else {}), verifier_sups="),
  ("L.BASE + '/lemmaS_cflux.json', '/home/claude/s1spec/step3nl_v33/verifier/rate_table_fin12f4_a0.00015.json', os.path.abspath(__file__), os.path.abspath(R3.__file__),",
   "L.BASE + '/lemmaS_cflux.json', CFG.RT, CFG.MODS, os.path.abspath(__file__), os.path.abspath(R3.__file__), L.CFG.__file__,"),
 ],
 't2_A/rigor/weights_t2.py': [
  ("        if j < 2:\n            fr = 2 - rtgt - used_c/amin_c; bc_c = Q_(math.floor(float(fr.lower())*100), 100)\n"
   "        else:\n            fr = 2 - used_c/amin_c; bc_c = Q_(math.floor(float(fr.lower())*1000), 1000)",
   "        if j < 2:\n            fr = 2 - rtgt - used_c/amin_c; bc_c = Q_(math.floor(float(fr.lower())*100), 100)\n"
   "        else:   # v36 (r_2): level 2 retains r_2 = 1/1000 (units amin_c) when V36_R2 = 1; v35: r_2 = 0\n"
   "            fr = 2 - A(R2Q) - used_c/amin_c; bc_c = Q_(math.floor(float(fr.lower())*1000), 1000)"),
  ("        rreq_c = (rr/amin_c) if j < 2 else arb(0); rreq_s = (rr/amin_s) if j < 2 else arb(0)",
   "        rreq_c = (rr/amin_c) if j < 2 else A(R2Q); rreq_s = (rr/amin_s) if j < 2 else A(R2Q)*amin_c/amin_s   # v36: absolute r_2 amin_c on [28, 350] at level 2"),
  ("    r_exp = [Q_(math.floor(float(out['levels'][str(j)]['core']['retained'].lower())*1000 - 1e-9), 1000) for j in range(2)] + [Q_(0)]",
   "    r_exp = [Q_(math.floor(float(out['levels'][str(j)]['core']['retained'].lower())*1000 - 1e-9), 1000) for j in range(2)] + [R2Q]\n"
   "    assert A(R2Q) <= out['levels']['2']['core']['retained'].lower()          # v36: the exported r_2 is below the certified core retained fraction"),
  ("import weights_rep as W0\n",
   "import weights_rep as W0\nR2Q = Q_(L.CFG.R2_EXPORT) if L.CFG.R2ON else Q_(0)   # v36 (ROUTE_DECISION.md G1/F2'): level-2 retained dissipation r_2 = 1/1000\n"),
 ],
 't2_A/rigor/export_t2.py': [],   # handled separately (larger insertion), see below
 't4_cert/t4lib.py': [
  ("T5 = ROOT + '/bulk/fix4_r3/t5_zone3'; T2 = ROOT + '/bulk/fix4_r3/t2_A'   # fix4_r3 (R3) copy: paths only",
   "T5 = ROOT + '/bulk/fix6_v36/t5_zone3'; T2 = ROOT + '/bulk/fix6_v36/t2_A'   # fix6_v36 (v36) copy\n" + CFGIMP),
  ("VERIF_ALPHA = '/home/claude/s1spec/step3nl_v33/verifier/modsups_fin12f4/final_alpha_points_5e-8.json'   # g13: fin12f4 (was step3nl/verifier fin11; six-digit values identical)",
   "VERIF_ALPHA = CFG.MODS   # v36: per V36_VER (g13 was fin12f4)"),
  ("             env={k: os.environ.get(k) for k in ('PYTHONHASHSEED', 'PYTHONDONTWRITEBYTECODE', 'OMP_NUM_THREADS')},",
   "             env={k: os.environ.get(k) for k in ('PYTHONHASHSEED', 'PYTHONDONTWRITEBYTECODE', 'OMP_NUM_THREADS') + CFG.ENVKEYS},"),
 ],
 't4_cert/u34_sweeps.py': [
  ("sys.path.insert(0, '/home/claude/s1spec/writeup/lemmaU/bulk/fix4_r3/r3'); import r3_forcing as R3   # fix4_r3 (R3, v36pre): DL = |lambda~ + p/4|",
   "sys.path.insert(0, '%s/r3'); import r3_forcing as R3   # fix4_r3 (R3, v36pre): DL = |lambda~ + p/4|; v36: gated by V36_R3" % N),
  ("cdU = A(Q_(13056, 10000))*(BETAB + A(R3.DL_Q)*ALPHA1 + P*ABAR*ALPHA1)/(P*CM)   # R3: betabar' = betabar + |lambda~ + p/4| alpha_1",
   "cdU = A(Q_(13056, 10000))*(BETAB + (A(R3.DL_Q)*ALPHA1 if L.CFG.R3ON else 0) + P*ABAR*ALPHA1)/(P*CM)   # R3: betabar' = betabar + |lambda~ + p/4| alpha_1 (gated)"),
 ],
 't5_zone3/zone3.py': [
  ("_VS = json.load(open('/home/claude/s1spec/step3nl_v33/verifier/modsups_fin12f4/final_alpha_points_5e-8.json'))['1.5e-4']   # M14 (M3): abar read from the verifier",
   "sys.path.insert(0, '%s'); import v36cfg as CFG\n" % N +
   "_VS = json.load(open(CFG.MODS))['1.5e-4']   # M14 (M3): abar read from the verifier; v36: per V36_VER"),
  ("assert abar.upper() <= arb('2.19103e-7').upper() and abar.lower() >= arb('2.19103e-7').lower()",
   "assert abar.upper() <= arb('2.19103e-7').upper()   # v36: <= the v35 value (fin13v4b: 2.14041e-7); q, cminus, Db are monotone in abar"),
 ],
 'r3/r3_forcing.py': [
  ("    BASE = '/home/claude/s1spec/writeup/lemmaU/bulk/fix4_r3/t2_A/base_ibp'\n"
   "    FINP = BASE + '/lemmaS_F10p_M0.52.json'; VA = '/home/claude/s1spec/step3nl_v33/verifier/modsups_fin12f4/final_alpha_points_5e-8.json'",
   "    BASE = '%s/t2_A/base_ibp'\n" % N + "    " + CFGIMP + "\n"
   "    FINP = BASE + '/lemmaS_F10p_M0.52.json'; VA = CFG.MODS\n"
   "    assert 'LAMT = -0.2561532;' in open(CFG.CONSTPY).read()          # v36: LAMT of the verifier in use equals LAMT_Q"),
  ("    ins = [os.path.abspath(__file__), FINP, VA]", "    ins = [os.path.abspath(__file__), FINP, VA, CFG.CONSTPY, CFG.__file__]"),
 ],
 'n1/n1_chart.py': [
  ("VER = '/home/claude/s1spec/step3nl_v33/verifier/rate_table_fin12f4_a0.00015.json'",
   CFGIMP + "\nVER = CFG.RT                                         # v36: verifier per V36_VER (was fin12f4)"),
  ("               label='certified (arb, outward); hypothesis |alpha| <= alpha_1 = 1.5e-4 (v3.3 bootstrap), verifier fin12f4 sups (alpha_1 row)',",
   "               label='certified (arb, outward); hypothesis |alpha| <= alpha_1 = 1.5e-4 (v3.3 bootstrap), verifier %s sups (alpha_1 row)' % CFG.VER,"),
  ("    dst = sys.argv[sys.argv.index('--out') + 1] if '--out' in sys.argv else f'{OV}/n1_chart.json'      # --out: test runs only",
   "    dst = sys.argv[sys.argv.index('--out') + 1] if '--out' in sys.argv else f'{CFG.N}/n1/out/n1_chart_{CFG.VER}.json'   # v36: never writes into step4"),
  ("                         inputs={p: sha(p) for p in [os.path.abspath(__file__), VER] + MODS}))",
   "                         inputs={p: sha(p) for p in [os.path.abspath(__file__), VER] + MODS}, v36_env=CFG.env_record()))"),
 ],
}

EXPORT_T2 = [
 ("        Dsrc = [R['info'][k]['D'] for k in range(2)]; F2 = [R['info'][k]['F2'] for k in range(2)]; coup = [R['info'][k]['coup_rate'] for k in range(2)]",
  "        NX = 3 if '2' in R['export'] else 2      # v36: level-2 export present iff r_2 > 0 (V36_R2 = 1)\n"
  "        Dsrc = [R['info'][k]['D'] for k in range(NX)]; F2 = [R['info'][k]['F2'] for k in range(NX)]; coup = [R['info'][k]['coup_rate'] for k in range(NX)]"),
 ("            Dsrc=vec(Dsrc), r_export=W['r_export'][:2], amin_c=num(acore, 'lo'), rhoJ=num(rhoJ), F2=vec(F2), coup=vec(coup),",
  "            Dsrc=vec(Dsrc), r_export=W['r_export'][:NX], amin_c=num(acore, 'lo'), rhoJ=num(rhoJ), F2=vec(F2), coup=vec(coup),"),
 ("            C0_ref=vec([exp[str(k)]['C0'] for k in range(2)]), C1_0=vec([exp0[str(k)]['C1'] for k in range(2)]), C1_ref=vec([exp[str(k)]['C1'] for k in range(2)]),",
  "            C0_ref=vec([exp[str(k)]['C0'] for k in range(NX)]), C1_0=vec([exp0[str(k)]['C1'] for k in range(NX)]), C1_ref=vec([exp[str(k)]['C1'] for k in range(NX)]),"),
 ("        out['res'][st] = res\n",
  "        if NX == 3:   # v36 (ROUTE_DECISION.md G1/F2'): level-2 window export, consumer = flux lemma v2 (Agmon trace / endpoint term)\n"
  "            res['r2_export'] = r2_block(R, R4, R0, exp, exp0, acore, W)\n"
  "        out['res'][st] = res\n"),
 ("def main():\n",
  "def r2_block(R, R4, R0, exp, exp0, acore, W):\n"
  "    \"\"\"v36: the level-2 dissipation export of (A) (k = 2).  Same Lemma D window form as k = 0, 1 (assembly_rep.export):\n"
  "    X_2(T) := int_{tau_a}^{tau_a+T} int_{ln 28}^{ln 350} y^-4 (28/y) z_4^2 ds dtau <= C0_2 + C1_2 T, any tau_a >= tau_1, T >= 0,\n"
  "    C0_2 = (EA_2 + D_2)/(r_2 amin_c),  C1_2 = (2 F_2^2/gp_2 + coup_2 + o_2 + D_2^slope)/(r_2 amin_c),  r_2 = 1/1000.\"\"\"\n"
  "    i2 = R['info'][2]\n"
  "    return dict(k=2, z='z_4 = d_s^4 (f - U), s = ln y', r_2=W['r_export'][2], amin_c=num(acore, 'lo'),\n"
  "                statement='int_{tau_a}^{tau_a+T} int_{ln28}^{ln350} y^-4 (28/y) z_4^2 ds dtau <= C0_2 + C1_2*T  (all tau_a >= tau_1, T >= 0)',\n"
  "                C0_ref=num(exp['2']['C0']), C1_ref=num(exp['2']['C1']), C1_0=num(exp0['2']['C1']), C0_ref4=num(R4['export']['2']['C0']), C1_ref4=num(R4['export']['2']['C1']),\n"
  "                C0_rule='(EA_2 + D_2)/(r_2 amin_c)', C1_rule='(2 F_2^2/gp_2 + coup_2 + o_2 + D_2^slope)/(r_2 amin_c);  C1_2(Nbar) = C1_0 + o_2/(r_2 amin_c)',\n"
  "                window_data=dict(EA_2=num(R['E'][2]), D_2=num(i2['D']), D_window_const=num(i2['D_window_const']), D_window_slope=num(i2['D_window_slope']),\n"
  "                                 F2_2=num(i2['F2']), gp_2=num(R['gp'][2], 'lo'), coup_rate_2=num(i2['coup_rate']), rho_w=num(i2['rho_w'])),\n"
  "                retained_core=num(W['levels']['2']['core']['retained'], 'lo'), Bc_core=W['levels']['2']['core']['Bc'],\n"
  "                note='r_2 is retained at level 2 on all of [28, 350] (weights_t2: core reference, zones [40,60], segments, transition [250,350]); '\n"
  "                     'every level-2 source in assembly_rep is a Tr (damped, window const, window slope), as at levels 0, 1; FOR THE FLUX-LEMMA DESIGNER: see V36_REPORT.md section r_2')\n"
  "\n"
  "def main():\n"),
]

def apply(path, edits):
    p = os.path.join(N, path); s = open(p).read()
    for old, new in edits:
        c = s.count(old)
        if c != 1: sys.exit('edit not unique/missing (%d) in %s: %r' % (c, path, old[:120]))
        s = s.replace(old, new)
    open(p, 'w').write(s); print('edited', path, len(edits))

for path, ed in EDITS.items():
    if ed: apply(path, ed)
apply('t2_A/rigor/export_t2.py', EXPORT_T2)
