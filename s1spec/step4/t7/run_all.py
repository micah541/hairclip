"""T7: the whole production pipeline in one command, on a given polished store.
  python run_all.py --base ../t1/polished_v2 --work run_v2 --jobs 2 [--steps all] [--mu 0.036] [--publish] [--t8]
  python run_all.py --base ../t1/polished_v2 --work run_v2se --se ../p2 [--prev run_v2] ...   (far constants on the (SE)-corrected
         base: se7.py widths from p2/out/se_corr.json; far7/nl7 --se, farg increments -> out_farg_se; --prev reuses the
         unchanged far cells of nl7 from an earlier run on the polished base)
Slabs: the whole store (its manifest.json; polished_v2: 0-180), or --slabs a-b (v1 store: 5-187).  Certificate intervals:
base7.schedule on the store's own tau axis (polished_v2: 26 intervals).
Steps (each restartable: units whose JSON exists are skipped):
  tip    tip7.py   E2 Garding, strip, ramp constants (per slab)
  far    far7.py   mu1, mu2, a_mu on the whole far curve -> lambda < 2 mu1 (per slab)
  cert   cert7.py  (C1) + (C2) local margins, reserved rank-3 budget (one JSON per certificate interval)
  c2r3   c2r3.py   rank-3 variation V per slab (verifies the reserves) and (C3)
  farg   farg7.py  modulation Gram matrices (C_gamma0, C_gamma2, G_x), one tc-ball per slab
  fargref farg7.py --refine 0.1 --nsub 4: 4 tc-balls on the slabs with sqrt(lmax Gamma0) > 0.1 (before outer insertions)
  nl     nl7.py    n2-n4: tip (box on e, e' only; e'' free, cubic m23), T5 far box, sharp far box
  merge  merge.py  -> WORK/constants.json;  closure_cert.py (certified rho_c at mu_e) -> WORK/closure_cert.json; merge again
  t8     (with --t8 and --publish; needs T6 on the same slabs, t6/residuals.json + t6/jets) T8's cconst.py (per-slab
         constants, recomputed: ~1 CPU-s per slab), merge.py, lambda < 2 mu1 with T8's certified rates (T7 merge re-run with
         --lambda_max, published), camp.py (window: trailing slabs with k >= 2.0e4); prints the Ybar requirement
Portable: subprocesses via sys.executable, os.path only."""
import os, sys, json, argparse, subprocess, time
HERE = os.path.dirname(os.path.abspath(__file__))
PY = sys.executable
SLABS = None


def run(args, log=None):
    print('>>', ' '.join(os.path.basename(a) if i == 1 else a for i, a in enumerate(args)), flush=True)
    t0 = time.time(); r = subprocess.run(args, cwd=HERE, stdout=log, stderr=subprocess.STDOUT if log else None)
    print(f'   exit {r.returncode} ({time.time() - t0:.0f}s)', flush=True); return r.returncode


def batch(script, out, base, jobs, logdir, extra=()):
    with open(os.path.join(logdir, script.replace('.py', '.log')), 'a') as lg:
        return run([PY, os.path.join(HERE, 'run_batch.py'), '--script', script, '--slabs', SLABS, '--jobs', str(jobs), '--out', out, '--base', base]
                   + (['--'] + list(extra) if extra else []), lg)


def certs(out, base, jobs, logdir, ids):
    groups = [ids[j::jobs] for j in range(jobs)]; groups = [g for g in groups if g]
    procs = []
    for j, g in enumerate(groups):
        lg = open(os.path.join(logdir, f'cert7_{j}.log'), 'a')
        cmd = [PY, '-u', os.path.join(HERE, 'cert7.py'), '--times', ','.join(map(str, g)), '--out', out, '--base', base, '--nbis', '6',
               '--slabs', SLABS]
        procs.append((subprocess.Popen(cmd, cwd=HERE, stdout=lg, stderr=subprocess.STDOUT), lg))
    for p, lg in procs: p.wait(); lg.close()
    return max(p.returncode for p, _ in procs)


SEFLAG = []


def t8_step(W, base, jobs):
    """T8 on the published T7 constants (T6 must cover the same slabs)"""
    T8 = os.path.normpath(os.path.join(HERE, '..', 't8')); T6 = os.path.normpath(os.path.join(HERE, '..', 't6'))
    a_, b_ = [int(v) for v in SLABS.split('-')]
    R6 = json.load(open(os.path.join(T6, 'residuals.json')))
    if list(R6['slabs_range']) != [a_, b_] or R6.get('missing'):
        print(f"t8 skipped: t6/residuals.json covers {R6['slabs_range']} (missing {R6.get('missing')}), not {SLABS}", flush=True); return 1
    R = json.load(open(os.path.join(W, 'closure_cert.json'))); q = json.load(open(os.path.join(HERE, 't8req.json')))
    q['P=2']['best']['muE'] = R['mu_e']; q['P=2']['best']['rho'] = R['rho_c_min']
    q['P=2']['best']['source_muE_rho'] = f'T7 certified closure ({os.path.relpath(W, HERE)}/closure_cert.json)'
    json.dump(q, open(os.path.join(HERE, 't8req.json'), 'w'), indent=1)
    sd = os.path.join(T8, 'constants', 'slabs')
    if os.path.isdir(sd) and os.listdir(sd):                     # per-slab constants of an earlier run: set aside, recompute
        os.replace(sd, sd + '_prev_' + time.strftime('%Y%m%d_%H%M%S'))
    run([PY, os.path.join(T8, 'cconst.py'), '--slabs', SLABS, '--out', sd, '--jobs', str(jobs)])
    run([PY, os.path.join(T8, 'merge.py'), '--in', sd, '--first', str(a_), '--last', str(b_)])
    C8 = json.load(open(os.path.join(T8, 'constants', 'constants.json')))
    if C8['missing']: print('t8: T8 per-slab constants missing for', C8['missing'], flush=True); return 1
    lam_max = max(p['lev'][j]['lam'] for p in C8['per_slab'] for j in p['lev'])
    run([PY, os.path.join(HERE, 'merge.py'), '--work', W, '--base', base, '--slabs', SLABS, '--lambda_max', repr(lam_max), '--publish'] + SEFLAG)
    run([PY, os.path.join(T8, 'camp.py')])
    C7 = json.load(open(os.path.join(HERE, 'constants.json'))); S7 = C7['slabs']
    mu1 = min(v['E_I']['mu1'] for v in S7.values() if 'E_I' in v)
    Cg0 = max(v['modulation']['Cgamma0'] for v in S7.values() if 'modulation' in v)
    a = json.load(open(os.path.join(T8, 'constants', 'amp.json'))); h = a['handoff']
    print(f"lambda < 2 mu1: T8 certified lambda_max = {lam_max:.4f} < 2 min mu1 = {2 * mu1:.4f} (margin {2 * mu1 / lam_max - 1:.3f})", flush=True)
    print(f"rho_c min = {C7['meta']['closure']['rho_c_min']} (mu_e = {C7['meta']['closure']['mu_e']}), K_e = {C7['meta']['K_e']['value']:.4g}, "
          f"C_gamma0 max = {Cg0:.4g}; T8 inputs: " + ', '.join(f"{k}={C8['t7'][k]}" for k in ('muE', 'rho', 'Ke', 'Cg0')), flush=True)
    print(f"T8 Ybar requirement: {h['Ymax']:.4e} (binding {h['binding']}), window slabs {h['window_slabs']} (k from {h['window_k_start'][0]:.4g})", flush=True)
    return 0


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--base', required=True); ap.add_argument('--work', required=True)
    ap.add_argument('--jobs', type=int, default=2); ap.add_argument('--steps', default='all'); ap.add_argument('--mu', default='0.036')
    ap.add_argument('--publish', action='store_true'); ap.add_argument('--t8', action='store_true')
    ap.add_argument('--se', default=None, help='p2 directory: far constants on the (SE)-corrected base')
    ap.add_argument('--prev', default=None, help='work dir of an earlier polished-base run (nl7 reuses its unchanged far cells)')
    ap.add_argument('--slabs', default=None, help='a-b: slabs of P2 (default: the whole store; v1 store: 5-187)')
    ap.add_argument('--cert_ids', default=None, help='(testing) certificate ids, e.g. 0-3 or 1,5 (default: all of the plan)')
    ap.add_argument('--nlclosure', default=None, help='T8 n3.json: nonlinear (Y) closure (T7 Lemma 4c) in the merge step')
    ap.add_argument('--ybox', default='1.7e-12', help='bootstrap box of the nonlinear closure (the n3.json run Y)')
    A = ap.parse_args()
    W = A.work if os.path.isabs(A.work) else os.path.normpath(os.path.join(HERE, A.work)); os.makedirs(W, exist_ok=True)
    logdir = os.path.join(W, 'logs'); os.makedirs(logdir, exist_ok=True)
    base = A.base if os.path.isabs(A.base) else os.path.normpath(os.path.join(HERE, A.base))
    import run_batch, base7
    global SLABS; SLABS = A.slabs or run_batch.store_range(base)
    a_, b_ = [int(v) for v in SLABS.split('-')]
    plan = base7.schedule(base, a_, b_)
    missing = [s for s in range(a_, b_ + 1) if not os.path.exists(os.path.join(base, 'slab_%03d.pkl' % s))]
    if missing: raise SystemExit(f'polished store {base} lacks slabs {missing[:5]}...')
    steps = ['sewidths', 'tip', 'far', 'cert', 'c2r3', 'farg', 'fargref', 'fargse', 'nl', 'merge', 't8'] if A.steps == 'all' else A.steps.split(',')
    SE = ['--se'] if A.se else []
    global SEFLAG; SEFLAG = SE
    if A.se and 'sewidths' in steps:
        run([PY, os.path.join(HERE, 'se7.py'), 'widths', '--p2', A.se if os.path.isabs(A.se) else os.path.normpath(os.path.join(HERE, A.se))])
    json.dump(dict(base=base, work=W, slabs=SLABS, plan=plan, steps=steps, mu=A.mu, started=time.ctime()), open(os.path.join(W, 'run_all.json'), 'w'), indent=1)
    print(f'store {base}: slabs {SLABS}, {len(plan)} certificate intervals, steps {steps}', flush=True)
    o = lambda d: os.path.join(W, d)
    if 'tip' in steps: batch('tip7.py', o('out_tip'), base, A.jobs, logdir)
    if 'far' in steps: batch('far7.py', o('out_far'), base, A.jobs, logdir, SE)
    ci = ([T['id'] for T in plan] if A.cert_ids is None else
          [int(v) for v in (range(int(A.cert_ids.split('-')[0]), int(A.cert_ids.split('-')[1]) + 1) if '-' in A.cert_ids else A.cert_ids.split(','))])
    if 'cert' in steps: certs(o('out'), base, A.jobs, logdir, ci)
    if 'c2r3' in steps: batch('c2r3.py', o('out_c2r3'), base, A.jobs, logdir)
    if 'farg' in steps: batch('farg7.py', o('out_farg'), base, A.jobs, logdir)
    if 'fargref' in steps: batch('farg7.py', o('out_farg'), base, A.jobs, logdir, ['--refine', '0.1', '--nsub', '4'])
    if A.se and 'fargse' in steps: run([PY, os.path.join(HERE, 'se7.py'), 'farg', '--in', o('out_farg'), '--out', o('out_farg_se')])
    prev_nl = (['--prev', os.path.join(A.prev if os.path.isabs(A.prev) else os.path.normpath(os.path.join(HERE, A.prev)), 'out_nl')]
               if A.se and A.prev else [])
    if 'nl' in steps: batch('nl7.py', o('out_nl'), base, A.jobs, logdir, SE + prev_nl)
    if 'merge' in steps:
        mg = [PY, os.path.join(HERE, 'merge.py'), '--work', W, '--base', base, '--slabs', SLABS] + SE
        try:                                  # keep T8's certified lambda* of an earlier t8 step (else merge's default 0.2527)
            lm = json.load(open(os.path.join(W, 'constants.json')))['meta'].get('lambda_max_T8')
            if lm: mg += ['--lambda_max', repr(lm)]
        except (OSError, ValueError, KeyError): pass
        run(mg)
        cc = [PY, os.path.join(HERE, 'closure_cert.py'), '--work', W] + ([] if A.mu == 'scan' else ['--mu', A.mu])
        run(cc)
        if A.nlclosure:                       # T7 Lemma 4c: nonlinear terms of (Y) on T8's far box (n3.json), then nonlinear mode
            run(mg)
            mu_ = json.load(open(os.path.join(W, 'closure_cert.json')))['mu_e']
            run([PY, '-I', os.path.join(HERE, 'nlterms7.py'), '--work', W, '--n3', A.nlclosure, '--Ybox', A.ybox])
            run([PY, os.path.join(HERE, 'closure_cert.py'), '--work', W, '--mu', repr(mu_), '--nl', os.path.join(W, 'nl_terms.json'),
                 '--out', 'closure_cert_nl.json'])
        run(mg + (['--publish'] if A.publish else []))
    if 't8' in steps and A.t8 and A.publish: t8_step(W, base, A.jobs)
    print('run_all done', time.ctime())


if __name__ == '__main__':
    main()
