#!/usr/bin/env python3
"""Reproduction kit B v3: the rest of step 4 (v3.2 chain) plus the 21 v3.3 hand-off stages, independent of kit A (starts from the reference P1 slab-95 checkpoint shipped in step4/t3/p1/ckpt).

    python run_all.py                 run (or RESUME) every stage; ends with check.py.  Safe to interrupt: run the same command again.
    python run_all.py --jobs 6        worker processes for the parallel stages (default 6; each process is single-threaded)
    python run_all.py --list          list the stages with their state
    python run_all.py --only a,b      run only these stages;   --from NAME  start at NAME;   --force a,b  redo these stages
    python run_all.py --status | --clean

Stages (sequential; inside a stage the work is spread over --jobs processes; a finished stage leaves kit_out/state/<stage>.done):
  inputs      SHA-256 of the shipped data inputs (P1 slabs 95/96, S1 CAPD profiles, step-3 hull data, Lemma O / S constants)
  t1_seed, t1_float, t1_polish, t1_se   base v2 (seed from P1 slab 95, float base, arb-128 polish, (SE) offsets)
  t6, rhoend, se_corr                   T6 residual enclosures of the polished base; far-end residuals; (SE) correction (used for slabs 0-1 only)
  t6c         B1: T6's engine on the (SE)-corrected base itself (residuals.json + jets, slabs 2-180)
  t7_a, t7_b  T7 certificates on the polished base (+T8), then on the corrected base (se7 on t6c; published constants)
  h1pts       the 16 exact H1 points (h1.py), the input of h1dt's summary           [kit fix 2]
  h1dt, winslab   slab-180 hand-off: 32 derivative balls, 32 window/shell balls (+ their summaries)
  tail1       round 1 of the final tail: the pre stages and the 'T8 ...' stages of pass A only (they make n3_b1, all t7_nl needs)   [kit fix 5]
  t7_nl       T7 republished with the nonlinear closure (Lemma 4c) on tail1's n3_b1
  tail2       the final tail again on the republished T7 (the reference chain: p2/out/FINAL_manifest.json), incl. the B8 check as its last stage
  b8          the B8 hull check (check_outA5.py) once more, standalone, in the kit's own s1_tight/consistency   [kit fix 6]
  v33_00 ... v33_20   the 21 stages of the v3.3 step-4 chain (reference: p2/out_v33/FINAL_manifest_v33.json of the proof tree): (H2f), tightened (H2),
              fixed-L cone family, B0, M9 (arb), N1 chart terms, f12_sliver.  Driver: step4/p2/v33/final_tail_v33.py (--from i --until i).
              Writes step4/p2/out_v33/FINAL_manifest_v33.json; check.py then checks the margins of STEP4_RESULT_v33_addendum and compares every output.
  v38         [v4] the 29-stage step-4 v38 chain (step4/p2/v38/final_tail_v38.py): carries v3.3 stages 0-7, 10-17 by hash from this run's v3.3 manifest, reruns
              8, 9, 18, 19, 20 on verifier fin13v4c / (A) chain v37, adds 21 flux v37, 22-26 C1_ke, GS_tau1, GS_boot, G1_check, U1_sym, 27 iface_compare, 28 v4c pins.
              Writes step4/p2/out_v38/FINAL_manifest_v38.json; check.py section v38.
Write guard [v4]: every process of the kit runs through kit_guard.py (writes / chdir / child cwd outside the kit raise).
Logs: kit_out/logs/<stage>.log"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); import kit_guard   # noqa: E402  [v4] write guard for this process and every child

FORCED_ENV = dict(OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1', NUMEXPR_NUM_THREADS='1',
                  VECLIB_MAXIMUM_THREADS='1', BLIS_NUM_THREADS='1', PYTHONHASHSEED='0', PYTHONUTF8='1', PYTHONIOENCODING='utf-8')
os.environ.update(FORCED_ENV)
os.environ.pop('PYTHONPATH', None)

import json, time, shutil, subprocess, platform, hashlib, glob

try: sys.stdout.reconfigure(encoding='utf-8', errors='replace'); sys.stderr.reconfigure(encoding='utf-8', errors='replace')
except Exception: pass

KIT = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(KIT, 'kit_out'); LOGS = os.path.join(OUT, 'logs'); STATE = os.path.join(OUT, 'state')
S4 = os.path.join(KIT, 'step4')
T1, T3, T6, T7, T8, P2 = (os.path.join(S4, d) for d in ('t1', 't3', 't6', 't7', 't8', 'p2'))
T8SE, T9H, POUT = os.path.join(P2, 't8se'), os.path.join(P2, 't9h'), os.path.join(P2, 'out')
PY = sys.executable
PINS = dict(flint='0.9.0', numpy='2.5.4', scipy='1.18.1', sympy='1.14.0', mpmath='1.3.0')
try: FM = json.load(open(os.path.join(KIT, 'reference', 'final_manifest.json'), encoding='utf-8'))      # the final tail of the proof tree (build.py)
except Exception as _e: print('ERROR: reference/final_manifest.json is missing or unreadable', repr(_e)); sys.exit(2)
S1TT = FM['h1env']['S1TT']

CFG = dict(iend=180, kend='2.31e4', ev='0,1,2,3,4,5,180', slabs='0-180', last=180, cert_ids=None)
JOBS = 6
IS_WIN = os.name == 'nt'


def hms(s):
    s = int(round(max(s, 0))); return f'{s // 3600}h{(s % 3600) // 60:02d}m' if s >= 3600 else f'{s // 60}m{s % 60:02d}s'
def now(): return time.strftime('%H:%M:%S')
def rel(p): return os.path.relpath(p, KIT)


# ------------------------------------------------------------------------------------------------------------ lock
_LOCK = None
def acquire_lock():
    global _LOCK
    os.makedirs(OUT, exist_ok=True)
    f = open(os.path.join(OUT, 'run.lock'), 'a+')
    try:
        if IS_WIN:
            import msvcrt; f.seek(0); msvcrt.locking(f.fileno(), msvcrt.LK_NBLCK, 1)
        else:
            import fcntl; fcntl.flock(f.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
    except OSError:
        print('Another run_all.py is already running in this folder (kit_out\\run.lock is locked). Not starting a second one.'); sys.exit(3)
    _LOCK = f


# ------------------------------------------------------------------------------------------------------------ process helpers
def child_env(extra=None):
    e = dict(os.environ); e.update(FORCED_ENV); e.pop('PYTHONPATH', None)
    if extra: e.update({k: str(v) for k, v in extra.items()})
    return e


def popen(cmd, cwd, env, logf):
    kw = {}
    if IS_WIN: kw['creationflags'] = getattr(subprocess, 'BELOW_NORMAL_PRIORITY_CLASS', 0)
    return subprocess.Popen([str(c) for c in cmd], cwd=cwd, env=child_env(env), stdout=logf, stderr=subprocess.STDOUT, **kw)


def last_line(path):
    try:
        with open(path, 'rb') as f:
            f.seek(0, 2); n = f.tell(); f.seek(max(0, n - 600)); txt = f.read().decode('utf-8', 'replace')
        lines = [l for l in txt.splitlines() if l.strip()]
        return lines[-1][:140] if lines else ''
    except OSError: return ''


class Stage:
    """execution context of one stage: log file, heartbeat, process pool"""
    def __init__(self, name):
        self.name = name; self.log = os.path.join(LOGS, name + '.log'); os.makedirs(LOGS, exist_ok=True); self.t0 = time.time(); self.progress = None
        self.last_beat = time.time()

    def header(self, cmd, cwd):
        with open(self.log, 'a', encoding='utf-8') as lf: lf.write(f'\n===== {time.strftime("%Y-%m-%d %H:%M:%S")}  {" ".join(str(c) for c in cmd)}  (cwd {rel(cwd)})\n')

    def beat(self, force=False):
        if force or time.time() - self.last_beat > 180:
            extra = self.progress() if self.progress else ''
            print(f'[{now()}]   {self.name}: running {hms(time.time() - self.t0)}  {extra or last_line(self.log)}', flush=True); self.last_beat = time.time()

    def run(self, cmd, cwd, env=None, ok_codes=(0,)):
        """one process, output appended to the stage log; returns the exit code"""
        self.header(cmd, cwd)
        with open(self.log, 'a', encoding='utf-8') as lf:
            p = popen(cmd, cwd, env, lf)
            try:
                while True:
                    try: rc = p.wait(timeout=10); break
                    except subprocess.TimeoutExpired: self.beat()
            except KeyboardInterrupt: interrupted(p)
        return rc

    def pool(self, tasks, jobs):
        """tasks: list of dict(cmd, cwd, env, skip=callable|None).  Up to `jobs` at once; returns number of failed tasks"""
        todo = [t for t in tasks if not (t.get('skip') and t['skip']())]
        running = []; failed = 0; total = len(todo); done = 0
        lf = open(self.log, 'a', encoding='utf-8')
        try:
            while todo or running:
                while todo and len(running) < jobs:
                    t = todo.pop(0); lf.write(f'\n===== {time.strftime("%Y-%m-%d %H:%M:%S")}  {" ".join(str(c) for c in t["cmd"])}  (cwd {rel(t["cwd"])})\n'); lf.flush()
                    running.append(popen(t['cmd'], t['cwd'], t.get('env'), lf))
                time.sleep(2)
                for p in list(running):
                    if p.poll() is not None:
                        running.remove(p); done += 1
                        if p.returncode != 0: failed += 1
                self.progress = lambda: f'{done}/{total} units done'
                self.beat()
        except KeyboardInterrupt:
            for p in running: p.terminate()
            interrupted(None)
        lf.close()
        return failed


def interrupted(p):
    print('\nInterrupted: stopping child processes. Everything finished so far is kept; run  python run_all.py  again to continue.', flush=True)
    if p is not None:
        try: p.wait(timeout=20)
        except Exception:
            p.terminate()
            try: p.wait(timeout=20)
            except Exception: p.kill()
    raise SystemExit(130)


def tail(path, n=25):
    try:
        with open(path, encoding='utf-8', errors='replace') as f: return ''.join(f.readlines()[-n:])
    except OSError: return '(no log)'


def split_range(a, b, n):
    """slabs a..b in n contiguous chunks 'x-y'"""
    tot = b - a + 1; n = max(1, min(n, tot)); out = []; s = a
    for i in range(n):
        m = tot // n + (1 if i < tot % n else 0); out.append((s, s + m - 1)); s += m
    return out


def ensure_dirs():
    for d in (LOGS, STATE, POUT, os.path.join(POUT, 'rhoend'), os.path.join(POUT, 't6se'), os.path.join(P2, 't6', 'jets'), os.path.join(T8SE, 'constants'),
              os.path.join(T8, 'constants'), os.path.join(T6, 'jets'), os.path.join(P2, 't6c', 'out'), os.path.join(P2, 't6c', 'jets'), os.path.join(POUT, 'h1slab')):
        os.makedirs(d, exist_ok=True)


# ------------------------------------------------------------------------------------------------------------ stages
def st_inputs(S):
    man = os.path.join(KIT, 'inputs', 'MANIFEST.sha256'); bad = []; n = 0
    with open(man, encoding='utf-8') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'): continue
            h, p = line.split('  ', 1); fp = os.path.join(KIT, *p.split('/')); n += 1
            if not os.path.exists(fp): bad.append((p, 'missing')); continue
            d = hashlib.sha256()
            with open(fp, 'rb') as g:
                for chunk in iter(lambda: g.read(1 << 20), b''): d.update(chunk)
            if d.hexdigest() != h: bad.append((p, 'SHA-256 differs'))
    msg = f'{n} input files checked, {len(bad)} problems'
    with open(S.log, 'a') as lf: lf.write(msg + '\n' + '\n'.join(f'{p}: {w}' for p, w in bad) + '\n')
    print(f'   {msg}' + ''.join(f'\n   PROBLEM {p}: {w}' for p, w in bad), flush=True)
    return not bad


def st_t1_seed(S): return S.run([PY, 'seed_p1.py', '95', '80'], T1) == 0 and os.path.exists(os.path.join(T1, 'seed_p1_95.pkl'))
def st_t1_float(S): return S.run([PY, 'run3.py', 'P11v2', CFG['kend'], 'seed_p1_95.pkl'], T1) == 0 and os.path.exists(os.path.join(T1, 'slabs_P11v2.pkl'))


def polish_count():
    try:
        with open(os.path.join(T1, 'polished_v2', 'manifest.json')) as f: return len(json.load(f))
    except Exception: return 0


def st_t1_polish(S):
    n = polish_count(); end = CFG['iend']
    if n >= end + 1: return True
    S.progress = lambda: f'slab {polish_count()} of {end + 1} polished'
    cmd = [PY, '-u', 'polish_v2.py', 'P11v2', str(end), CFG['ev'], '--seed', 'seed_p1_95.pkl', '--store', 'polished_v2']
    if n > 0:
        print(f'   resuming the polish at slab {n} (stored slabs 0..{n - 1} are reloaded losslessly)', flush=True); cmd += ['--resume', str(n)]
    return S.run(cmd, T1) == 0 and polish_count() >= end + 1


def st_t1_se(S):
    n = polish_count(); tasks = []; parts = []
    for a, b in split_range(0, n - 1, JOBS * 2):
        fn = os.path.join(T1, f'se_polished_v2_{a}_{b + 1}.json'); parts.append((a, fn))
        tasks.append(dict(cmd=[PY, 'se_v2.py', 'polished_v2', f'{a}:{b + 1}'], cwd=T1, skip=(lambda fn=fn: os.path.exists(fn))))
    if S.pool(tasks, JOBS): return False
    recs = []
    for a, fn in sorted(parts):
        with open(fn) as f: recs += json.load(f)
    recs.sort(key=lambda r: r['slab']); assert [r['slab'] for r in recs] == list(range(n))
    with open(os.path.join(T1, 'se_polished_v2.json'), 'w') as f: json.dump(recs, f, indent=1)
    return True


def st_t6(S):
    a, b = 0, CFG['last']
    if S.run([PY, '-u', 'run.py', '--slabs', f'{a}-{b}', '--out', 'out_v2', '--jobs', str(JOBS), '--poldir', os.path.join('..', 't1', 'polished_v2')], T6): return False
    return S.run([PY, '-u', 'merge.py', '--in', 'out_v2', '--first', str(a), '--last', str(b)], T6) == 0 and os.path.exists(os.path.join(T6, 'residuals.json'))


def st_rhoend(S):
    tasks = []
    for a, b in split_range(0, CFG['last'], JOBS * 3):
        tasks.append(dict(cmd=[PY, 'rho_end.py', f'{a}-{b}'], cwd=P2, skip=(lambda a=a, b=b: all(os.path.exists(os.path.join(POUT, 'rhoend', f'slab_{i:03d}.json')) for i in range(a, b + 1)))))
    return S.pool(tasks, JOBS) == 0


def st_se_corr(S): return S.run([PY, 'se_corr.py'], P2) == 0 and os.path.exists(os.path.join(POUT, 'se_corr.json'))


def t7_cmd(work, extra):
    c = [PY, '-u', 'run_all.py', '--base', os.path.join('..', 't1', 'polished_v2'), '--work', work, '--jobs', str(JOBS)] + extra
    if CFG['cert_ids']: c += ['--cert_ids', CFG['cert_ids']]
    return c


def t7_failed_steps(S, since):
    """T7's run_all.py returns 0 even when an inner step fails: scan its log (this session) for '   exit N' with N != 0"""
    try:
        with open(S.log, encoding='utf-8', errors='replace') as f: f.seek(since); txt = f.read()
    except OSError: return ['log unreadable']
    import re
    return [m.group(0) for m in re.finditer(r'^   exit [1-9]\d*.*$', txt, re.M)]


def run_t7(S, work, extra):
    since = os.path.getsize(S.log) if os.path.exists(S.log) else 0
    rc = S.run(t7_cmd(work, extra), T7)
    bad = t7_failed_steps(S, since)
    if bad: print(f'   T7 inner step(s) failed: {bad[:3]}', flush=True)
    return rc == 0 and not bad and os.path.exists(os.path.join(T7, 'constants.json'))


def st_t7_a(S):
    p = os.path.join(T7, 'constants.json')
    if os.path.exists(p): os.remove(p)                                        # a stale published file must not mask a failure
    return run_t7(S, 'run_v2', ['--publish'] + ([] if CFG['cert_ids'] else ['--t8']))      # quick mode: no hand-off window for T8's camp


def st_t7_b(S):
    # [kit fix 1] se7.py (source t6c) reads t7/se_widths_secorr.json for the slabs outside P2 (0-1): the widths of the se_corr source, kept under that name in the
    # proof tree (T7_REPORT.md:318), not shipped.  Regenerated here from this run's p2/out/se_corr.json:  se7.py widths --source se_corr, then renamed.
    # (se7.py widths writes t7/se_widths.json; the sewidths step of T7's run_all below overwrites it with the t6c-source widths.)
    wsc = os.path.join(T7, 'se_widths_secorr.json')
    if S.run([PY, 'se7.py', 'widths', '--source', 'se_corr', '--p2', P2], T7) != 0 or not os.path.exists(os.path.join(T7, 'se_widths.json')): return False
    os.replace(os.path.join(T7, 'se_widths.json'), wsc)
    a, b = os.path.join(T7, 'run_v2'), os.path.join(T7, 'run_v2se')
    os.makedirs(b, exist_ok=True)
    for d in ('out', 'out_c2r3', 'out_tip', 'out_farg'):                  # tip-region units are those of the polished-base run (T7_REPORT)
        if os.path.isdir(os.path.join(a, d)): shutil.copytree(os.path.join(a, d), os.path.join(b, d), dirs_exist_ok=True)
    return run_t7(S, 'run_v2se', ['--se', os.path.join('..', 'p2'), '--prev', 'run_v2', '--steps', 'sewidths,far,fargse,nl,merge', '--publish'])


def st_t6c(S):
    t6c = os.path.join(P2, 't6c'); a, b = 2, CFG['last']
    S.progress = lambda: f"{len(glob.glob(os.path.join(t6c, 'out', 'slab_*.json')))} of {b - a + 1} slabs"
    if S.run([PY, '-u', 'runc.py', '--slabs', f'{a}-{b}', '--jobs', str(JOBS)], t6c): return False
    if len(glob.glob(os.path.join(t6c, 'out', 'slab_*.json'))) < b - a + 1: print('   t6c: some slabs are missing', flush=True); return False
    return S.run([PY, '-u', 'mergec.py'], t6c) == 0 and os.path.exists(os.path.join(t6c, 'residuals.json'))


def st_b8(S):
    # [kit fix 6] B8 hull check, once more and standalone (the tail's own last stage runs it too, since [kit fix 7]); output: s1_tight/consistency/check_outA5.json
    d = os.path.join(KIT, 's1_tight', 'consistency')
    return S.run([PY, 'check_outA5.py'], d) == 0 and os.path.exists(os.path.join(d, 'check_outA5.json'))


# ---- slab-180 hand-off chain (p2/t9h): the run_h1slab.sh / run_winslab.sh of the proof tree, in Python
H1A = (1.04655427e-4, 0.99999491111, -1.17386340e-4); H1B = (1.00712248e-4, 0.99999510274, -1.12963590e-4)      # starting points of the Krawczyk iteration (run_h1slab.sh)


def h1_env(**kw):
    e = dict(COHERENT=1, S1PROF=os.path.join(KIT, 's1_tight', 'ext_tip', 'prof_tip.txt'), S1TT=S1TT); e.update(kw); return e


def tk_of(k): return -1 + (2 * k + 1) / 16


def st_h1pts(S):
    tasks = []
    for k in range(16):
        tk = tk_of(k); f = (tk + 1) / 1.9990234375; ref = ','.join(repr(x + f * (y - x)) for x, y in zip(H1A, H1B))
        out = os.path.join(POUT, f't9_h1_p{tk}_coh_s180.json')
        tasks.append(dict(cmd=[PY, 'h1.py', f'p{tk}'], cwd=T9H, env=h1_env(H1TAG='_s180', H1REF=ref), skip=(lambda out=out: os.path.exists(out))))
    return S.pool(tasks, JOBS) == 0 and all(os.path.exists(os.path.join(POUT, f't9_h1_p{tk_of(k)}_coh_s180.json')) for k in range(16))


DT_FINAL = os.environ.get('KIT_DT_POINTS') == 'final'
# The proof tree's dt balls (h1dt.py) were computed from the point enclosures that h1.py produced BEFORE the B5 rerun with the exact u4_cell hull
# (p2/out/h1slab/pre_u4); the points were rerun afterwards.  Both enclose the same zero p(t_k).  The kit reproduces the reference: dt reads the shipped
# earlier enclosures (reference/h1_pre_u4, SHA-256 in inputs/MANIFEST).  KIT_DT_POINTS=final uses the kit's own (later) points instead (not bit-identical).


def st_h1dt(S):
    tasks = []; outs = []
    for k in range(16):
        tk = tk_of(k); fpt = f't9_h1_p{tk}_coh_s180.json'
        pt = os.path.join('..', 'out', fpt) if DT_FINAL else os.path.join(KIT, 'reference', 'h1_pre_u4', fpt)
        for sg in (-1, 1):
            tc = tk + sg * 0.03125; fo = os.path.join(POUT, 'h1slab', f'dt_{tc}_0.03125.json'); outs.append(fo)
            tasks.append(dict(cmd=[PY, 'h1dt.py', str(tc), '0.03125', pt], cwd=T9H, env=h1_env(TREF=tk), skip=(lambda fo=fo: os.path.exists(fo))))
    if S.pool(tasks, JOBS) or not all(os.path.exists(o) for o in outs): return False
    return S.run([PY, 'h1slab_summaries.py', 'h1', POUT], T9H) == 0


def st_winslab(S):
    cmds = os.path.join(POUT, 'h1slab', 'win_cmds.txt')
    with open(cmds, 'w', encoding='utf-8') as f:
        S.header([PY, 'win_slab.py'], T9H)
        p = subprocess.run([PY, 'win_slab.py'], cwd=T9H, env=child_env(), stdout=f, stderr=subprocess.PIPE)
    if p.returncode: print(p.stderr.decode('utf-8', 'replace')[-600:]); return False
    lines = [l.strip() for l in open(cmds, encoding='utf-8') if l.strip().startswith('T9TC=')]
    if len(lines) != 32: print(f'   win_slab printed {len(lines)} commands, expected 32'); return False
    tasks = []
    for l in lines:
        env = dict(kv.split('=', 1) for kv in l.split()); fo = env['WINOUT']
        tasks.append(dict(cmd=[PY, 'win.py'], cwd=T9H, env=env, skip=(lambda fo=fo: os.path.exists(fo))))
    if S.pool(tasks, JOBS): return False
    return S.run([PY, 'h1slab_summaries.py', 'win', POUT], T9H) == 0


# ---- the final tail (reference/final_manifest.json = p2/out/FINAL_manifest.json of the proof tree)
def tr(v):
    v = str(v)
    import re   # [kit fix 3] single pass, so a KIT under /home/claude/s1spec is not rewritten twice
    M = {'/home/claude/s1spec/step4': S4.replace('\\\\', '/'), '/home/claude/s1_tight': os.path.join(KIT, 's1_tight').replace('\\\\', '/'), '/home/claude/s1spec': KIT.replace('\\\\', '/')}   # provenance strings of the manifest -> this kit
    return re.sub('/home/claude/s1spec/step4|/home/claude/s1_tight|/home/claude/s1spec', lambda m: M[m.group(0)], v)


def tail_run(S, st, drop_env=()):
    """one manifest stage: command, cwd, env, log file exactly as in the proof tree (log at its manifest path)"""
    cmd = [PY]; a = st['cmd']
    for i, x in enumerate(a):
        cmd.append(str(JOBS) if i and a[i - 1] == '--jobs' else tr(x))
    env = {k: tr(v) for k, v in st['env'].items() if k not in drop_env}
    if 'H1JOBS' in env: env['H1JOBS'] = str(JOBS)
    cwd = os.path.join(S4, st['cwd'][len('step4/'):]) if st['cwd'].startswith('step4/') else os.path.join(S4, st['cwd'])
    if st['cwd'].startswith('../../s1_tight/'): cwd = os.path.join(KIT, 's1_tight', *st['cwd'].split('/')[3:])      # [kit fix 7] the proof tree's /home/claude/s1_tight is the kit's s1_tight (the old
    cwd = os.path.normpath(cwd)                                                                                    #   lookup left the kit and ran B8 in, and overwrote the output of, the proof tree)
    assert os.path.realpath(cwd).startswith(os.path.realpath(KIT) + os.sep), ('stage cwd outside the kit', cwd)
    logp = os.path.join(S4, st['log']); os.makedirs(os.path.dirname(logp), exist_ok=True)
    S.header(cmd, cwd)
    with open(logp, 'w', encoding='utf-8') as lf:
        p = popen(cmd, cwd, env, lf)
        try:
            while True:
                try: rc = p.wait(timeout=10); break
                except subprocess.TimeoutExpired: S.beat()
        except KeyboardInterrupt: interrupted(p)
    with open(S.log, 'a', encoding='utf-8') as g: g.write(f'   {st["stage"]}: exit {rc}  (log {rel(logp)})\n')
    if rc: print(f'   tail stage failed: {st["stage"]} (exit {rc}); log {rel(logp)}', flush=True); print(tail(logp, 12), flush=True)
    return rc == 0


def sha_file(p):
    try:
        with open(p, 'rb') as f: return hashlib.sha256(f.read()).hexdigest()
    except OSError: return None


def run_tail(S, final):
    st = FM['stages']; grp = lambda g: [x for x in st if x['group'] == g]
    kch = os.path.join(POUT, 't9_cone_kch.json')
    for x in grp('pre'):
        if final and x['stage'].startswith(('H1 exact points', 'H1 slab summary')): continue      # independent of T7/T8: done in round 1
        if not tail_run(S, x): return False
    if not final:                                                            # [kit fix 4+5] round 1 only makes n3_b1 for t7_nl: the T8 steps of pass A;
        for x in grp('A'):                                                   # n3boxes and t9 need T7's tip box (t7_nl), so they run in tail2 only
            if x['stage'].startswith('T8 '):
                if not tail_run(S, x, drop_env=('T9JSON',)): return False
            else: print('   tail: round 1 skips ' + x['stage'] + ' (run in tail2)', flush=True)
        return True
    if not os.path.exists(kch):                                              # bootstrap: the per-slab values t9 reads do not depend on T9JSON
        print('   tail: bootstrap pass without T9JSON', flush=True)
        for x in grp('A'):
            if not tail_run(S, x, drop_env=('T9JSON',)): return False
    prev = sha_file(kch); passes = 0
    while True:                                                              # pass A, pass B, ... until t9's output is a fixed point (the manifest has A, B)
        for x in grp('A') if passes == 0 else grp('B'):
            if not tail_run(S, x): return False
        passes += 1; cur = sha_file(kch)
        print(f'   tail: pass {passes}: t9_cone_kch.json {cur[:12]} ({"unchanged" if cur == prev else "changed"})', flush=True)
        if passes >= 2 and cur == prev: break
        if passes >= 5: print('   tail: no fixed point after 5 passes'); return False
        prev = cur
    for x in grp('post'):
        if not tail_run(S, x): return False
    return True


def st_tail1(S):
    return run_tail(S, final=False) and os.path.exists(os.path.join(T8SE, 'constants', 'n3_b1.json'))


def st_t7_nl(S):
    n3 = os.path.join('..', 'p2', 't8se', 'constants', 'n3_b1.json')
    with open(os.path.join(T8SE, 'constants', 'n3_b1.json'), encoding='utf-8') as f: ybox = str(json.load(f)['args']['Y'])
    snap = os.path.join(OUT, 'tail1_n3_b1.json'); shutil.copyfile(os.path.join(T8SE, 'constants', 'n3_b1.json'), snap)      # round-1 n3_b1: compared with round 2 by check.py
    if run_t7(S, 'run_v2se', ['--se', os.path.join('..', 'p2'), '--prev', 'run_v2', '--steps', 'merge', '--nlclosure', n3, '--ybox', ybox, '--publish']) is False: return False
    for d in ('slabs_b1',):                                                   # round 2 recomputes every T8 product from the republished T7
        shutil.rmtree(os.path.join(T8SE, 'constants', d), ignore_errors=True)
    for f in ('constants_b1.json',): 
        try: os.remove(os.path.join(T8SE, 'constants', f))
        except OSError: pass
    return True


def drop_args(o):
    """[kit fix] the n3_b1 files of the two rounds differ in their `args` block (argv, paths, hashes of the inputs); the numbers are what must agree"""
    if isinstance(o, dict): return {k: drop_args(v) for k, v in o.items() if k != 'args'}
    if isinstance(o, list): return [drop_args(v) for v in o]
    return o


def st_tail2(S):
    if not run_tail(S, final=True): return False
    snap = os.path.join(OUT, 'tail1_n3_b1.json'); cur = os.path.join(T8SE, 'constants', 'n3_b1.json')
    try:
        with open(snap, encoding='utf-8') as f: a = drop_args(json.load(f))
        with open(cur, encoding='utf-8') as f: b = drop_args(json.load(f))
        same = a == b
    except (OSError, ValueError): same = False
    with open(os.path.join(OUT, 'tail_rounds.json'), 'w') as f: json.dump(dict(n3_b1_unchanged=same, compared='all keys except args', sha256_equal=sha_file(snap) == sha_file(cur)), f)
    print(f'   n3_b1 of round 1 (T7 without the nonlinear closure) {"equals" if same else "DIFFERS FROM"} round 2 (republished T7), args dropped', flush=True)
    return True


# ---- the 21 stages of the v3.3 step-4 chain (reference: p2/out_v33/FINAL_manifest_v33.json of the proof tree; 366 paths, 39 scripts)
V33 = os.path.join(P2, 'v33'); OUTV = os.path.join(P2, 'out_v33')
V33_STAGES = [('h2f_ift', 0.2), ('h2f_balls', 8), ('t9_h2f', 0.5), ('h2f_camp', 0.2), ('h2t', 0.2), ('cone', 0.5), ('h4loc', 0.5), ('lemmaO_h1', 1.5), ('ratio', 0.2), ('notes', 0.2),
              ('b0_dprof', 0.5), ('h2u', 0.2), ('m9_balls', 16), ('h3_arb', 0.3), ('h2t_arb', 0.3), ('h2u_arb', 0.3), ('h2f_camp_arb', 0.3), ('lemmaO_arb', 1.5), ('h4far', 0.5),
              ('n1_chart', 0.5), ('f12_sliver', 1.5)]


def v33_stage(i):
    """stage i of final_tail_v33.py; earlier stages are reused from kit's out_v33/FINAL_manifest_v33.partial.json (written by the previous v33 stage of this kit)"""
    def fn(S):
        rc = S.run([PY, '-B', '-X', 'utf8', os.path.join(V33, 'final_tail_v33.py'), '--from', str(i), '--until', str(i)], V33)
        if rc: return False
        return os.path.exists(os.path.join(OUTV, 'FINAL_manifest_v33.json' if i == len(V33_STAGES) - 1 else 'FINAL_manifest_v33.partial.json'))
    return fn


V38 = os.path.join(P2, 'v38'); OUTV38 = os.path.join(P2, 'out_v38')
def st_v38(S):
    """[v4] the whole v38 driver (29 stages; carried stages are re-hashed, the rest rerun); a stale partial manifest is removed first"""
    for f in ('FINAL_manifest_v38.partial.json', 'FINAL_manifest_v38.json'):
        if os.path.exists(os.path.join(OUTV38, f)): os.remove(os.path.join(OUTV38, f))
    rc = S.run([PY, '-B', '-X', 'utf8', os.path.join(V38, 'final_tail_v38.py')], V38)
    return rc == 0 and os.path.exists(os.path.join(OUTV38, 'FINAL_manifest_v38.json'))


# name, function, estimate of wall minutes on the sandbox at 6 jobs (for the ETA only)
STAGES = [('inputs', st_inputs, 0.5), ('t1_seed', st_t1_seed, 0.2), ('t1_float', st_t1_float, 8), ('t1_polish', st_t1_polish, 40), ('t1_se', st_t1_se, 6),
          ('t6', st_t6, 12), ('rhoend', st_rhoend, 6), ('se_corr', st_se_corr, 1), ('t6c', st_t6c, 15), ('t7_a', st_t7_a, 75), ('t7_b', st_t7_b, 25),
          ('h1pts', st_h1pts, 20), ('h1dt', st_h1dt, 35), ('winslab', st_winslab, 3), ('tail1', st_tail1, 12), ('t7_nl', st_t7_nl, 2), ('tail2', st_tail2, 14),
          ('b8', st_b8, 0.3)] + [(f'v33_{i:02d}_{n}', v33_stage(i), e) for i, (n, e) in enumerate(V33_STAGES)] + [('v38', st_v38, 20)]
NAMES = [s[0] for s in STAGES]


def done_file(n): return os.path.join(STATE, n + '.done')
def is_done(n): return os.path.exists(done_file(n))


def write_env():
    info = dict(python=sys.version, platform=platform.platform(), processor=platform.processor(), cpu_count=os.cpu_count(), jobs=JOBS, started=time.strftime('%Y-%m-%d %H:%M:%S'))
    ok = True
    try:
        import flint, numpy, scipy, sympy, mpmath
        v = dict(flint=flint.__version__, numpy=numpy.__version__, scipy=scipy.__version__, sympy=sympy.__version__, mpmath=mpmath.__version__)
        info['versions'] = v; info['pinned'] = PINS
        for k, want in PINS.items():
            if v[k] != want: print(f'  note: {k} is {v[k]}, the reference run used {want}')
    except Exception as e:
        info['import_error'] = repr(e); ok = False
        print('  ERROR: cannot import the required packages (' + repr(e) + '). Run: python -m pip install -r requirements.txt')
    with open(os.path.join(OUT, 'env.json'), 'w', encoding='utf-8') as f: json.dump(info, f, indent=1)
    return ok


def main():
    global JOBS, CFG
    a = sys.argv[1:]
    if '-h' in a or '--help' in a: print(__doc__); return 0
    def opt(name, default=None):
        if name in a:
            i = a.index(name); return a[i + 1]
        return default
    JOBS = int(opt('--jobs', 6)); quick = False
    ensure_dirs()
    only = opt('--only'); frm = opt('--from'); force = (opt('--force') or '').split(',') if '--force' in a else []
    stages = [s for s in STAGES if not (quick and s[0] in QUICK_SKIP)]
    if only: stages = [s for s in stages if s[0] in only.split(',')]
    if frm: stages = stages[[s[0] for s in stages].index(frm):]
    if '--list' in a or '--status' in a:
        for n, _, e in STAGES: print(f'  {n:10s} {"done" if is_done(n) else "-"}   (about {e:g} min on the sandbox)')
        return 0
    if '--clean' in a:
        shutil.rmtree(OUT, ignore_errors=True); print('kit_out\\state and logs deleted (generated files under step4\\ stay; delete the kit folder for a pristine start).'); return 0
    acquire_lock()
    print(f'Kit 4B  |  {time.strftime("%Y-%m-%d %H:%M:%S")}  |  Python {sys.version.split()[0]}  |  {platform.platform()}  |  jobs {JOBS}' + ('  |  QUICK TEST MODE' if quick else ''), flush=True)
    if not write_env(): return 1
    t_start = time.time(); est_done = 0.0; act_done = 0.0
    remaining = [s for s in stages if (s[0] in force) or not is_done(s[0])]
    tot_est = sum(s[2] for s in remaining)
    for idx, (name, fn, est) in enumerate(stages, 1):
        if name not in force and is_done(name):
            print(f'[{now()}] [{idx}/{len(stages)}] {name}: already done, skipping', flush=True); continue
        ratio = (act_done / est_done) if est_done > 0 else 1.7
        rem_est = sum(s[2] for s in remaining if s[0] == name or NAMES.index(s[0]) > NAMES.index(name))
        print(f'[{now()}] [{idx}/{len(stages)}] {name}: start   (ETA for the whole remainder ~{hms(rem_est * 60 * min(max(ratio, 0.5), 4))})', flush=True)
        S = Stage(name); t0 = time.time()
        try: ok = fn(S)
        except SystemExit: raise
        except Exception as e:
            ok = False; print('   exception:', repr(e))
        dt = time.time() - t0
        if not ok:
            print(f'[{now()}] [{idx}/{len(stages)}] {name}: FAILED after {hms(dt)}.  Last log lines ({rel(S.log)}):\n{tail(S.log)}\nRunning  python run_all.py  again retries from this stage; if it fails again paste the log tail.', flush=True)
            return 1
        with open(done_file(name), 'w') as f: f.write(time.strftime('%Y-%m-%d %H:%M:%S') + f'  wall {dt:.0f}s\n')
        est_done += est; act_done += dt / 60
        print(f'[{now()}] [{idx}/{len(stages)}] {name}: ok ({hms(dt)})', flush=True)
    print(f'\n[{now()}] all stages finished in {hms(time.time() - t_start)}. Running check.py ...\n', flush=True)
    return subprocess.call([PY, 'check.py'], cwd=KIT, env=child_env())


if __name__ == '__main__':
    sys.exit(main())
