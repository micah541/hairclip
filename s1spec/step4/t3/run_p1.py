"""P1 production driver: resumable validated chain from the certified Gamma0 (t2/gamma0_final.json) to k_max >= target.

Run / resume (identical command; it continues from the last complete checkpoint; safe to kill at any time):
    mkdir -p /home/claude/s1spec/step4/t3/p1 && cd /home/claude/s1spec/step4/t3
    nohup python3 -u -I run_p1.py --target-k 12 >> p1/run.log 2>&1 &        (runs at nice 10 itself)
Status:  cat p1/status.json        Checkpoints: p1/ckpt/slab_NNNNN.{json,npz}   (slab 0 = Gamma0)
Extend later (e.g. to k = 20): rerun with --target-k 20 (the stop test is re-evaluated, the chain continues).

Per slab: graded rule dt = rho_g t (t from the curve's float time 0.001) or h = c/theta_max^4, whichever is smaller,
shortened until the float base's Chebyshev decay ratio r satisfies r^K >= 1e24; K = 24; block M_a = #{lam h <= 8};
M += 80 while the top 40 centre modes exceed 1e-26 (not in the graded phase), M -= 80 while the top 120 are below
1e-31; intra-slab absorption cap 1.5x the last measured eta/max(w); a failed slab is retried with h/2 and a 10x cap
(4 tries).  Slab 1 uses the slaved start for lam h > 26 and a preset cap.  The rigorous k_max enclosure is computed on
every slab once k > 0.8 target; the chain stops when its lower end reaches the target.
"""
import os, sys
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
import json, time, argparse, glob, math
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import slab as S
import run_proto as R

TF0 = 0.001            # float time label of Gamma0 (graded rule only; the certified clock starts at 0)


# ---------------------------------------------------------------- exact (de)serialisation
def a2j(x):
    """arb -> [man, exp, rad_man, rad_exp] exactly"""
    m, e = x.mid().man_exp(); r = x.rad(); rm, re_ = r.mid().man_exp()
    return [str(int(m)), int(e), str(int(rm)), int(re_)]


def j2a(v):
    mid = S.arb(int(v[0])) * S.arb(2) ** int(v[1])
    rad = S.arb(int(v[2])) * S.arb(2) ** int(v[3])
    return S.arb(mid, rad) if int(v[2]) != 0 else mid


def atomic_write(path, writer):
    tmp = path + '.tmp'
    writer(tmp)
    try:
        with open(tmp, 'rb+') as f: os.fsync(f.fileno())
    except OSError:
        pass
    os.replace(tmp, path)


def _dump(obj, indent=None):
    def w(p):
        with open(p, 'w') as f: json.dump(obj, f, indent=indent)
    return w


def save_ckpt(out, k, st, meta):
    base = os.path.join(out, 'ckpt', f'slab_{k:05d}')
    def wnpz(p):
        with open(p, 'wb') as f: np.savez(f, G=st.G, w=np.asarray(st.w, float), box=np.asarray(st.box, float))
    atomic_write(base + '.npz', wnpz)
    rec = dict(meta, slab=k, T=float(st.T), tau=float(st.tau), M=int(len(st.xbar)), d=int(st.G.shape[1]),
               xbar=[list(x.mid().man_exp()) for x in st.xbar_arb], L=a2j(st.L), t=a2j(st.t), X0=a2j(st.X0),
               Lint=(None if st.Lint is None else a2j(st.Lint)), complete=True)
    rec['xbar'] = [[str(int(m)), int(e)] for m, e in rec['xbar']]
    atomic_write(base + '.json', _dump(rec))


def load_ckpt(out):
    """latest complete checkpoint -> (k, state, meta) or None"""
    for f in sorted(glob.glob(os.path.join(out, 'ckpt', 'slab_*.json')), reverse=True):
        try:
            rec = json.load(open(f)); z = np.load(f[:-5] + '.npz')
            if not rec.get('complete'): continue
            xa = [S.arb(int(m)) * S.arb(2) ** int(e) for m, e in rec['xbar']]
            st = S.SlabState(np.array([float(x.mid()) for x in xa]), z['G'], z['w'], z['box'], rec['T'], rec['tau'],
                             None if rec['Lint'] is None else j2a(rec['Lint']), xa, j2a(rec['L']), j2a(rec['t']), j2a(rec['X0']))
            return rec['slab'], st, rec
        except Exception as ex:                         # torn/partial file: fall back to the previous checkpoint
            print(f'  skipping unreadable checkpoint {f}: {ex}', flush=True)
    return None


def gamma0_state(path, M):
    g = json.load(open(path))
    bfull = np.array([float.fromhex(x) for x in g['b_hex']]); j = int(g['j'])
    bj = S.arb(int(g['b_j_center_exact']['man'])) * S.arb(2) ** int(g['b_j_center_exact']['exp'])
    hw = float(g['b_j_halfwidth']); L0 = float.fromhex(g['L0_hex']); x0 = float.fromhex(g['X0_hex'][0])
    # the closing centre has up to 192 bits: use its 128-bit midpoint and add the rounding to the b_j interval
    hw += float(bj.rad().upper().mid()) * 2 + 1e-300; bj = S.arb(bj.mid())
    xa = [S.arb(float(v)) for v in bfull[:M]]; xa[j - 1] = bj
    b = bfull[:M].copy(); b[j - 1] = float(bj.mid()); fm = S.FM(M)
    h = 0.5 * TF0 / (L0 / 2 / np.pi) ** 4; Ma = int(np.sum(fm.lam * h <= 8.0))
    st = S.initial_state(b, Ma, 4 * M, 1e-30, fm)
    st.w[j - 1] += hw * (1 + 1e-12); st.box[Ma:M] += np.abs(st.G[Ma:, j - 1]) * hw * (1 + 1e-12)
    m2 = min(len(bfull), 4 * M); st.box[M:m2] += np.abs(bfull[M:m2]) * (1 + 1e-12)
    st.xbar_arb = xa; st.L = S.arb(L0); st.t = S.arb(0); st.X0 = S.arb(x0)
    return st


# ---------------------------------------------------------------- monitors (float; for the error budget)
def norms(st):
    """|e|_j bounds (j = 0..4) of the set: hull of the parallelepiped + box + tail"""
    M = len(st.xbar); Nc = len(st.box); n = np.arange(1, Nc + 1, dtype=float)
    col = np.abs(st.G) @ np.asarray(st.w, float)                    # sum_i |G_ni| w_i, n <= M
    out = []
    for j in range(5):
        tail = st.T * Nc ** (j + 1 - S.QT) / (S.QT - j - 1)
        out.append(float((n[:M] ** j * col).sum() + (n ** j * st.box).sum() + tail))
    return out


def budget(st, fm):
    e = norms(st); L = float(st.L.mid()); q = 2 * np.pi / L
    th = float(fm.thmax(st.xbar[None])[0]); kap = q * th
    c2 = q * e[1] + kap * e[0] + kap * float(st.L.rad()) / L
    c5 = q ** 4 * e[4]
    c5full = sum(math.comb(4, j) * kap ** (4 - j) * q ** j * e[j] for j in range(5))
    return dict(e=e, C2=c2, C5_lead=c5, C5_crude=c5full, kmax_float=kap, theta=th)


# ---------------------------------------------------------------- one slab with retries
def run_slab(st, fm, h, K, k, meta, target, args):
    t_abs = TF0 + float(st.t.mid()); Lc = float(st.L.mid())
    def plan(uh, h=h):
        hn = hrule(uh, t_abs + h * (Lc / 2 / np.pi) ** 4, Lc, args)
        return hn, int(np.sum(fm.lam * hn <= 8.0))
    kf = 2 * np.pi / Lc * float(fm.thmax(st.xbar[None])[0])
    kng = 2048 if kf > 0.8 * target else 0
    cap = meta.get('etacap', args.cap1 if k == 0 else 1e-2)
    for attempt in range(4):
        S.ETACAP = cap; S.validate_slab._depth = 0
        try:
            c0 = time.process_time()
            st_new, info = S.validate_slab(st, fm, h, K=K, Ma=st.G.shape[1], verbose=False, plan_next=plan,
                                           kmax_ng=kng, lamcut=26.0 if k == 0 else None)
            info['cpu'] = time.process_time() - c0; info['attempt'] = attempt; info['h_used'] = h
            return st_new, info
        except (RuntimeError, ValueError, np.linalg.LinAlgError) as ex:
            print(f'  slab {k + 1} attempt {attempt} failed ({ex}); retry with h/2, cap x10', flush=True)
            h *= 0.5; cap *= 10
    raise RuntimeError('slab failed 4 times')


def hrule(x, t_abs, Lf, args):
    th = R.thmax_any(x)
    return min(args.rho_g * t_abs / (Lf / 2 / np.pi) ** 4, args.c / th ** 4)


def hcheck(st, fm, h, K):
    rreq = 10 ** (24.0 / K)
    for it in range(4):
        Ms = int(np.sum(fm.lam * h <= 8.0 * K * K))
        U = S.base_float(fm, st.xbar, h, K, np.tile(st.xbar, (K + 1, 1)), min(max(int(np.sum(fm.lam * h <= 8)), 60), 140), Ms)
        cm = np.abs(S.vinv_f(K) @ U).max(axis=1)
        r = (cm[1] / cm[6]) ** 0.2
        if r >= 1.1 * rreq: return h, r
        h = h * r / (1.25 * rreq)
    return h, r


def mrule(x, Mmax, graded):
    """+80 while the top 40 centre modes exceed 1e-26 (truncation residual ~1e12 |x_M| <~ 1e-14); not during the
    graded phase, where the n > M content of Gamma0 (<= 1e-23) only enters the box and decays in the first
    sub-interval; -80 while the top 120 modes are below 1e-31"""
    M = len(x); a = np.abs(x)
    if a[M - 40:].max() > (1e-20 if graded else 1e-26): return min(Mmax, M + max(80, 40 * int(round(0.15 * M / 40))))  # Nov 10: step max(80, 0.15 M) (T3 fix a)
    if M > 360 and a[M - 120:].max() < 1e-31: return M - 80
    return M


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--target-k', type=float, default=12.0)
    ap.add_argument('--out', default=os.path.join(HERE, 'p1'))
    ap.add_argument('--gamma0', default='/home/claude/s1spec/step4/t2/gamma0_final.json')
    ap.add_argument('--M0', type=int, default=1000)
    ap.add_argument('--Mmax', type=int, default=2400)
    ap.add_argument('--K', type=int, default=24)
    ap.add_argument('--rho-g', type=float, default=0.5)
    ap.add_argument('--c', type=float, default=2.0)
    ap.add_argument('--max-slabs', type=int, default=10 ** 6, help='stop after this many slabs in total (testing)')
    ap.add_argument('--nice', type=int, default=10)
    ap.add_argument('--cap1', type=float, default=30.0, help='slab-1 absorption cap relative to max w (b_2 halfwidth)')
    args = ap.parse_args()
    if args.nice: os.nice(args.nice)
    os.makedirs(os.path.join(args.out, 'ckpt'), exist_ok=True)
    got = load_ckpt(args.out)
    if got is None:
        st = gamma0_state(args.gamma0, args.M0); k = 0
        meta = dict(gamma0=args.gamma0, h_next=None, cpu_total=0.0, history=[])
        save_ckpt(args.out, 0, st, meta); print(f'started from {args.gamma0}: M={args.M0}, block {st.G.shape[1]}', flush=True)
    else:
        k, st, meta = got; print(f'resumed at slab {k}: t = {st.t.str(12)}, M = {len(st.xbar)}, block {st.G.shape[1]}', flush=True)
    fm = S.FM(len(st.xbar))
    while True:
        b = budget(st, fm)
        last = meta['history'][-1] if meta['history'] else {}
        klo = (last.get('kmax') or [0, 0])[0]
        done = klo >= args.target_k
        status = dict(slab=k, done=done, target_k=args.target_k, t=st.t.str(20), t_label=TF0 + float(st.t.mid()),
                      L=st.L.str(20), X0=st.X0.str(20), kmax_float=b['kmax_float'], kmax_enclosure=last.get('kmax'),
                      theta_max=b['theta'], M=len(st.xbar), block=int(st.G.shape[1]), h_next=meta.get('h_next'),
                      rho_max=float(np.max(st.w)), e_norms=b['e'], C2=b['C2'], C5_lead=b['C5_lead'], C5_crude=b['C5_crude'],
                      cpu_total_h=meta['cpu_total'] / 3600, last_slab=last, updated=time.strftime('%Y-%m-%d %H:%M:%S'))
        atomic_write(os.path.join(args.out, 'status.json'), _dump(status, 1))
        if done:
            print(f'target k = {args.target_k} reached at slab {k}: k_max in {klo}..; C2 ~ {b["C2"]:.2e}, C5 ~ {b["C5_lead"]:.2e}', flush=True)
            break
        if k >= args.max_slabs:
            print(f'stopping after slab {k} (--max-slabs)', flush=True); break
        t_abs = TF0 + float(st.t.mid()); Lc = float(st.L.mid())
        h = meta.get('h_next') or hrule(st.xbar, t_abs, Lc, args)
        try:
            h2, rdec = hcheck(st, fm, h, args.K)
        except Exception as ex:
            h2, rdec = 0.5 * h, float('nan'); print(f'  hcheck failed ({ex}); h/2', flush=True)
        print(f'--- slab {k + 1}: h {h:.3e} -> {h2:.3e} (decay ratio {rdec:.1f}), t = {t_abs:.6g}, theta_max {b["theta"]:.2f}, '
              f'M {len(st.xbar)}, block {st.G.shape[1]}, cap {meta.get("etacap", args.cap1 if k == 0 else 1e-2):.2g}', flush=True)
        st_new, info = run_slab(st, fm, h2, args.K, k, meta, args.target_k, args)
        print(R.summary(info, Lc)); print(R.qline(info, st_new), flush=True)
        q = info.get('quad') or {}
        rec = dict(slab=k + 1, h=info['h_used'], cpu=info['cpu'], attempt=info['attempt'], M=info['M'], block=[info['d'], info['a1']],
                   resid_max=info['resid_max'], w_out=info['wnorm_out'], eta_intra=info['eta_intra'], hull1=info['hull1_out'],
                   box1=info['box1_out'], box3=info['box_out_w3'], T=info['T_out'], theta_end=info['theta_end'],
                   t=q.get('t'), kmax=q.get('kmax'))
        meta['history'].append(rec); meta['cpu_total'] += info['cpu']
        meta['etacap'] = min(1.0, max(1e-2, 1.5 * info['eta_intra'] / max(info['w_out'], 1e-300)))   # w_out = max radius
        meta['h_next'] = info['h_next']
        st = st_new
        graded = args.rho_g * (TF0 + float(st.t.mid())) / (float(st.L.mid()) / 2 / np.pi) ** 4 < args.c / R.thmax_any(st.xbar) ** 4
        Mn = mrule(st.xbar, args.Mmax, graded)
        if Mn != len(st.xbar):
            st = S.resize_M(st, fm, Mn); fm = S.FM(Mn); print(f'  M resized -> {Mn}', flush=True)
        k += 1
        save_ckpt(args.out, k, st, meta)


if __name__ == '__main__':
    main()
