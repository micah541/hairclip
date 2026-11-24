"""v3.3 fix B0, stage 1: certified bounds on the profile difference between the two deviation variables of v3.3,
    v3.2 shell / verifier variable  w' = f - alpha psi3^'(x) - alpha^2 h2^'(x)        (x = y - SIG_TIP, step-4 win.py 'shell' mode)
    v3.3 energy variable            f - U(alpha, y),  U(alpha, y) = 1.3056 alpha y^p,  p = 1.0246
    d := (f - U) - w' = alpha g(y) + alpha^2 h(y),   g = psi3^'(x) - 1.3056 y^p,   h = h2^'(x).
psi3^, h2^ are the certified C^4 profile functions of step 3 (step3nl/agent_int/v3/c4fun.py, the module step-4 win.py and ift.py use),
evaluated with their own arb jets; SIG_TIP from p2/t9h/s1chart.py (as win.py).  All arithmetic in arb (prec 128), outward rounding.
Per cell Y = [y0, y1] (width 1/64) on y in [14, 60]:
  D = y d/dy (= d_s, s = ln y);  D^k h = sum_m S(k,m) y^m h^(m+1)(x)  (naive on the ball);
  D^k g = sum_m S(k,m) y^m psi^(m+1)(x) - 1.3056 p^k y^p, bounded by the first-order form  D^k g(y_c) + [-r, r] sup_Y |D^(k+1) g / y|
  (d/dy D^k g = D^(k+1) g / y; this avoids the cancellation psi3^' ~ 1.3056 y^p).
  x-derivatives (verifier convention): y^m d^m/dy^m = D (D - 1) ... (D - m + 1), so |d^m g| <= y^-m sum_k |s(m,k)| |D^k g|.
Tables (k = 0..3, m = 0..3), for |alpha| <= a in {a_t1, alpha_1}:
  E_k[I] = int_I (28/y) (D^k d)^2 ds   (ds = dy/y), I = [14, 60], [14, 28], [28, 60]: sum over cells meeting I of sup_cell (|a| G_k + a^2 H_k)^2 * 28 (1/y0 - 1/y1);
  sup_k[I] = max over cells of |a| G_k + a^2 H_k;
  W_m (ramp, verifier weight) = int_{14}^{28} W(x) (d^m d / dy^m)^2 dy,  W = (1 + x^2)^-2,  x = y + X0  (X0: the verifier's 2.2490776515 and -SIG_TIP, union).
  a_t1 = max |alpha| over the 32 widened slab-180 balls p2/out_v33/h2f/ift_m.json (true curve at t_1, h2f_ift.py);
  alpha_1 = 1.5e-4 (the bootstrap bound |alpha| <= alpha_1 of v3.3, Lemma dU; a hypothesis of step 3, not of step 4).
Also coefficient-only tables (sqrt E of g and h separately) for the record.  The chart allowance 2e-7 of win.py (|theta1 - Theta| <= B(14),
level 0) is recorded as a convention value (label 'step-4 convention'), not added here.
Writes p2/out_v33/b0_dprof.json.   usage: python3 -I -B b0_dprof.py"""
import sys, os, json, math, time, hashlib
sys.path.insert(0, '/home/claude/s1spec/step4/p2/t9h'); sys.path.insert(0, '/home/claude/s1spec/step3nl/agent_int/v3'); sys.path.insert(0, '/home/claude/s1spec/step3nl')
from flint import arb, ctx
ctx.prec = 128
import s1chart
import c4fun as C4
P2 = '/home/claude/s1spec/step4/p2'; OV = f'{P2}/out_v33'
sha = lambda p: hashlib.sha256(open(p, 'rb').read()).hexdigest()
U = lambda v: math.nextafter(float(abs(v).upper()), math.inf)
STIP = s1chart.SIG_TIP
P = arb('1.0246'); C13 = arb('1.3056')
X0V = arb('2.2490776515'); X0 = X0V.union(-STIP)
STIR = {0: {0: 1}, 1: {1: 1}, 2: {1: 1, 2: 1}, 3: {1: 1, 2: 3, 3: 1}, 4: {1: 1, 2: 7, 3: 6, 4: 1}}
S1ST = {0: {0: 1}, 1: {1: 1}, 2: {1: -1, 2: 1}, 3: {1: 2, 2: -3, 3: 1}}          # y^m d^m = sum_k s(m,k) D^k (signed Stirling, 1st kind)
YLO, YHI, NPER = 14, 60, 64
IFTW = [f'{OV}/h2f/ift_{m}.json' for m in range(32)]
C4F = ['/home/claude/s1spec/step3nl/agent_int/v3/c4fun.py', '/home/claude/s1spec/step3nl/psi3fun.py', '/home/claude/s1spec/step3nl/h2fun.py',
       f'{P2}/t9h/s1chart.py']


def Dk(J, y, k, extra=None):
    """sum_m S(k,m) y^m J[m+1]"""
    return sum((arb(c) * y ** m * J[m + 1] for m, c in STIR[k].items()), arb(0))


def cell(i):
    y0 = arb(YLO) + arb(i) / NPER; y1 = arb(YLO) + arb(i + 1) / NPER
    yc = (y0 + y1) / 2; r = arb(1) / (2 * NPER); Y = yc + arb(0, 1) * r      # Y contains [y0, y1]
    X = Y - STIP; Xc = yc - STIP
    Jp = C4.PSI.jets_any(X, float(X.mid()), 6); Jh = C4.H2.jets_any(X, float(X.mid()), 6)
    Jpc = C4.PSI.jets_any(Xc, float(Xc.mid()), 5)
    G = []; H = []
    for k in range(4):
        gc = Dk(Jpc, yc, k) - C13 * P ** k * yc ** P                              # point value (ball of STIP width only)
        g1 = (Dk(Jp, Y, k + 1) - C13 * P ** (k + 1) * Y ** P) / Y                 # d/dy D^k g on the cell
        G.append(arb(U(gc + arb(0, 1) * r * arb(U(g1)))))
        H.append(arb(U(Dk(Jh, Y, k))))
    Gx = [arb(U(sum((abs(arb(c)) * G[k] for k, c in S1ST[m].items()), arb(0)) / y0 ** m)) for m in range(4)]
    Hx = [arb(U(sum((abs(arb(c)) * H[k] for k, c in S1ST[m].items()), arb(0)) / y0 ** m)) for m in range(4)]
    Wb = (1 + (Y + X0) ** 2) ** -2
    return dict(y0=y0, y1=y1, G=G, H=H, Gx=Gx, Hx=Hx, W=arb(U(Wb)), om=arb(28) * (1 / y0 - 1 / y1), wd=y1 - y0)


def tables(cells, a):
    a = arb(a); out = {}
    for nm, lo, hi in (('14_60', 14, 60), ('14_28', 14, 28), ('28_60', 28, 60)):
        sel = [c for c in cells if c['y1'].upper() > lo and c['y0'].lower() < hi]
        E = [sum(((a * c['G'][k] + a * a * c['H'][k]) ** 2 * c['om'] for c in sel), arb(0)) for k in range(4)]
        S = [max((a * c['G'][k] + a * a * c['H'][k] for c in sel), key=lambda v: float(v.upper())) for k in range(4)]
        out['E_' + nm] = [U(e) for e in E]; out['sup_' + nm] = [U(s) for s in S]
    sel = [c for c in cells if c['y0'].lower() < 28]
    out['W_ramp_14_28'] = [U(sum((c['W'] * (a * c['Gx'][m] + a * a * c['Hx'][m]) ** 2 * c['wd'] for c in sel), arb(0))) for m in range(4)]
    return out


if __name__ == '__main__':
    t0 = time.time()
    n = (YHI - YLO) * NPER
    cells = [cell(i) for i in range(n)]
    assert cells[0]['y0'] == arb(YLO) and (cells[-1]['y1'] - arb(YHI)) == 0
    al_t1 = max(max(abs(json.load(open(p))['alpha'][0]), abs(json.load(open(p))['alpha'][1])) for p in IFTW)
    al_t1 = math.nextafter(al_t1, math.inf)
    coef = dict(sqrtE_g_14_60=[U(sum((c['G'][k] ** 2 * c['om'] for c in cells), arb(0)).sqrt()) for k in range(4)],
                sqrtE_h_14_60=[U(sum((c['H'][k] ** 2 * c['om'] for c in cells), arb(0)).sqrt()) for k in range(4)],
                sup_g_14_60=[max(U(c['G'][k]) for c in cells) for k in range(4)], sup_h_14_60=[max(U(c['H'][k]) for c in cells) for k in range(4)])
    out = dict(item='B0: profile difference d = (f - U(alpha,y)) - w\' = alpha (psi3^\'(x) - 1.3056 y^p) + alpha^2 h2^\'(x) on y in [14, 60]', label='certified (arb, outward)',
               definitions=dict(D='y d/dy = d_s', energy='int (28/y) (D^k d)^2 ds', W_ramp='int_14^28 (1 + (y + X0)^2)^-2 (d_y^m d)^2 dy', levels_k=[0, 1, 2, 3], orders_m=[0, 1, 2, 3]),
               alpha_t1=al_t1, alpha_t1_source='max |alpha| over p2/out_v33/h2f/ift_m.json (widened, true curve at t_1)',
               alpha_1='1.5e-4', alpha_1_label='hypothesis of step 3 (v3.3 bootstrap |alpha| <= alpha_1, Lemma dU)',
               at_alpha_t1=tables(cells, al_t1), at_alpha_1=tables(cells, arb('1.5e-4')), coefficients=coef,
               chart_allowance_level0=dict(value=2e-7, label='step-4 convention (win.py shell mode Bs: |theta1 - Theta| <= B(14)); not included above'),
               cells=n, cell_width=f'1/{NPER}',
               args=dict(argv=sys.argv[1:], env={}, inputs={p: sha(p) for p in [os.path.abspath(__file__)] + C4F + IFTW}))
    tmp = f'{OV}/b0_dprof.json.tmp'; json.dump(out, open(tmp, 'w'), indent=1); os.replace(tmp, f'{OV}/b0_dprof.json')
    for key in ('at_alpha_t1', 'at_alpha_1'):
        T = out[key]; print(key, 'E_14_60', ['%.4e' % v for v in T['E_14_60']], 'E_28_60', ['%.3e' % v for v in T['E_28_60']],
                             'sup_14_28', ['%.3e' % v for v in T['sup_14_28']], 'sup_28_60', ['%.3e' % v for v in T['sup_28_60']], 'W_ramp', ['%.3e' % v for v in T['W_ramp_14_28']])
    print('coef sqrtE g', ['%.3e' % v for v in coef['sqrtE_g_14_60']], 'sqrtE h', ['%.3e' % v for v in coef['sqrtE_h_14_60']], f'alpha_t1 {al_t1:.6e} ({time.time() - t0:.1f}s)')
