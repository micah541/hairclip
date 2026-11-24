"""Informational (not a hypothesis of step 3): the hand-off curve self-intersects.  On each of NT t_c-balls covering slab SLAB, the TRUE arm
gamma = Gamma~ + e N crosses the real axis at an interior point of the arm (base arclength x >= XMIN, away from the tip and the far end);
the full curve is the arm and its mirror image (reflection symmetry), so it passes twice through that point.
Positions from the T6 engine (p2/t9h/basegeo.py, as t9_lemmaO.py), correlated in time: Gamma(s, t) = Gamma(s, 0) + t Gamma_t(s, xi),
|t| <= 1 (t_c), with Gamma_t enclosed over the whole slab; on each element the exact points s_j = -1 + j/NG (j = 0..2 NG).
Sign certificate at a point: Im Gamma~ > |e| (resp. < -|e|), |e| <= c_e Ybar + sqrt(F0) (y^3/84)^(1/2), y = max(x - X0, 28) (T7 Lemma 4b
constant 1.6352e4 and LEMMA_F C4 with the amp's sqrt(F0), as t9_lemmaO.py), so Im gamma has that sign there.  A crossing on a t_c-ball:
two grid points of one element, consecutive among the certified ones, with opposite certified signs (the true arm is continuous in s).
For each ball the first crossing along the arm is reported (element, s, x hull, Re mid, |e|, the certified Im bounds).
usage: python3 -I t9_axis_crossing.py AMP_JSON [SLAB]   (env AX_NT: t_c-balls, default 16; AX_NG: s-grid spacing 1/AX_NG on [-1, 1], default 32)
output: p2/out/t9_axis_crossing.json"""
import sys, os, json, math
sys.path.insert(0, '/home/claude/s1spec/step4/p2/t9h')
from flint import arb, acb, arb_poly, ctx
import basegeo
ctx.prec = 128
SLAB = int(sys.argv[2]) if len(sys.argv) > 2 else 180; NT = int(os.environ.get('AX_NT', '16')); NG = int(os.environ.get('AX_NG', '32'))
amp = json.load(open(sys.argv[1])); sF0 = arb(repr(amp['handoff']['at_Ybar']['sqrtF'][0])); YB = arb('1e-12'); X0 = arb('2.2490776515')
XMIN = 150.0
U = lambda v: math.nextafter(float(v.upper()), math.inf); Lo = lambda v: math.nextafter(float(v.lower()), -math.inf)
def at_dt(g, i, s):                                              # t_c-derivative of the position at fixed s (as t9_lemmaO.py)
    el = g.els[i]; T_ = g.T; sh = arb_poly([arb(T_.mid()), arb(1)]); dl = arb(0, T_.rad())
    tv = lambda tp: tp(sh)(dl); tvd = lambda tp: tp.derivative()(sh)(dl)
    Gr, Gi = el['G'].at(s); Pr, Pi = el['Pg'].at(s)
    return acb(tvd(Pr), tvd(Pi)) * acb(tv(Gr), tv(Gi)) + acb(tv(Pr), tv(Pi)) * acb(tvd(Gr), tvd(Gi))
def ebox(xhi):
    y = max(arb(xhi) - X0, arb(28)); return arb('1.6352e4') * YB + sF0 * (y ** 3 / 84).sqrt()
g0 = basegeo.Geo(SLAB, 4.0e5, arb(0)); gB = basegeo.Geo(SLAB, 4.0e5, arb(0, 1))
SP = [arb(-1) + arb(j) / NG for j in range(2 * NG + 1)]          # exact points
pts = {}                                                          # (element, j) -> (Gamma(s, 0), Gamma_t)
for i in range(len(g0.els)):
    if float(gB.xe[i].upper() + gB.hh[i].upper()) < XMIN: continue
    for j, s in enumerate(SP): pts[(i, j)] = (g0.at(i, s)['Gam'], at_dt(gB, i, s))
TB = [arb(-1) + arb(2 * k + 1) / NT + arb(0, 1) / NT for k in range(NT)]       # NT balls of radius 1/NT covering t_c in [-1, 1]
XEND = gB.xe[-1].lower()                                         # the far end of the arm lies beyond the start of its last element
cross = []
for k, T_ in enumerate(TB):
    found = None
    for i in range(len(g0.els)):
        if found: break
        if (i, 0) not in pts: continue
        xe, hh = g0.els[i]['x'](T_), g0.els[i]['h'](T_)
        last = None
        for j, s in enumerate(SP):
            G0, Gt = pts[(i, j)]; im = G0.imag + T_ * Gt.imag; x = xe + hh * (1 + s) / 2
            if float(x.lower()) < XMIN: continue
            e = ebox(x.upper()); sg = 1 if im > e else (-1 if im < -e else 0)
            if sg == 0: continue
            if last is not None and last[0] == -sg:
                ja, sa, xa, ima, ea = last[1:]; assert x.upper() < XEND   # interior point of the arm (not the far end)
                re = G0.real + T_ * Gt.real
                found = dict(tc=[Lo(T_), U(T_)], element=i, s=[float(SP[ja].mid()), float(s.mid())], x=[Lo(xa), U(x)], Re=float(re.mid()),
                             e=max(U(ea), U(e)), Im=[[Lo(ima), U(ima)], [Lo(im), U(im)]])
                break
            last = (sg, j, s, x, im, e)
    cross.append(found)
assert TB[0].lower() <= -1 and TB[-1].upper() >= 1
out = dict(args=dict(argv=sys.argv[1:], slab=SLAB, tballs=NT, grid=f'-1 + j/{NG}', x_min=XMIN, err='T7 4b + C4 core (as t9_lemmaO.py)',
                     note='Im of the base half-curve changes sign beyond |e| (true curve), so the arm crosses the axis; the full curve (arm + mirror) self-intersects there'),
           all_t=all(c is not None for c in cross), crossings=cross)
json.dump(out, open('/home/claude/s1spec/step4/p2/out/t9_axis_crossing.json', 'w'), indent=1)
print(f"axis crossing on {sum(c is not None for c in cross)}/{NT} t_c-balls: " + ', '.join(
    f"el {c['element']} s {c['s']} x [{c['x'][0]:.0f}, {c['x'][1]:.0f}] |e| <= {c['e']:.3f}" for c in cross if c)[:600])
assert out['all_t']
