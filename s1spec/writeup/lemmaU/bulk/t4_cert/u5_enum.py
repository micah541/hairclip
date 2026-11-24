#!/usr/bin/env python3
"""Task 4, U5: cross-check of P2's direct ramp classes (as rerun in u2_K.py) with the class enumerator of (A)
(t2_A/rigor/ramp_rep.py: level_terms, fold; its own tables WX.build_tables() and its own A-coefficient bounds arlib.at_eval), on [28, 56], weight Theta*phi.

  PYTHONHASHSEED=0 nice -n 19 python3 -B u5_enum.py [--u2 out/u2_K.json]  -> out/u5_enum.json

Imported read-only (no file of task 2 is written: -B, its out/ dir exists).  ramp_rep.SUPLOW is set IN MEMORY to Lemma U's folding sups (0.11, 0.534).
Two runs of the A-enumerator, Lemma U's ramp box (0.11, 0.52), U = 1/80:
  (a) DIRECT: the enumerator's code path with every pair taken as it stands (weight (pair, d = 0, 1); no integration by parts), i.e. Lemma U's route.
      Every class (pair or folded jet class) must be present in P2's list, and its weight sup_ramp sum |c| |rem_m(Theta phi)| must agree with P2's within 5%
      (A's entries are summed before the sup, P2 folds per unfolded key after the sup, so A <= P2 up to b&b tolerance is expected).
      Pair columns: K_jk (k != j) and the pair part of K_jj compared with u2's K / K_pair_only plus the folded entries.
  (b) NATIVE (route J of (A)): its integration-by-parts route for the pairs (a, j+1); the extra classes it creates are listed as explained
      ('created by the IBP of (A)'s route, not used by Lemma U, whose ramp terms carry no integration by parts).
Weights use rampphi.ramp_q (same formulas as u2) and ulib.bb_sup (certified upper endpoints over a cover)."""
import sys, os, json, math, time, argparse, inspect, textwrap
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import t4lib as L
import rampphi as RP
sys.path.insert(0, L.T2 + '/rigor')
import ramp_rep as RR                      # imports arlib (A's machinery; chdir restored)
from flint import arb, ctx
from t4lib import Q_, fup
import ulib as U
ctx.prec = 128
ap = argparse.ArgumentParser(); ap.add_argument('--u2', default=os.path.join(L.OUTD, 'u2_K.json')); ap.add_argument('--tol', type=float, default=0.003)
ap.add_argument('--maxeval', type=int, default=150000); args = ap.parse_args()
T0 = time.time()
U2 = json.load(open(args.u2))
RP.G_LO, RP.G_HI = (Q_(x) for x in U2['g_range'])
M0R, M1R, UR = Q_(11, 100), Q_(52, 100), Q_(1, 80)
RR.SUPLOW = {0: RR.A(M0R), 1: RR.A(Q_(U2['global_box']['M1z']))}
RR.PM3_DIRECT = False
TA = RR.WX.build_tables(); AtA = RR.L.at_eval(M0R, M1R, UR); AbA = RR.L.ab_honest(AtA)

# DIRECT variant of level_terms: same code, every pair weight ((a, b), 0, 1)
src = textwrap.dedent(inspect.getsource(RR.level_terms))
old = """        if b <= j: wts = [((a, b), 0, Q_(1))]"""
assert src.count(old) == 1
src = src.replace(old, """        if True: wts = [((a, b), 0, Q_(1))]""").replace('def level_terms(', 'def level_terms_direct(')
ns = dict(RR.__dict__); exec(src, ns); level_terms_direct = ns['level_terms_direct']

def bb(f):
    ub, lb, nev, left = U.bb_sup(f, RP.ramp_boxes(), tol=args.tol, maxeval=args.maxeval); return dict(ub=ub, lb=lb, nev=nev)
def weight(entries):
    """sup over the ramp of sum |c| |rem_m| (q = 0 entries only)"""
    coef = {}
    for c, (m, q) in entries:
        assert q == 0
        coef[m] = coef.get(m, arb(0)) + arb(abs(c).upper())
    return bb(RP.mkramp(coef, None))

res = dict(direct={}, native={}, compare={})
maxrel = 0.0; unmatched = []
for j in (0, 1, 2):
    n = j + 2
    jf, cls = level_terms_direct(j, TA, AtA, AbA, 'J')
    # pair columns, Lemma U weights w_k (1 for (k,k), 1/2 per index of a mixed pair) -- folded jet classes landing on pairs added like P2's augment
    col = {}
    for (tag, pr), lst in jf.items():
        a, b = pr
        for k in set(pr):
            w = Q_(1) if a == b else Q_(1, 2)
            col.setdefault(k, []).extend([(c*RR.A(w), kk) for c, kk in lst])
    jetcls = {}
    for (tag, I, ft), lst in cls.items():
        if len(I) == 2:
            a, b = I
            for k in set(I):
                w = Q_(1) if a == b else Q_(1, 2)
                col.setdefault(k, []).extend([(e['c']*RR.A(w), e['k']) for e in lst])
        jetcls[str(tuple(I))] = [(e['c'], e['k']) for e in lst]
    rowA = {}
    for k, ent in sorted(col.items()):
        wv = weight(ent); rowA[str(k)] = wv['ub']
    # P2 side: K (k != j) / K_pair_only (k = j) plus folded classes landing on (k, .)
    rowP = {}
    for k in rowA:
        base = U2['K'][str(j)].get(k, 0.0) if int(k) != j else U2['K_pair_only'][str(j)][k]
        fold = sum((fc['weight']*(1 if fc['to'][0] == fc['to'][1] else 0.5)) for fc in U2['folded_classes'] if fc['level'] == j and int(k) in fc['to'])
        rowP[k] = base + fold
    cmpk = {k: dict(A_enum=rowA[k], P2=rowP.get(k), rel=(rowA[k]/rowP[k] - 1) if rowP.get(k) else None) for k in rowA}
    for k, v in cmpk.items():
        if v['P2'] is None: unmatched.append(f'level {j} column {k}')
        else: maxrel = max(maxrel, abs(v['rel']))
    cub = {}
    for key, ent in jetcls.items():
        if len(eval(key)) == 2: continue
        wv = weight(ent)['ub']; p2 = U2['cubic_remaining'].get(f'{j}|{key}')
        cub[key] = dict(A_enum=wv, P2=p2, rel=(wv/p2 - 1) if p2 else None)
        if p2 is None: unmatched.append(f'level {j} cubic {key}')
        else: maxrel = max(maxrel, abs(cub[key]['rel']))
    p2cub = [k.split('|')[1] for k in U2['cubic_remaining'] if k.startswith(f'{j}|')]
    for k in p2cub:
        if k not in cub: unmatched.append(f'level {j} cubic {k} (in P2, not in A-enumerator)')
    res['direct'][str(j)] = dict(columns=cmpk, cubic=cub, jet_classes=sorted(jetcls))
    # native route J of (A)
    jfN, clsN = RR.level_terms(j, TA, AtA, AbA, 'J')
    natkeys = sorted({(str(tuple(I)), ft) for (tag, I, ft) in clsN}, key=str)
    dirkeys = sorted(jetcls)
    extra = [k + ('' if ft is None else ' ' + ft) for k, ft in natkeys if k not in dirkeys]
    natkeys = [k + ('' if ft is None else ' ' + ft) for k, ft in natkeys]
    res['native'][str(j)] = dict(classes=natkeys, not_in_direct=extra, jetfree_pairs=sorted(str(pr) for (tag, pr) in jfN),
                                 explanation='classes created by the integration by parts of pairs (a, j+1) in (A)\'s route J (derivatives of A_i^(r)); '
                                             'Lemma U keeps these pairs as data (columns k = j+1 of K) and creates no such class')
    print(f'level {j}: columns ' + ', '.join(f'k={k}: A {v["A_enum"]:.6g} P2 {v["P2"]:.6g} ({v["rel"]:+.2e})' for k, v in cmpk.items()), flush=True)
    for k, v in cub.items(): print(f'   cubic {k}: A {v["A_enum"]:.6g} P2 {v["P2"]} ({v["rel"]})')
    print(f'   native-only classes: {extra}', flush=True)
res['compare'] = dict(max_rel_diff=maxrel, unmatched=unmatched, ok=bool(maxrel <= 0.05 and not unmatched),
                      note='A-enumerator sums entries before the sup; P2 takes per-key sups then folds: A <= P2 expected')
res['args'] = L.args_block([args.u2, os.path.abspath(__file__), os.path.join(HERE, 'rampphi.py'), os.path.join(HERE, 't4lib.py'),
                            L.T2 + '/rigor/ramp_rep.py', L.T2 + '/rigor/arlib.py', RR.L.BASE + '/rigor/weights_exact.py'], dict(SUPLOW=['11/100', U2['global_box']['M1z']]))
L.dump(res, os.path.join(L.OUTD, 'u5_enum.json'))
print('max rel diff', maxrel, 'unmatched', unmatched, ' wall %.1f s' % (time.time() - T0))
