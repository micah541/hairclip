"""[kit step 4b] t9_ift_tight.json for win.py: the hull of alpha, xc, lam, Theta over the NB t_c-balls of W' = [1 - 1/512, 1] (out/t9_h1_K_coh_tip32.json,
the step-3 IFT with prof_tip + ext813, covering the whole sub-window).  In the reference proof tree this file was the IFT at the W' midpoint with
prof_tight, widened by the drift over W'; here it is the certified hull of the ball runs (tighter in Theta)."""
import os, json
HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(HERE, '..', 'out')
NB = int(os.environ.get('H1NB', '32'))
R = [json.load(open(os.path.join(OUT, f't9_h1_{K}_coh_tip32.json'))) for K in range(NB)]
hull = lambda key: [min(r[key][0] for r in R), max(r[key][1] for r in R)]
out = dict(alpha=hull('alpha'), xc=hull('xc'), lam=hull('lam'), Theta=hull('Theta'), tc_ball=[1 - 1 / 512, 1.0],
           note=f'kit 4b: hull of the {NB} t_c-ball IFT runs of W\' (h1.py, prof_tip + ext813)')
json.dump(out, open(os.path.join(OUT, 't9_ift_tight.json'), 'w'), indent=1); print(json.dumps(out))
