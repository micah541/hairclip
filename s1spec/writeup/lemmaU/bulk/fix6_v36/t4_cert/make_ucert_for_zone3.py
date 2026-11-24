"""Write the zone3 --ucert input from an assembled U_cert: copies ufine.Bhat to supersolution.Bhat (zone3.py's loader key).
usage: python3 -B make_ucert_for_zone3.py <U_cert.json> <out.json>    (replaces the earlier hand-made copy step)"""
import sys, os, json, hashlib, platform
src, dst = sys.argv[1], sys.argv[2]
U = json.load(open(src))
out = dict(label='cert', provisional=False,
           note='U-fine Bhat (ufine.Bhat of %s), stored under supersolution.Bhat for zone3.py --ucert' % os.path.basename(src),
           supersolution={'Bhat': U['ufine']['Bhat']},
           args=dict(argv=sys.argv, env={k: os.environ.get(k) for k in ('PYTHONHASHSEED', 'PYTHONDONTWRITEBYTECODE')}, python=platform.python_version(),
                     inputs={os.path.abspath(src): hashlib.sha256(open(src, 'rb').read()).hexdigest(),
                             os.path.abspath(__file__): hashlib.sha256(open(__file__, 'rb').read()).hexdigest()}))
json.dump(out, open(dst, 'w'), indent=1, sort_keys=True)
