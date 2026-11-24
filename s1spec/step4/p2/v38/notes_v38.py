"""v38 copy of p2/v33/notes_v33.py (zone3_interface of the (A) chain v37, output out_v38; nothing else changed).  v3.3 task 6: covering notes (no computation beyond asserts).  Writes p2/out_v38/notes_v38.json.
1. Far family of v3.3 (task 5, zone3_interface.json cone_family.covering.far_members: L = 113, X = 2, centres >= 490) is a subset of the
   v3.2 far family (centres >= 376) whose hand-off items are certified and pinned in FINAL_manifest.json: (H5) far data t9_far_v32.json
   (x >= 150), H4 h4loc_v32.json (FAR_M = 0.15), Lemma O H1 arc by arc t9_lemmaO_far.json (every centre at base arclength >= 375.125, ok).
   A member centred at r_c >= 490 has its cylinder in [r_c - 226, r_c + 226], r >= 264 >= 150: inside the region of t9_far's bounds.
   So no far stage is rerun.
2. t9_h6_v32 (v3.2 H6 data check) is superseded by task 5's H6 (zone3_interface.json cone_family.H6, min margin recorded below);
   its pinned output stays in FINAL_manifest.json unchanged and is not an input of v3.3.
3. The (H3) windows, alpha_1, sqrt(Lambda), the (H5) far data and every other v3.2 hand-off item are unchanged (Z5_REPORT section 4)."""
import sys, os, json, hashlib
sha = lambda p: hashlib.sha256(open(p, 'rb').read()).hexdigest()
O = '/home/claude/s1spec/step4/p2/out'; Z = '/home/claude/s1spec/writeup/lemmaU/bulk/fix13_v37chain/t5_zone3/out/zone3_interface.json'
z = json.load(open(Z)); cov = z['cone_family']['covering']; h6 = z['cone_family']['H6']
lf = json.load(open(f'{O}/t9_lemmaO_far.json')); fv = json.load(open(f'{O}/t9_far_v32.json'))
far_min_centre = 490; v32_min_centre = 376
assert 'centres >= 490' in cov['far_members'] and far_min_centre >= v32_min_centre and far_min_centre - 2 * 113 >= 150
assert lf.get('ok') is True
out = dict(far_family=dict(v33='L = 113, X = 2, centres >= 490 (' + cov['far_members'] + ')', v32='centres >= 376', subset=True,
                           cylinder_base_min_r=far_min_centre - 2 * 113, t9_far_region='x >= 150', far_covers_from=cov['far_covers_from'],
                           pinned=dict(t9_far_v32=f'{O}/t9_far_v32.json', h4loc_v32='/home/claude/s1spec/step4/p2/t8se/constants/h4loc_v32.json',
                                       t9_lemmaO_far=f'{O}/t9_lemmaO_far.json'), t9_lemmaO_far_ok=lf['ok'], rerun=False),
           h6=dict(t9_h6_v32='superseded by task 5 H6 (zone3_interface.json cone_family.H6)',
                   task5_at=('at_Bhat_ref' if 'at_Bhat_ref' in h6 else None), task5_min_margin_at_Bhat_ref=h6.get('at_Bhat_ref', {}).get('min_margin')),
           unchanged='(H3) windows, alpha_1, sqrt(Lambda), (H5) far data, H4 far, far Lemma O H1: v3.2 values, FINAL_manifest.json',
           args=dict(argv=sys.argv[1:], env={}, inputs={p: sha(p) for p in (os.path.abspath(__file__), Z, f'{O}/t9_lemmaO_far.json', f'{O}/t9_far_v32.json')}))
json.dump(out, open('/home/claude/s1spec/step4/p2/out_v38/notes_v38.json', 'w'), indent=1); print(json.dumps(out)[:600])
