"""Scratch [est, float]: kc2_model's patched U_model2 head (window-uniform data), with the chi ramp moved:
 outer term DT6 -> Delta^R for [r_chi, 379.5] times D(r_chi) phi; present centres up to r_chi/1.5; margins against (0.15, 0.52)."""
import sys, os, json, math
sys.dont_write_bytecode = True
SRC='/home/claude/s1spec/writeup/lemmaU/phase2a/U_ramp/U_model2.py'; EXP='/home/claude/s1spec/writeup/lemmaU/phase2a/A_repair/export_for_U.json'
def build(rchi, DT6, YC):
    src=open(SRC).read(); head=src.split('# ================= runs =================')[0]
    for a,b in (("        return (math.sqrt(S_X*Dk) + math.sqrt(eps28*S_E*p**(2*k)*dU[0]))**2","        return (math.sqrt(S_X*Dk + cfg.get('C1', (0., 0., 0.))[k - 2]*SIG_X) + math.sqrt(eps28*S_E*p**(2*k)*dU[0]))**2"),
                ("    D = cfg['D']\n","    D = cfg['D']; SIG_X = sweep(Dr*yR/28, yc, mu, nu, a)\n"),
                ("    return out, dict(Gw=Gw, S_in=S_in, S_E=S_E, S_X=S_X, Xi=Xi","    return out, dict(Gw=Gw, S_in=S_in, S_E=S_E, S_X=S_X, SIG_X=SIG_X, Xi=Xi"),
                ("(1 + g/190.**4)", f"(1 + g/{rchi}**4)"),
                ("YC = (56., 60., 70., 84., 100., 126.6, 150., 200., 253.)", f"YC = {YC!r}")):
        assert head.count(a)==1,a; head=head.replace(a,b)
    head=head.replace("ap.parse_args()","ap.parse_args(['--inputs', EXP])").replace("os.path.abspath(__file__)","SRC")
    g=dict(__name__='kc2',__file__=SRC,EXP=EXP,SRC=SRC)
    import io, contextlib
    with contextlib.redirect_stdout(io.StringIO()): exec(compile(head,SRC,'exec'),g)
    g['DT6']=DT6
    return g
X=json.load(open(EXP)); EA=tuple(X['EA']); gp=tuple(X['gp']); C0=tuple(X['D_from_any_tau_a_with_E_tau_a_le_EA']); C1=tuple(X['D_slope'])
def run(tag, rchi, DT6, KS):
    ymax=rchi/1.5; YC=tuple([56.,60.,70.,84.,100.,126.6]+[y for y in (140.,150.,160.) if y<ymax]+([ymax] if ymax>126.7 else []))
    g=build(rchi,DT6,YC); cfg=g['cfg_new'](dict(EA=EA,gp=gp,D=C0),'data56','widened (0.15,0.52), U .09, M0z .16',KS=KS); cfg['C1']=C1
    g['cfg_box_is_wide'][0]=True
    lt=g['late'](cfg); rows=g['present'](cfg)
    print(f"{tag} r_chi={rchi} DT6_2={DT6[2]} KS={KS}: LATE M1 x{lt['M1']:.3f} M0(.15) x{lt['M0w']:.2f} B2 {lt['B'][2]:.4f} outer {lt['ev'][2]['outer']:.4f}")
    for r in rows:
        print(f"     yc {r['yc']:7.2f}: U1(1.5yc) {1.3056*1.5e-4*(1.5*r['yc'])**1.0246:.4f}  f {r['f']:.4f} (M0(.15) x{r['M0w']:.2f})  yf' {r['yfp']:.4f} (M1 x{r['M1']:.2f})  B=({r['B'][0]:.2e},{r['B'][1]:.2e},{r['B'][2]:.4f})")
for KS in (1.33, 2):
    run('r1 ', 190., (9.86e-08, 1.04e-06, 0.0161), KS)
    run('r2 ', 253., (1.005e-07, 1.085e-06, 0.01745), KS)
