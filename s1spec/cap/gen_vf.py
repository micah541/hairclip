# regenerate evans_seg.cpp from the template
s = open('evans_seg.tmpl.cpp').read()
W = ["a1r","a1i","a2r","a2i","a3r","a3i","b2r","b2i","b3r","b3i","c3r","c3i"]
LR = "(lra+t*dr)"; LI = "(lia+t*di)"
def fields(scaled):
    s2 = "*iS2" if scaled else ""
    A = f"(0.25-(3*m*m+k*p){s2})"; B1 = f"(5*k*m{s2}+0.25*q)"; K2 = f"k*k{s2}"
    prof = ["k*q", f"-1-k*p{s2}", "-m", "-0.25*p"]
    return prof + ["-a2r", "-a2i", "-(a3r+b2r)", "-(a3i+b2i)",
     f"-(-{B1}*a1r - {K2}*a2r + b3r)", f"-(-{B1}*a1i - {K2}*a2i + b3i)", "-b3r", "-b3i",
     f"-(-{A}*a1r + ({LR}*a1r - {LI}*a1i) - {K2}*b2r + c3r)", f"-(-{A}*a1i + ({LR}*a1i + {LI}*a1r) - {K2}*b2i + c3i)",
     f"-(-{A}*a2r + ({LR}*a2r - {LI}*a2i) + {B1}*b2r)", f"-(-{A}*a2i + ({LR}*a2i + {LI}*a2r) + {B1}*b2i)"]
funs = fields(False); sfun = fields(True)
damp = funs[:4] + [f"({f}) - (0.25*q)^(1/3)*{w}" for f, w in zip(funs[4:], W)]
var = "var:p,q,k,m," + ",".join(W) + ",t;"
VF = 'static const char* VF = "par:lra,lia,dr,di;' + var + 'fun:' + ",".join(funs + ["0"]) + ';";'
VFB = 'static const char* VFB = "par:lra,lia,dr,di;' + var + 'fun:' + ",".join(damp + ["0"]) + ';";'
VFS = 'static const char* VFS = "par:lra,lia,dr,di,iS2;' + var + 'fun:' + ",".join(sfun + ["0"]) + ';";'
open('evans_seg.cpp', 'w').write(s.replace('@VF@', VF).replace('@VFB@', VFB).replace('@VFS@', VFS))
