import re,sys
t=open(sys.argv[1]).read()
body=re.sub(r'(?<!\\)%.*','',t)
body=body.replace(r'\{','').replace(r'\}','')
d=0;bad=0
for c in body:
    d+= c=='{'; d-= c=='}'
    if d<0: bad=1;d=0
print('braces net',d,'neg',bad)
st=[];err=[]
for m in re.finditer(r'\\(begin|end)\{([^}]*)\}',body):
    if m.group(1)=='begin': st.append(m.group(2))
    else:
        if not st or st[-1]!=m.group(2): err.append(m.group(0))
        else: st.pop()
print('env open left',st,'mismatch',err)
labels=set(re.findall(r'\\label\{([^}]*)\}',body))
refs=set(re.findall(r'\\(?:ref|eqref)\{([^}]*)\}',body))
print('undefined refs',sorted(refs-labels))
print('dollars even per paragraph:',all(p.replace(r'\$','').count('$')%2==0 for p in re.split(r'\n\s*\n',body)))
print('PASS' if d==0 and not bad and not st and not err and not (refs-labels) else 'FAIL')
