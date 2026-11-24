from ibp_sym import *
g = Y[0]*dC[4]
print(dC[4])
x=g
for s in range(4):
    x=Dx(x); print(s+1, sorted(map(str,x.free_symbols & set(a))), len(x.args))
