"""Lemma N3' (LEMMA_F.md 6.1): symbolic check of the integration-by-parts identities used by n3far.py.
For each level the difference of the two integrands is checked to be an exact x-derivative (Euler operator = 0 with respect to
every function involved): the weak form at j = 1, 2, 3 (with the level's alpha_1) and the simplified Psi_{j,j+1}, Psi_{j,j+2};
the level-0 form; and the transport identity 2 T w_j (g w_{j+1} + (j+1) g' w_j) = ((2j+1) g' T - g T') w_j^2 + (T g w_j^2)'.
usage: python3 -I n3wcheck.py   (about a minute)"""
import sympy as sp
x = sp.symbols('x')
T, A, B, a0, ab, W, G = [sp.Function(n)(x) for n in ('T', 'A', 'B', 'a0', 'ab', 'W', 'G')]
d = lambda f, k=1: sp.diff(f, x, k)
def euler(dens, funcs, order=8):
    out = []
    for f in funcs:
        e = sum(((-1) ** k) * d(sp.diff(dens, d(f, k)) if k else sp.diff(dens, f), k) for k in range(order + 1))
        out.append(sp.simplify(sp.expand(e)))
    return out
# generic level j >= 1: lhs = 2 (T W)'' (A W'' + a1 W' + a0 W + ab); W = w_j
for name, a1 in (('j3', 2 * d(A) + B), ('j2', d(A) + B), ('j1', B)):
    lhs = 2 * d(T * W, 2) * (A * d(W, 2) + a1 * d(W) + a0 * W + ab)
    phi = 4 * d(T) * a1 - 2 * d(T, 2) * A - 2 * d(T) * d(A)
    P45 = 2 * T * a1; P35 = 2 * T * a0 + 2 * d(T, 2) * A - phi
    P34 = 4 * d(T) * a0 + 2 * d(T, 2) * a1 - d(phi); P33 = 2 * d(T, 2) * a0
    rhs = 2 * T * A * d(W, 2) ** 2 + P45 * d(W) * d(W, 2) + P35 * W * d(W, 2) + P34 * W * d(W) + P33 * W ** 2 \
        + 2 * (T * d(W, 2) + 2 * d(T) * d(W) + d(T, 2) * W) * ab
    ok = all(e == 0 for e in euler(sp.expand(lhs - rhs), [W, T, A, B, a0, ab]))
    # the simplified forms used in n3far
    if name == 'j3': simp = d(T) * (4 * a0 - 6 * d(A, 2) - 4 * d(B)) - 2 * d(T, 2) * B + 2 * d(T, 3) * A
    if name == 'j2': simp = d(T) * (4 * (a0 - d(B)) - 2 * d(A, 2)) + 2 * d(T, 2) * (d(A) - B) + 2 * d(T, 3) * A
    if name == 'j1': simp = d(T) * (4 * a0 - 4 * d(B) + 2 * d(A, 2)) + d(T, 2) * (4 * d(A) - 2 * B) + 2 * d(T, 3) * A
    s35 = {'j3': 2 * T * a0 + 4 * d(T, 2) * A - d(T) * (6 * d(A) + 4 * B), 'j2': 2 * T * a0 + 4 * d(T, 2) * A - d(T) * (2 * d(A) + 4 * B),
           'j1': 2 * T * a0 + 4 * d(T, 2) * A - d(T) * (4 * B - 2 * d(A))}[name]
    print(name, 'weak form identity:', ok, ' Psi_{j,j+1} simplified:', sp.simplify(P34 - simp) == 0, ' Psi_{j,j+2}:', sp.simplify(P35 - s35) == 0)
# level 0: lhs = 2 T W (A W''' + B W'' + C)' with W = w0  (w3 = W''', w2 = W''),  C a function
C = sp.Function('C')(x)
lhs0 = 2 * T * W * d(A * d(W, 3) + B * d(W, 2) + C)
P12 = 2 * T * (d(A) - B); P02 = 4 * d(T, 2) * A + 2 * d(T) * (2 * d(A) - B); P01 = 2 * d(T, 3) * A + 4 * d(T, 2) * d(A) + 2 * d(T) * d(A, 2)
rhs0 = 2 * T * A * d(W, 2) ** 2 + P12 * d(W) * d(W, 2) + P02 * W * d(W, 2) + P01 * W * d(W) - 2 * (T * d(W) + d(T) * W) * C
print('j0 weak form identity:', all(e == 0 for e in euler(sp.expand(lhs0 - rhs0), [W, T, A, B, C])))
# strong form of a transport g w0 at level j: 2 T w_j (g w0)^(j+1); top two terms -> ((2j+1) g' T - g T') w_j^2
for j in range(4):
    lhs = 2 * T * d(W, 0) * sum(sp.binomial(j + 1, k) * d(G, k) * d(W, 1 - k) for k in range(2))   # w_j = W, w_{j+1} = W'
    rhs = ((2 * j + 1) * d(G) * T - G * d(T)) * W ** 2
    lhs = 2 * T * W * (G * d(W) + (j + 1) * d(G) * W)
    print('transport j=%d:' % j, all(e == 0 for e in euler(sp.expand(lhs - rhs), [W, T, G])))
