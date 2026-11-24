"""Monte-Carlo sanity test of the componentwise Neumann bound  |x_i| <= t_i + tau (N d)_i/(1-nu)  with random G-hat = G + E, |E_kj| <= Ebound, random R."""
import json, numpy as np
rng = np.random.default_rng(1)
L = json.load(open('mod_lemma_eta1e-2.json'))
G = np.array([[(a + b) / 2 for a, b in row] for row in L['G']])
H = np.linalg.inv(G).T; Hab = np.abs(H)
I = json.load(open('lem_integrals.json'))
AM, W1, AR = 1e-3, 1e-3, 1e-3
Eb = np.array([[I[f'e1_{j}'] * AM + I[f'Pq1_{j}'] * W1 + I[f'Pq0_{j}'] * AM * AR for j in (1, 2, 3)],
               [I[f'e2_{j}'] * AM + I[f'Pc1_{j}'] * W1 + I[f'Pc0_{j}'] * AM * AR for j in (1, 2, 3)],
               [2 * np.sqrt(171.83119) * n * AM for n in (np.sqrt(G[0, 0]), np.sqrt(G[1, 1]), 0.019476463)]])
N = np.array([[sum(Hab[i, j] * Eb[k, j] for j in range(3)) for k in range(3)] for i in range(3)])
ev, V = np.linalg.eig(N); v = np.abs(V[:, np.argmax(ev.real)].real); d = v + 1e-3 * v.max()
nu = max((N @ d) / d); print("nu", nu)
worst = 0
for _ in range(20000):
    E = np.where(rng.random((3, 3)) < 0.5, Eb * np.sign(rng.standard_normal((3, 3))), Eb * rng.uniform(-1, 1, (3, 3)))
    Rb = rng.uniform(0, 1, 3) * 10 ** rng.uniform(-6, 0, 3); R = Rb * rng.choice([-1, 1], 3)
    x = np.linalg.solve((G + E).T, R)
    t = Hab @ Rb; tau = max(t / d)
    bound = t + tau * (N @ d) / (1 - nu)
    worst = max(worst, max(np.abs(x) / bound))
print("max |x_i|/bound over trials:", worst)
