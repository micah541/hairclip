"""Illustrative figure for main.tex (floating point; not part of any certificate).

Panel (a): the initial curve Gamma0 of step 4 (t2/gamma0_final.json), evaluated in doubles.
Panel (b): zoom on the bend of Gamma0.
Panel (c): the shrinker S1 near its tip, from the centres of the certified enclosures
           s1_tight/ext_tip/prof_tip.txt (sigma in [sigma_tip, 8]) and s1_tight/ext813/ext813_prof.txt
           (sigma in [8, 13]), reflected in the symmetry line; dashed: the asymptotic rays.
Usage: python3 -I make_fig.py   (writes fig_curves.pdf next to this script)
"""
import json
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
G0 = "/home/claude/s1spec/step4/t2/gamma0_final.json"
PTIP = "/home/claude/s1_tight/ext_tip/prof_tip.txt"
P813 = "/home/claude/s1_tight/ext813/ext813_prof.txt"


def gamma0(N=2**18):
    d = json.load(open(G0))
    b = np.array([float.fromhex(h) for h in d["b_hex"]])
    lo, hi = (float(x) for x in d["b_j_closing_interval"])
    b[d["j"] - 1] = 0.5 * (lo + hi)                     # b_2: centre of its certified interval
    L0 = float.fromhex(d["L0_hex"])
    X0 = float.fromhex(d["X0_hex"][0]) + 1j * float.fromhex(d["X0_hex"][1])
    c = np.zeros(N, dtype=complex)
    c[1:len(b) + 1] = b
    s = (np.fft.ifft(c) * N).imag                       # sum_n b_n sin(n alpha_k), alpha_k = 2 pi k / N
    alpha = 2 * np.pi * np.arange(N) / N
    theta = -np.pi / 2 + alpha + s
    e = np.exp(1j * theta)
    h = 2 * np.pi / N
    e_ext = np.append(e, e[0])
    X = X0 + (L0 / (2 * np.pi)) * np.concatenate([[0], np.cumsum(0.5 * h * (e_ext[1:] + e_ext[:-1]))])
    return X, abs(X[-1] - X[0])


def rows(path):
    out = []
    with open(path) as f:
        for line in f:
            if line.startswith("END"):
                break
            v = [float(x) for x in line.split()]
            out.append(v[:13])
    return np.array(out)


def s1_branch():
    a = rows(PTIP)                                      # tau in [0, tau_tip], sigma = 8 - tau
    b = rows(P813)                                      # tau in [-5, 0]
    def pts(r):
        tc = r[:, 2]
        p = 0.5 * (r[:, 3] + r[:, 4]); q = 0.5 * (r[:, 5] + r[:, 6]); th = 0.5 * (r[:, 11] + r[:, 12])
        return 8 - tc, np.exp(1j * th) * (q + 1j * p)
    sa, za = pts(a); sb, zb = pts(b)
    sig = np.concatenate([sa, sb]); z = np.concatenate([za, zb])
    o = np.argsort(sig)
    sig, z = sig[o], z[o]
    P = z[0]                                            # tip (smallest sigma)
    rot = -abs(P) / P                                   # tip -> (-D, 0), symmetry line -> real axis
    return sig, z * rot


def main():
    X, gap = gamma0()
    sig, z = s1_branch()
    fig = plt.figure(figsize=(7.2, 4.0))
    ax1 = fig.add_axes([0.07, 0.70, 0.90, 0.26])
    ax2 = fig.add_axes([0.07, 0.07, 0.42, 0.52])
    ax3 = fig.add_axes([0.57, 0.07, 0.40, 0.52])
    idx = np.linspace(0, len(X) - 1, 40000).astype(int)
    ax1.plot(X.real[idx], X.imag[idx], lw=0.8, color="k")
    ax1.set_aspect("equal"); ax1.set_title("(a)", fontsize=8, loc="left")
    ax1.tick_params(labelsize=7)
    # zoom: around the left-most part of the curve (the bend), chosen from the data
    xmin = X.real.min()
    Z = X.copy()
    Z[X.real >= xmin + 3.0] = np.nan                    # break the line outside the zoom window
    ax2.plot(Z.real, Z.imag, lw=0.8, color="k")
    ax2.set_xlim(xmin - 0.1, xmin + 2.9)
    ax2.set_aspect("equal"); ax2.set_title("(b)", fontsize=8, loc="left")
    ax2.tick_params(labelsize=7)
    ax3.plot(z.real, z.imag, lw=0.9, color="k")
    ax3.plot(z.real, -z.imag, lw=0.9, color="k")
    ang = np.deg2rad(106.67)
    R = np.array([0, 16])
    for sgn in (1, -1):
        ax3.plot(R * np.cos(ang), sgn * R * np.sin(ang), ls="--", lw=0.6, color="0.5")
    ax3.plot([0], [0], marker="+", color="0.3", ms=6)
    ax3.set_xlim(-6.5, 3.0); ax3.set_ylim(-5.5, 5.5)
    ax3.set_aspect("equal"); ax3.set_title("(c)", fontsize=8, loc="left")
    ax3.tick_params(labelsize=7)
    fig.savefig(os.path.join(HERE, "fig_curves.pdf"))
    print("closure defect (float):", gap, " S1 sigma range:", sig.min(), sig.max(), " tip:", z[0])


if __name__ == "__main__":
    main()
