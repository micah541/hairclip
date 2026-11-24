"""Rigorous zero count of W01 along the inward path from MON lines (enclosures over subintervals).
On a run of subintervals where the W01-enclosure contains 0, W01 is strictly monotone if the W02-enclosure
(dW01/dtau = -W02) has one fixed sign on the whole run; then the run contains at most one zero, and exactly one
iff the definite signs before and after the run differ."""
import sys
def sgn(lo, hi):
    return 1 if lo > 0 else (-1 if hi < 0 else 0)
def count(lines):
    segs = []
    for l in lines:
        s = l.split()
        if s[0] != "MON": continue
        a, b = float(s[2]), float(s[3]); w01 = (float(s[4]), float(s[5])); w02 = (float(s[6]), float(s[7]))
        segs.append((a, b, sgn(*w01), sgn(*w02), w01, w02))
    zeros, last, run, ok = [], None, [], True
    for seg in segs:
        if seg[2] != 0:
            if run:
                s2 = {r[3] for r in run}
                if 0 in s2 or len(s2) != 1:
                    ok = False; print("UNRESOLVED run", run[0][0], run[-1][1], s2)
                if last is not None and last != seg[2]:
                    zeros.append((run[0][0], run[-1][1]))
                elif last is not None and last == seg[2]:
                    print("run without sign change (consistent only if no zero):", run[0][0], run[-1][1])
                run = []
            elif last is not None and last != seg[2]:
                zeros.append((seg[0], seg[0]))       # sign change between adjacent definite subintervals
            last = seg[2]
        else:
            run.append(seg)
    if run:
        ok = False; print("path ends inside an unresolved run")
    return zeros, ok, segs
if __name__ == "__main__":
    lines = open(sys.argv[1]).read().splitlines()
    zeros, ok, segs = count(lines)
    print(f"subintervals: {len(segs)}; zeros of W01 (sigma intervals): {zeros}; resolved: {ok}")
    print("first/last sign:", segs[0][2], segs[-1][2])
