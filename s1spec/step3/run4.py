import sys, time, gen_hs2 as G
tag = sys.argv[1]
if tag.startswith("Z"):
    gi = tag[1:]; path = f"in4{tag}.txt"; G.build_ZMi(path, gi)
else:
    gj, gi = tag[1], tag[2]; path = f"in4{tag}.txt"; G.build_Qji(path, gj, gi)
out = G.run(path); open(f"out4{tag}.txt", "w").write(out)
print(tag, "RESULT OK" if "RESULT OK" in out else out[-300:])
