import sys, time, gen_hs2 as G
mode = sys.argv[1]; t0 = time.time()
if mode == "ZM": G.build_ZM("in3ZM.txt"); out = G.run("in3ZM.txt")
else: j = int(mode[2:]); G.build_QS(f"in3{mode}.txt", j); out = G.run(f"in3{mode}.txt")
open(f"out3{mode}.txt", "w").write(out)
print(("RESULT OK" if "RESULT OK" in out else out[-500:]), f"({time.time()-t0:.1f}s)")
