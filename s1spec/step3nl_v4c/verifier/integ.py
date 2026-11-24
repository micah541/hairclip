"""Recompute the time-integral bounds (closure.integrals) from a saved rate_table_*.json.  Usage: ADD=abcd THETA=0.3 python3 integ.py rate_table_X.json"""
import sys, json, os
import closure
res = json.load(open(sys.argv[1])); th = float(os.environ.get("THETA", 0.3))
lf = open(sys.argv[1].replace("rate_table", "integ_log").replace(".json", ".txt"), "w")
def log(s): print(s, flush=True); lf.write(s + "\n"); lf.flush()
res["integrals"] = closure.integrals(res, th, closure.Oracle(5e-4, 0.005, log), log)
json.dump(res, open(sys.argv[1], "w"), indent=1)
