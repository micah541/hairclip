"""Run only the unconditional integral pass on a saved rate_table (checkpoint).  Usage: ADD=abcd FIX=1 THETA=0.3 DELTA=0.2 python3 unc_run.py rate_table_X.json"""
import sys, json, os
import closure
res = json.load(open(sys.argv[1])); th = float(os.environ.get("THETA", 0.3))
lf = open(sys.argv[1].replace("rate_table", "unc_log").replace(".json", ".txt"), "w")
def log(s): print(s, flush=True); lf.write(s + "\n"); lf.flush()
res["unc"] = closure.unconditional(res, th, float(os.environ.get("DELTA", 0.2)), closure.Oracle(5e-4, 0.005, log), log)
json.dump(res, open(sys.argv[1], "w"), indent=1)
