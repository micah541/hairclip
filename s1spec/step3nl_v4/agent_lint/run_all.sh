#!/bin/bash
cd /home/claude/s1spec/step3nl/agent_lint
for c in "3 12" "2.5 12" "2.5 10" "2 12" "2 10" "1.5 12" "1.5 10" "1.5 8"; do
  set -- $c; echo "== PX=$1 XA=$2"; python lint_cert_p.py --PX $1 --XA $2
done
echo ALLDONE
