#!/bin/bash
# usage: lane.sh "arm,W,w arm,W,w ..."
for cfg in $1; do
  IFS=, read arm W w <<< "$cfg"
  tag="s_L${arm}_W${W}_w${w}"
  ARM=$arm WID=$W SLOT=$w AMOLL=0.025 python3 run_hairpin.py $tag 1e-2 16384 450 > /dev/null 2>&1
  echo "$tag done" >> lanes.log
done
