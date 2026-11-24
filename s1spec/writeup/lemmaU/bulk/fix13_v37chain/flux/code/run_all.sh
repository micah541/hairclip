#!/bin/bash
# X8 driver: full production pass X1..X7 into out/ and texts/ (one process at a time, nice 19, PYTHONHASHSEED=0, python3 -B).
# usage: F2_V36RUN=<final v36 snapshot dir> bash code/run_all.sh <tag>      (copies out/ texts/ to snap/<tag>/ afterwards)
set -e
D=/home/claude/s1spec/writeup/lemmaU/bulk/fix13_v37chain/flux
tag=$1; [ -n "$tag" ] || { echo "tag?"; exit 1; }
[ ! -e $D/snap/$tag ] || { echo "snapshot $tag exists"; exit 1; }
unset F2_OUT F2_TX F2_DRYRUN
for s in flux2_forms.py flux2_coef.py chart4_shell.py inputs_manifest.py flux2_data.py flux2_assemble.py flux2_sanity.py gen_texts.py; do
  echo "== $s $(date -Is)"
  bash $D/code/run_stage.sh $s > /dev/null || { echo "STAGE FAILED: $s"; exit 2; }
  grep -E "Maximum resident|Elapsed" $D/out/${s%.py}.time | tr '\n' ' '; echo
done
mkdir -p $D/snap/$tag; cp -a $D/out $D/texts $D/snap/$tag/
echo "done $tag $(date -Is)"
