#!/usr/bin/env bash
# usage: raw.sh CWFlow/Curve/X.lean   (unguarded compile with own watchdog; see scratch_curve/rawlean.py)
cd /mnt/c/Users/micahw/agents/lean/CWFlow || exit 2
f=${1//\//\\}
cmd.exe /c "cd /d C:\\Users\\micahw\\agents\\lean\\CWFlow && set RAW_MIN_MB=60&& python scratch_curve\\rawlean.py $f" 2>&1 | tr -d '\r'
