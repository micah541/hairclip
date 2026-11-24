#!/usr/bin/env bash
# usage: lj.sh [--measure] relative/path/File.lean   (relative to the CWFlow repo root, forward slashes)
# Runs the repo's ../tools/leanjob.py (same as scripts/compile1.sh) with Windows python, LEAN_NUM_THREADS=4.
cd /mnt/c/Users/micahw/agents/lean/CWFlow || exit 2
args=""
for a in "$@"; do args="$args ${a//\//\\}"; done
LINES=${LINES_SHOWN:-80}
timeout ${LJ_TIMEOUT:-3000} cmd.exe /c "cd /d C:\\Users\\micahw\\agents\\lean\\CWFlow && set LEAN_NUM_THREADS=${LJ_THREADS:-4}&& set LINES_SHOWN=$LINES&& set LEAN_START_COMMIT_MB=${LEAN_START_COMMIT_MB:-0}&& set LEAN_MIN_COMMIT_MB=${LEAN_MIN_COMMIT_MB:-1200}&& set LEAN_START_WAIT_S=${LEAN_START_WAIT_S:-600}&& python ..\\tools\\leanjob.py $args" 2>&1 | tr -d '\r'
exit ${PIPESTATUS[0]}
