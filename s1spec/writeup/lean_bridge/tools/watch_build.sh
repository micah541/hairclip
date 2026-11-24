#!/usr/bin/env bash
# Waits until Windows has enough free commit, then compiles the blocked files (one leanjob at a time).
# Stage A (needs ~8 GB free commit): N1Basic, LiftV, Lift3   (Scalar/Weinstein closure, ~6.4 GB commit)
# Stage B (needs ~11.5 GB free commit): N1, Lift, Criterion  (IsFlow closure, ~9.8 GB estimated)
# usage: watch_build.sh [hours]   logs: lean_bridge/build_stage{A,B}.log, lean_bridge/watch.log
H=${1:-12}
LB=/home/claude/s1spec/writeup/lean_bridge
LJ=$LB/tools/lj.sh
end=$(( $(date +%s) + H*3600 ))
A=0; B=0
log(){ echo "$(date '+%F %T') $*" >> $LB/watch.log; }
log "watcher start (hours=$H)"
while [ $(date +%s) -lt $end ]; do
  free=$(powershell.exe -NoProfile -Command "[int]((Get-CimInstance Win32_OperatingSystem).FreeVirtualMemory/1024)" 2>/dev/null | tr -d '\r')
  if [ -n "$free" ]; then
    if [ $A -eq 0 ] && [ "$free" -ge 8000 ]; then
      log "free=$free MB: stage A"
      LJ_THREADS=1 LINES_SHOWN=200 $LJ CWFlow/Curve/N1Basic.lean CWFlow/Curve/LiftV.lean CWFlow/Curve/Lift3.lean > $LB/build_stageA.log 2>&1
      rc=$?; log "stage A rc=$rc"; grep -q "exit 137\|NOT STARTED" $LB/build_stageA.log || A=1
    elif [ $A -eq 1 ] && [ $B -eq 0 ] && [ "$free" -ge 11500 ]; then
      log "free=$free MB: stage B"
      LJ_THREADS=1 LINES_SHOWN=200 $LJ CWFlow/Curve/N1.lean CWFlow/Curve/Lift.lean CWFlow/Curve/Criterion.lean > $LB/build_stageB.log 2>&1
      rc=$?; log "stage B rc=$rc"; grep -q "exit 137\|NOT STARTED" $LB/build_stageB.log || B=1
    fi
  fi
  [ $A -eq 1 ] && [ $B -eq 1 ] && { log "both stages attempted to completion"; exit 0; }
  sleep 120
done
log "watcher timeout"
