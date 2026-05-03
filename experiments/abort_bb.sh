#!/usr/bin/env bash
# abort_bb.sh — kill the Bayesian NB2 head-to-head fit (06_bb_replication.py).
# Safe to run repeatedly. Reports what it killed.

killed=0
for name in nb2_glm.exe model.exe; do
  mapfile -t pids < <(tasklist 2>/dev/null | awk -v n="$name" '$1 == n {print $2}')
  for pid in "${pids[@]}"; do
    [ -n "$pid" ] || continue
    echo "killing $name PID $pid"
    taskkill //F //PID "$pid" >/dev/null 2>&1 && killed=$((killed + 1))
  done
done

mapfile -t fit_pids < <(wmic process where "name='python.exe'" get ProcessId,CommandLine /format:csv 2>/dev/null \
  | grep -i "06_bb_replication" \
  | awk -F',' '{print $(NF)}' | tr -d '\r' | grep -E '^[0-9]+$')
for pid in "${fit_pids[@]}"; do
  [ -n "$pid" ] || continue
  echo "killing fit launcher python PID $pid"
  taskkill //F //PID "$pid" >/dev/null 2>&1 && killed=$((killed + 1))
done

if [ "$killed" -eq 0 ]; then
  echo "no BB-fit processes found"
else
  echo "aborted $killed process(es)"
fi
