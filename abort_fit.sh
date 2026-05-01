#!/usr/bin/env bash
# abort_fit.sh — kill any running Collatz Stan fit + its python launcher.
# Safe to run repeatedly; reports what it killed.
set -u

killed=0

# Kill model.exe (the Stan executable, runs all chains in one process)
mapfile -t model_pids < <(tasklist 2>/dev/null | awk '/^model\.exe/ {print $2}')
if [ "${#model_pids[@]}" -gt 0 ]; then
  echo "killing model.exe: ${model_pids[*]}"
  for pid in "${model_pids[@]}"; do
    taskkill //F //PID "$pid" >/dev/null 2>&1 && killed=$((killed + 1))
  done
fi

# Kill the parent fit.py python.exe by command line match
mapfile -t fit_pids < <(wmic process where "name='python.exe'" get ProcessId,CommandLine /format:csv 2>/dev/null \
  | grep -i "fit\.py" \
  | awk -F',' '{print $(NF)}' | tr -d '\r' | grep -E '^[0-9]+$')
if [ "${#fit_pids[@]}" -gt 0 ]; then
  echo "killing fit.py launcher: ${fit_pids[*]}"
  for pid in "${fit_pids[@]}"; do
    taskkill //F //PID "$pid" >/dev/null 2>&1 && killed=$((killed + 1))
  done
fi

if [ "$killed" -eq 0 ]; then
  echo "no Collatz fit processes found"
else
  echo "aborted $killed process(es)"
fi
