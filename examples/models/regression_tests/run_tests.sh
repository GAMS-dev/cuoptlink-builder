#!/bin/bash
# Runs all self-checking regression models with the cuOpt solver.
# Each model aborts (non-zero GAMS return code) if a result deviates from the reference.
cd "$(dirname "$0")"
fails=0
for model in *.gms; do
    name="${model%.gms}"
    workdir=$(mktemp -d)
    if gams "$PWD/$model" curdir="$workdir" o="$workdir/$name.lst" lo=2 lf="$workdir/$name.log" > /dev/null 2>&1; then
        echo "[PASS] $name"
        rm -rf "$workdir"
    else
        echo "[FAIL] $name (see $workdir/$name.lst)"
        fails=$((fails + 1))
    fi
done
echo "$fails failure(s)"
exit $fails
