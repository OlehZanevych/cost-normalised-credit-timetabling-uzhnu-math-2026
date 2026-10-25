#!/usr/bin/env bash
# Run every experiment of this article and compare the new records with the published ones.
#
#   scripts/run_all.sh              full design into results-rerun/, then compare with results/
#   scripts/run_all.sh --quick      reduced design into results-quick/ (a smoke test, minutes)
#   scripts/run_all.sh --overwrite  full design into results/ itself (replaces the published records)
#
# The full design takes about 1.5 hours on two cores (each subcommand is one process; running them
# in parallel shortens it).
set -euo pipefail
cd "$(dirname "$0")/.."

# Always configure and build: an executable left over from older sources would otherwise be run
# silently.  When nothing has changed this takes a second.
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release > /dev/null
cmake --build build

QUICK=0
OVERWRITE=0
for arg in "$@"; do
  [ "$arg" = "--quick" ] && QUICK=1
  [ "$arg" = "--overwrite" ] && OVERWRITE=1
done

# Before any experiment: the generator, the instance file format, the two independent evaluators and
# Proposition 1 (selftest), the statistics (check_stats), and the shipped instance files
# (check_data).  A failure here would make every number below meaningless.
mkdir -p build/selftest
./build/experiments selftest --results build/selftest
python3 scripts/check_stats.py > /dev/null && echo "statistics checks passed"
scripts/check_data.sh

# Never write into results/ unless asked: it holds the records the paper was built from.
RES=results-rerun
[ $QUICK -eq 1 ] && RES=results-quick
[ $OVERWRITE -eq 1 ] && RES=results
mkdir -p "$RES"
run() { ./build/experiments "$@" --results "$RES"; }

if [ $QUICK -eq 1 ]; then
  run twoarm   --seeds 5
  run fixedq   --seeds 200
  run instancestats --seeds 3
  run ablation --seeds 3 --budget 300000
  run covered  --seeds 3 --budget 300000
  run credit   --seeds 3 --budget 300000
  echo
  echo "results written to $RES/; to see the tables: python3 scripts/analyse.py --results $RES --partial --out $RES"
  exit 0
fi

# Every solver run: one worker, a fixed budget of 10 000 000 bucket recomputations.  These are the
# commands that produced every record in results/ and every number in the paper.
run twoarm   --seeds 50
run fixedq   --seeds 2000
run instancestats --seeds 10
run ablation --seeds 10 --budget 10000000
run covered  --seeds 10 --budget 10000000
run credit   --seeds 10 --budget 10000000

echo
if [ "$RES" != results ]; then
  # Everything must agree with the published records except elapsed time and the time-normalised
  # runs, whose credit divides by measured microseconds.
  python3 scripts/compare_results.py results "$RES"
fi
echo "results written to $RES/; tables from them: python3 scripts/analyse.py --results $RES --out $RES"
