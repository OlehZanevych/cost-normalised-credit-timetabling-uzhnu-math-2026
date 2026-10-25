#!/usr/bin/env bash
# Check the shipped instance files: (1) each decompresses to the content recorded in
# data/SHA256SUMS, and (2) the current generator reproduces every one of them byte for byte.
# The experiments generate their instances in memory from (size, seed); the files are shipped so that
# a solver sharing no code with this one can read the same instances.
set -euo pipefail
cd "$(dirname "$0")/.."
EXE=${EXE:-build/experiments}
tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT
mkdir -p "$tmp/shipped" "$tmp/regen"
for f in data/*.txt.gz; do gunzip -c "$f" > "$tmp/shipped/$(basename "$f" .gz)"; done
(cd "$tmp/shipped" && sha256sum --quiet -c "$OLDPWD/data/SHA256SUMS")
echo "shipped files match data/SHA256SUMS ($(ls "$tmp/shipped" | wc -l) files)"
"$EXE" instances --data "$tmp/regen" > /dev/null
diff -rq "$tmp/shipped" "$tmp/regen"
echo "the generator reproduces every shipped file"
