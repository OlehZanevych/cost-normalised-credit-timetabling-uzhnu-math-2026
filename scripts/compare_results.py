#!/usr/bin/env python3
"""Compare two result directories record by record, ignoring what is allowed to differ.

    python3 scripts/compare_results.py results results-rerun [--subset]

Every CSV present in both directories is compared.  A record is identified by its key columns (the
columns before the first measured one: n, seed, arm / credit, kappa, cost_ratio, rep, ...), and every
other column must agree exactly, apart from

  * `elapsed_ms`, wall-clock time, which no rerun reproduces; and
  * the whole record of a time-normalised run (`credit == time`), whose credit divides by measured
    microseconds, so its trajectory depends on the machine and its load.  These are counted and
    reported, not compared.

With --subset the second directory may hold only part of the design (a spot check); without it,
both must hold the same set of records.  Exits non-zero on any disagreement.
"""
import csv
import os
import sys

IGNORED = {"elapsed_ms"}
# The first measured column of each file; the columns before it identify a record.
FIRST_MEASURED = {
    "ablation.csv": "hard", "covered.csv": "hard", "credit.csv": "hard",
    "twoarm.csv": "work_budget", "fixedq.csv": "n2", "instances.csv": "classes",
}


def load(path, first):
    with open(path, newline="") as f:
        r = csv.DictReader(f)
        cols = r.fieldnames
        keys = cols[:cols.index(first)]
        rows = {}
        for row in r:
            rows[tuple(row[k] for k in keys)] = row
    return cols, keys, rows


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    subset = "--subset" in sys.argv
    if len(args) != 2:
        print(__doc__)
        return 2
    a_dir, b_dir = args
    bad = 0
    for name in sorted(FIRST_MEASURED):
        pa, pb = os.path.join(a_dir, name), os.path.join(b_dir, name)
        if not (os.path.exists(pa) and os.path.exists(pb)):
            if os.path.exists(pa) != os.path.exists(pb) and not subset:
                print(f"{name}: present in only one directory")
                bad += 1
            continue
        cols_a, keys, ra = load(pa, FIRST_MEASURED[name])
        cols_b, _, rb = load(pb, FIRST_MEASURED[name])
        if cols_a != cols_b:
            print(f"{name}: different columns")
            bad += 1
            continue
        missing = [k for k in rb if k not in ra]
        extra = [k for k in ra if k not in rb]
        if missing or (extra and not subset):
            print(f"{name}: {len(missing)} records only in {b_dir}, {len(extra)} only in {a_dir}")
            bad += 1
        same = differ = timed = 0
        for k, row in rb.items():
            if k not in ra:
                continue
            if row.get("credit") == "time":
                timed += 1
                continue
            diff = [c for c in cols_a if c not in IGNORED and ra[k][c] != row[c]]
            if diff:
                differ += 1
                if differ <= 3:
                    print(f"{name}: record {k} differs in {diff[:6]}")
            else:
                same += 1
        if differ:
            bad += 1
        note = f", {timed} time-normalised records not compared" if timed else ""
        print(f"{name}: {same} identical, {differ} different{note}")
    print("all compared records identical" if bad == 0 else f"{bad} file(s) disagree")
    return 0 if bad == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
