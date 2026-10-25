# Provenance of the records in `results/`

These are the records every table, figure and number of the article is generated from
(`python3 scripts/analyse.py`). They were all produced on 25 October 2026.

## The records

| file | records | design | SHA-256 |
|---|---:|---|---|
| `ablation.csv` | 200 | 10 arms x 2 sizes x 10 instances, W = 10 000 000 | `37410abdfdc4e27f52451e362d0e968cf4571cbdcb34bbf590e07db07a876462` |
| `covered.csv` | 120 | 3 arms x 2 credit rules x 2 sizes x 10 instances | `ae8b4fe4e89186773cd8f6cc5c76ccba6e42a20532b7aea5533001070f72c15f` |
| `credit.csv` | 100 | 5 regimes x 2 sizes x 10 instances | `f2a19411bc7d04944af38597e06f8f98d561319b2a791503f8c2150aa128edcf` |
| `twoarm.csv` | 1 500 | 2 kappa x 3 credit modes x 5 cost ratios x 50 repetitions | `370aaf5e082e217ce2d766f5c564a5fb8a41f2a0d5c6654edc90f4c8a7662334` |
| `fixedq.csv` | 40 000 | 2 omega2 x 5 q x 2 budgets x 2 000 repetitions | `6562ce689d1411a80ef143a23645cb65311f1533e163d69a1eb259d6f7d63792` |
| `instances.csv` | 20 | 2 sizes x 10 instances | `8e2b03ddcee7dbc0f7046719274353c0754024d88ba4668df8884418629dcb24` |

Commands: those listed in `README.md` (the article's design, which is also each subcommand's
default); the exact jobs, their order, the host and the verification output are in `logs/README.md`.

## Where and how they were produced

- **Host:** a Linux container, two cores of an Intel Xeon at 2.1 GHz, Ubuntu 24.04, GCC 13.3,
  CMake Release build (`-O3`), x86-64. Two experiments ran concurrently, one per core; every solver
  run uses one worker.
- **Order and source versions.** `twoarm`, `fixedq` and `credit` were recorded first. After the second
  pre-submission review two operator corrections were made (see `README.md`), and the experiments
  they can affect were recorded again with the corrected source: `covered` (the `winfix` correction)
  and `ablation` (the selector correction for `no-cluster`); `instances.csv` was added then. All other
  records of those two files reproduced byte for byte apart from elapsed time.
- **Later changes to the source** (input/output and command line only: the results directory is
  created when missing and write failures stop the program; unknown options are rejected; the
  article's design is the default; `--sizes` and `--seedfrom` split a design; the self-test writes its
  scratch files to the system temp directory; `-ffp-contract=off` and ISO C++ in the build) do not
  touch any computation.

## Verification of the complete set (25 October 2026, evening)

1. **Full re-run, x86-64.** Every experiment of the design was run again from scratch into a separate
   directory and compared with `scripts/compare_results.py`:

   ```
   ablation.csv: 200 identical, 0 different
   covered.csv: 120 identical, 0 different
   credit.csv: 80 identical, 0 different, 20 time-normalised records not compared
   fixedq.csv: 40000 identical, 0 different
   instances.csv: 20 identical, 0 different
   twoarm.csv: 1500 identical, 0 different
   ```

   "Identical" means every column except `elapsed_ms`. The final source (with the input/output
   changes above) reproduces `twoarm.csv`, `fixedq.csv` and `instances.csv` byte for byte and a
   sample of `covered` records exactly.
2. **arm64.** The final source built on the Mac's Linux VM (Apple arm64, Ubuntu 22.04, GCC 11.4,
   `-std=c++20 -O3 -ffp-contract=off`) reproduces `twoarm.csv`, `fixedq.csv` and `instances.csv` byte
   for byte, and 20 solver records (10 ablation arms, 6 `covered` and 4 `credit` runs of instance 2 at
   n = 800) exactly. A GNU-mode build without `-ffp-contract=off` also reproduced the 6 `covered`
   records; the flag is a safeguard, not a correction.
3. **Input data.** `scripts/check_data.sh`: the 100 shipped instance files match `data/SHA256SUMS`, and
   the generator reproduces each of them byte for byte.
4. **Code checks.** `selftest` (four checks) and `check_stats.py` pass; a build with
   `-Wall -Wextra -Wshadow -Wpedantic` gives no warning; AddressSanitizer and
   UndefinedBehaviorSanitizer report nothing on every subcommand (reduced designs); clang-tidy
   (bugprone, performance) reports only benign signed/unsigned character comparisons of day indices.

## What does not reproduce exactly

The 20 **time-normalised** records (`credit == time`) divide each operator's reward by measured
microseconds, so their trajectories depend on the machine and its load. In the re-run none of the 20
objectives was identical to the published one (ratios 0.54--1.41, median 1.01). The article's
conclusion about this regime at n = 3200 reproduced (worse than per-application credit on all ten
instances, adjusted p = 0.008; median paired ratio 1.31 against the published 1.20); at n = 800 it is
not significantly different from per-application credit in either run, but which regime has the lowest
median there changed (per-application credit in the published records, time normalisation in the
re-run). `elapsed_ms` is wall-clock time and never reproduces.
