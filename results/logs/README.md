# Logs of the recorded run and of its verification

## The recorded run (25 October 2026)

Host: a Linux container, two cores of an Intel Xeon at 2.1 GHz, Ubuntu 24.04.5, GCC 13.3, CMake
Release build (`-O3`), x86-64. Two shell jobs ran concurrently, one per core; each solver run uses one
worker. The commands, in the order run (`build/experiments` built from the source of the time):

```
job A:  ./build/experiments ablation --seeds 10 --budget 10000000 --results results
job B:  ./build/experiments covered  --seeds 10 --budget 10000000 --results results
        ./build/experiments credit   --seeds 10 --budget 10000000 --results results
        ./build/experiments twoarm   --seeds 50 --results results
        ./build/experiments fixedq   --seeds 2000 --results results
```

After the second pre-submission review two operators were corrected (see the README), and the
experiments they can affect were recorded again, one job each, with the corrected source:

```
job C:  ./build/experiments covered  --seeds 10 --budget 10000000 --results results   (winfix fix)
job D:  ./build/experiments ablation --seeds 10 --budget 10000000 --results results   (no-cluster fix)
        ./build/experiments instancestats --seeds 10 --results results
```

In both re-recordings every record that the correction cannot affect reproduced byte for byte apart
from elapsed time (80 of 120 `covered` records, 180 of 200 `ablation` records), and the new records
replaced the old files; no earlier result file is kept. Each job's console output was one line per
subcommand (`ablation done`, `covered done`, ...), followed by `DONE`.

## Verification of the complete set (25 October 2026, evening)

See `../PROVENANCE.md` for the context. Outputs, verbatim:

Full re-run of the design on the same host, compared with `scripts/compare_results.py results <rerun>`:

```
ablation.csv: 200 identical, 0 different
covered.csv: 120 identical, 0 different
credit.csv: 80 identical, 0 different, 20 time-normalised records not compared
fixedq.csv: 40000 identical, 0 different
instances.csv: 20 identical, 0 different
twoarm.csv: 1500 identical, 0 different
all compared records identical
```

arm64 (Apple-silicon Mac, Linux VM, Ubuntu 22.04, GCC 11.4, `-std=c++20 -O3 -ffp-contract=off`),
instance 2 at n = 800, compared with `--subset`:

```
ablation.csv: 10 identical, 0 different
covered.csv: 6 identical, 0 different
credit.csv: 4 identical, 0 different, 1 time-normalised records not compared
all compared records identical
```

and `twoarm.csv`, `fixedq.csv`, `instances.csv` byte-identical to the files here.

Self-test, statistics and data checks (`scripts/run_all.sh --quick`, final source):

```
planted schedules feasible ......... ok
instance files round-trip .......... ok
evaluators agree term by term ...... ok
winfix proposals within move ....... ok (21436 proposals, 0 outside, identity case ok)

self-test passed
statistics checks passed
shipped files match data/SHA256SUMS (100 files)
the generator reproduces every shipped file
```
