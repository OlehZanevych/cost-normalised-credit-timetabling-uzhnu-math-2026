# Cost-normalised credit can favour cheap operators

Code, input data and per-run records for the article *Cost-normalised credit can favour cheap
operators: an operator-portfolio study for university course timetabling* (O. B. Zanevych,
V. M. Kukharskyy, Ivan Franko National University of Lviv), prepared for the *Scientific Bulletin of
Uzhhorod University. Series of Mathematics and Informatics*, 2026. Repository:
<https://github.com/OlehZanevych/cost-normalised-credit-timetabling-uzhnu-math-2026>.

The article makes two elementary statements and measures their consequences. For a fixed mixture of
a cheap and an expensive operator, the number of expensive draws a work budget buys is bounded by
`(1-q)W/A <= E N2(W) < (1-q)(W+omega_max)/A` (A the mean work per draw, q the cheap share), its
long-run rate `(1-q)/A` falls as q grows, and the expected work until the expensive operator first
proposes a state the cheap one cannot is `A/((1-q)p)`. And, for constant per-draw costs, dividing an
operator's credit by the work of its draws is per-application credit on a reward divided by the cost
(for variable costs the credit tends to mean reward over mean cost); in a deterministic
unclipped mean-field model of the weight update, the article gives the exact condition under which
this raises the cheap operator's share. Whether the shift harms the search is an empirical question, and
the solver experiments measure it at a fixed work budget. Coverage (an operator whose every proposal
another operator can also make) does not by itself decide whether an operator should be kept.

In the recorded study, cost normalisation lowers the share of draws of the two large neighbourhoods
from 75.5 % to 45.7 % at n = 3200, but only time normalisation and uniform selection are
significantly worse than per-application credit at that size; a cheap covered operator (`winfix`) is
significantly worse in one of four size--credit cells and better in none. All numbers in the article
are generated from `results/` by `scripts/analyse.py`.

## Requirements

- a C++20 compiler (GCC 11 or newer has been used) and, optionally, CMake 3.16 or newer;
- Python 3.8 or newer, standard library only, for the analysis and the checks;
- for the article itself (not part of this repository): pdfLaTeX with the `cm-super` fonts.

## Experiments

| subcommand | what it measures | article |
|---|---|---|
| `ablation` | leave-one-out over seven operators, plus the per-operator cost and reward | Table 2, Sect. 8.3 |
| `covered` | adding `winfix` (covered by `move`) or `dayfix` (aimed; coverage not established), under per-application and cost-normalised credit | Table 4, Sect. 8.5 |
| `credit` | five credit regimes: none, per draw, per bucket, per 30 buckets, per microsecond | Table 3, Sect. 8.4 |
| `twoarm` | the allocation effect in a two-armed model, against the mean-field prediction | Table 1, Fig. 1, Sect. 8.2 |
| `fixedq` | Theorem 1 at a fixed mixture: draws bought by a budget, work to the first exclusive proposal | Sect. 8.1 |
| `instancestats` | the realised composition of the generated instances (class kinds, domain sizes) | Sect. 7 |
| `instances` | writes the instance files of `data/` and stops; runs no search | --- |

Every solver run uses **one worker**, no elite adoption, and a fixed budget of **10 000 000 bucket
recomputations** (`--budget`), counted by the evaluator and including construction, the escape phases
and draws that produced no candidate. A run stops at the first check after the budget is reached;
the work actually spent is recorded beside the budget. A single-worker run under a non-time credit
rule is a deterministic function of its seeds; time-normalised runs are not, because their credit
divides by measured microseconds, so their trajectories depend on the machine and its load even
though the work counter itself is exact.

`twoarm` and `fixedq` need no timetabling instance and run in seconds.

Each subcommand writes one CSV to `results/` (or `--results DIR`), with a header row and one line
per run: the size and instance,
the budget and work spent, the validator's counts (hard, soft, objective, unplaced, a flag for
breaches of the hard filters), the escape counts, weight updates and how many ended at a bound, and
per operator the draws, applications, improvements, work, gain and final weight. Nothing is
aggregated by the C++ code.

```bash
./build/experiments twoarm   --seeds 50
./build/experiments fixedq   --seeds 2000
./build/experiments instancestats --seeds 10
./build/experiments ablation --seeds 10 --budget 10000000
./build/experiments covered  --seeds 10 --budget 10000000
./build/experiments credit   --seeds 10 --budget 10000000
```

Without `--seeds` and `--budget` each subcommand runs exactly this design; an unknown option is an
error rather than silently ignored. A solver design can be split across processes with
`--sizes 800` / `--sizes 3200` and `--seedfrom S` (instances `S..seeds`): each run depends only on its
own size and seeds, so the parts concatenated are the same records.

On two cores of an Intel Xeon at 2.1 GHz (Ubuntu 24.04, GCC 13.3, `-O3`), with two experiments
running concurrently, one solver run takes 10--15 s and the full design about 1.5 hours.

### Seeds

Everything is a deterministic function of these seeds (xoshiro256**, `src/core/rng.hpp`):

| record | instance | search / random stream |
|---|---|---|
| `ablation` | generator seed = instance 1..10, size 800 or 3200 | 20260912 + instance, same for every arm |
| `covered` | as above | 777 + instance, same for every arm and credit rule |
| `credit` | as above | 5150 + instance, same for every regime |
| `twoarm` | --- | 9001 + 31 rep + 7 ratio + kappa, shared by the three credit modes |
| `fixedq` | --- | 424242 + 101 rep + 7 omega2 + 1000 q + W |
| `instances` | as above | --- |

## Checking the implementation first

```bash
./build/experiments selftest        # scratch files go to the system temp directory
```

Four checks, each a fraction of a second. The planted schedule of every generated instance is re-scored by
the independent validator and must be feasible, so that a residual hard cost is never the
generator's fault. An instance written to a file and read back must score exactly as the original
did, so that the files shipped in `data/` carry the instances the experiments used. And the incremental
evaluator is compared term by term with the independent validator along a random walk of
placements, so that two implementations of the objective that share no code cannot be wrong in the
same way. Finally, every proposal of `winfix` is checked to be one that `move` can also make (a
single movable class moved to a slot of its domain in the room it holds, never the unchanged
schedule): on a hand-built configuration in which an earlier version returned the unchanged schedule,
and along random walks on generated instances. This is Proposition 1 of the article, checked on the
code. `scripts/run_all.sh` runs it before anything else, together with

```bash
python3 scripts/check_stats.py
```

which checks the statistics in `scripts/ttstats.py` against their definitions rather than against
another library: the exact Wilcoxon signed-rank test is compared with an enumeration of its own null
distribution over several hundred random samples, the sign test and the effect size against values
that can be computed by hand, the chi-squared survival function against the tabulated critical
values, Friedman's statistic (with its correction for ties) against its closed form, and the
distribution-free interval for a median against its order statistics. `ttstats.py` uses the standard library
alone, so the analysis has no dependency to install and no version that can drift; the price is that
its tests have to be written, and this is them.

## Building

Requires a C++20 compiler and CMake 3.16 or newer. Nothing else: no external libraries, no package
manager, no network access.

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build
```

Without CMake, the same build is

```bash
mkdir -p build && g++ -std=c++20 -O3 -ffp-contract=off -Isrc src/main.cpp src/core/*.cpp -o build/experiments -pthread
```

The build uses ISO C++ and `-ffp-contract=off`, so that no compiler fuses a multiply and an add into
one differently rounded instruction; with it, the records reproduce across architectures (checked:
x86-64 with GCC 13.3 and arm64 with GCC 11.4 give identical records, see `results/PROVENANCE.md`).

## Reproducing the paper

```bash
python3 scripts/analyse.py  # the published records in results/ -> tables/*.tex, figures/*.dat
scripts/run_all.sh          # re-run every experiment into results-rerun/ and compare with results/
```

`analyse.py` writes `tables/` and `figures/` inside this folder; when the folder sits inside the
article's own directory (a `main.tex` one level up), it writes them beside the article instead, which
is how the article was built (`pdflatex main.tex`, twice; the journal class needs pdfLaTeX and the
`cm-super` fonts). Every number in the article is produced this way: no value in any table or figure
is typed by hand, and the numbers quoted in the prose reach the text as macros in
`tables/facts.tex`.

`run_all.sh` builds, runs the self-test, `check_stats.py` and `scripts/check_data.sh`, runs the full
design into `results-rerun/` and then `scripts/compare_results.py results results-rerun`, which must
report every record identical apart from elapsed time and the time-normalised runs (their credit
divides by measured microseconds). `--overwrite` writes into `results/` instead, replacing the
published records.

`run_all.sh` accepts `--quick` to run a reduced design (fewer seeds, smaller budgets) that is enough
to check that everything works end to end; it takes a few minutes, against hours for the full
design. It writes to `results-quick/` rather than `results/`, so a smoke test cannot
overwrite the records the paper was built from --- point `analyse.py --results results-quick --partial`
at it to see the tables a reduced design produces (`--partial` is needed because the analysis checks
the records against the paper's declared design and a reduced design does not match it).

## Input data

The experiments generate their instances in memory, since an instance is a deterministic function of
its size and seed. The same instances are shipped in `data/` as plain-text files (gzipped), for a
solver that shares no code with this one: all 20 instances the experiments use (n = 800 and 3200,
seeds 1--10), plus n = 200, 400 and 1600, each with its planted feasible schedule
(`*.planted.txt.gz`). `data/SHA256SUMS` holds the checksums of the uncompressed files, and

```bash
scripts/check_data.sh
```

checks them and regenerates every file with the current generator, which must reproduce each one
byte for byte. The two models (`twoarm`, `fixedq`) have no input data: their parameters are in
`src/main.cpp` and their random streams are seeded as in the table above.

## Analysis checks

`analyse.py` refuses to read a results file in any other directory than the published `results/`
that was written in the last two minutes (an experiment appends to its CSV as it runs, so a fresh file
may be unfinished; an analysis of a re-run started immediately after it ends is therefore refused
once, and succeeds two minutes later). The published records are always read, whatever their
modification time after a clone or copy. It then checks
the records against the declared design -- every size, arm and repetition count the paper reports,
no cell outside it, no duplicate key, no non-finite value -- and stops otherwise. `--partial` turns these stops into warnings, and
`TT_MIN_AGE=0` lifts the age check, while iterating on the analysis itself; neither is for the paper.

`analyse.py --out DIR` writes the tables and figure data elsewhere, for example
`python3 scripts/analyse.py --results results-rerun --out results-rerun` for a re-run.

## Layout

```
src/core/        the solver kernel, shared with the other applications in this series
  model.*        instances: entities, domains, the compressed tick axis
  state.*        the schedule, its (entity, day) bucket decomposition, incremental evaluation
  search.*       construction, acceptance, the nine neighbourhoods, the weight update, the escape
  validate.*     an INDEPENDENT scorer, sharing no code with state.cpp -- see below
  generate.*     instances built backwards around a hidden feasible schedule
  io.*           the plain-text instance format
  flow.hpp       Dinic's maximum flow, for the matching bounds
  hungarian.hpp  the assignment problem, for the permutation operator
  mask.hpp       128-bit tick masks
  rng.hpp        xoshiro256**, so that a seed reproduces a run exactly
  selftest.*     the four self-test checks
src/             this application's experiment drivers (main.cpp)
scripts/         run_all.sh, analyse.py, ttstats.py (statistics, standard library only),
                 check_stats.py, check_data.sh, compare_results.py
data/            the instances and planted schedules, gzipped, with SHA256SUMS
results/         one CSV per experiment, one line per run, PROVENANCE.md, and logs/ (commands,
                 host and verification output of the recorded run)
LICENSE          MIT
```

## Two things that make the numbers trustworthy

**An independent validator.** `src/core/validate.cpp` re-reads a schedule from scratch and counts
every penalty term and every hard filter by the most direct method the definitions allow -- sorting
and scanning, never a bitmask. It deliberately shares no code with the incremental evaluator in
`state.cpp`, so that an error in one cannot hide behind an identical error in the other. Every
number in the paper comes from the validator, never from the search's own counters.

**Instances with a known answer.** The generator builds a valid schedule first, walking the week slot
by slot and placing each class into resources that are free at that moment and reachable in the time
available, and then reads the instance off it. So a schedule with no hard violation provably exists
for every instance, and anything the search cannot reach is a property of the search rather than of
the data. The planted schedule is a feasibility witness, not a good timetable.

## Budgets

Every solver budget is counted in **bucket recomputations**, never in seconds. The count is
deterministic and machine-independent, so a single-worker run under a non-time credit rule
reproduces exactly on other hardware; time-normalised runs do not, since their credit divides by
measured microseconds. It is a reproducible
proxy for effort, not a constant amount of CPU time: the cost of one recomputation grows with the
occupancy of the bucket. Wall-clock is recorded in every results file alongside it.

## Corrections after the second review

Two operators' behaviour was corrected after the second pre-submission review, and the affected
experiments were re-run (all other records reproduced byte for byte): `winfix` now moves only a class
taught in the week whose gap it closes, to a slot taught in that week, and never returns the
unchanged schedule (checked by the self-test); and without the community partition `ruin` and the
perturbations now draw uniformly among their remaining selectors instead of drawing an empty one.

## What git tracks

Everything needed to rebuild the article's evidence: the sources, the build files, the scripts, all
of `data/` (gzipped instances, planted schedules, `SHA256SUMS`) and all of `results/` (the per-run
records, `PROVENANCE.md`, `logs/`). `.gitignore` excludes only what is regenerated: build output,
Python bytecode, `tables/` and `figures/` written by the analysis, a reader's `results-rerun/` and
`results-quick/`, and uncompressed copies of the instance files.

## Licence

Released under the MIT License; see [LICENSE](LICENSE).
