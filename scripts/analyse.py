#!/usr/bin/env python3
"""Tables, figure data and inline numbers for the operator-portfolio and credit-assignment paper.

Reads the five result files of `experiments` (ablation, covered, credit, twoarm, fixedq) and writes
tables/*.tex, tables/facts.tex and figures/*.dat.  No number in the paper is typed by hand.

Quality is compared lexicographically by (unplaced, validator breaches, H, f): feasibility first,
then the objective.  The two are never folded into one number.  A paired comparison uses
  * the exact Wilcoxon signed-rank test on f when every pair is feasible on both sides, and
  * the exact sign test on the lexicographic key otherwise (marked with a dagger in the tables);
the paired effect is the median over instances of f_arm / f_base, with the distribution-free
order-statistic interval for that median (for ten pairs [x(2), x(9)], coverage 97.9 %), and it is
reported only when every pair is feasible.
"""

import argparse
import collections
import math
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ttstats as t  # noqa: E402

OPS = ["move", "swap", "chain", "kempe", "ruin", "repack", "kopt", "dayfix", "winfix"]
ABL_ARMS = ["full", "no-move", "no-swap", "no-chain", "no-kempe", "no-ruin", "no-repack",
            "no-kopt", "no-cluster", "simple"]
COV_ARMS = ["base", "+winfix", "+dayfix"]
CREDITS = ["none", "application", "work", "work-scaled", "time"]
DAGGER = "$^\\dagger$"

# The design, declared rather than read off the records: an experiment whose file lacks a whole
# size, a whole arm or an equal number of repetitions in every cell must fail the checks below,
# which a check derived from the data alone would let pass.
SIZES = [800, 3200]
SOLVER_REPS = 10
TWOARM_KAPPAS = ["300", "60"]
TWOARM_RATIOS = ["2", "5", "10", "20", "40"]
TWOARM_REPS = 50
FIXEDQ_OMEGAS = ["10", "40"]
FIXEDQ_QS = ["0", "0.25", "0.5", "0.75", "0.9"]
FIXEDQ_BUDGETS = ["100", "100000"]
FIXEDQ_REPS = 2000


class Incomplete(Exception):
    pass


# ---------------------------------------------------------------------------------------------
#  Helpers
# ---------------------------------------------------------------------------------------------


def feasible(r):
    return t.num(r, "unplaced") == 0 and t.num(r, "hard") == 0 and t.num(r, "clean") == 1


def lex(r):
    """The lexicographic quality key; smaller is better."""
    return (t.num(r, "unplaced"), 1 - t.num(r, "clean"), t.num(r, "hard"), t.num(r, "objective"))


def f_of(r):
    return t.num(r, "objective")


def big(x, digits=3):
    """A number for a table cell: integers below a million as they are, larger ones as
    $a\\cdot10^{k}$."""
    if x != x:
        return "--"
    if abs(x) >= 1e6:
        e = int(math.floor(math.log10(abs(x))))
        m = x / 10 ** e
        ms = f"{m:.{digits - 1}f}"
        return f"${ms}\\cdot10^{{{e}}}$"
    return t.fmt(round(x), 0)


def thin(n):
    """10000000 -> 10\\,000\\,000."""
    return "{:,}".format(int(n)).replace(",", "\\,")


def pct(x, digits=1):
    return t.fmt(100.0 * x, digits)


def tag(s):
    """A string usable inside a macro name."""
    return s.replace("+", "plus").replace("-", "")


def by_seed(rows):
    return {r["seed"]: r for r in rows}


def compare(base, arm):
    """Paired comparison of `arm` against `base` (lists of rows, paired on the seed).

    Returns a dict: p (two-sided, uncorrected), test ('W' signed-rank on f, or 'S' sign test on
    the lexicographic key), ratio / lo / hi (median of f_arm/f_base and its interval; NaN unless
    every pair is feasible), better / ties / worse (paired counts by the lexicographic key)."""
    a, b = by_seed(base), by_seed(arm)
    keys = sorted(set(a) & set(b), key=int)
    if len(keys) < 5:
        return None
    pa, pb = [a[k] for k in keys], [b[k] for k in keys]
    allfeas = all(feasible(x) for x in pa + pb)
    better, ties, worse = t.wins_losses([lex(x) for x in pb], [lex(x) for x in pa])
    out = {"better": better, "ties": ties, "worse": worse, "pairs": len(keys)}
    if allfeas:
        _, p = t.wilcoxon_signed_rank([f_of(x) for x in pb], [f_of(x) for x in pa])
        ratios = [f_of(y) / max(1e-9, f_of(x)) for x, y in zip(pa, pb)]
        lo, hi, cov = t.median_ci(ratios)
        out.update(p=p, test="W", ratio=t.median(ratios), lo=lo, hi=hi, coverage=cov)
    else:
        n = better + worse
        k = min(better, worse)
        p = min(1.0, 2.0 * sum(math.comb(n, i) for i in range(k + 1)) / 2 ** n) if n else 1.0
        out.update(p=p, test="S", ratio=float("nan"), lo=float("nan"), hi=float("nan"),
                   coverage=float("nan"))
    return out


def ratio_cell(c):
    if c is None or c["ratio"] != c["ratio"]:
        return "--"
    return f"{c['ratio']:.2f}\\,{{\\scriptsize[{c['lo']:.2f}, {c['hi']:.2f}]}}"


def p_cell(p, test):
    return t.fmt_p(p) + (DAGGER if test == "S" else "")


def share(rows, ops):
    """Pooled share of draws of the operators `ops` among all adaptive draws of `rows`."""
    num = sum(t.num(r, f"{o}_draws") for r in rows for o in ops)
    den = sum(t.num(r, f"{o}_draws") for r in rows for o in OPS)
    return num / den if den else float("nan")


def clip_frac(rows):
    up = sum(t.num(r, "weight_updates") for r in rows)
    hits = sum(t.num(r, "ceiling_hits") + t.num(r, "floor_hits") for r in rows)
    return hits / up if up else float("nan")


def check_cells(name, rows, keyfields, expected_cells, reps_field, partial, expected_reps):
    """Every declared cell present with exactly the declared number of runs, no cell outside the
    declared design, no duplicate keys, every numeric field finite."""
    problems = []
    seen = collections.Counter(tuple(r[k] for k in keyfields + [reps_field]) for r in rows)
    dup = [k for k, c in seen.items() if c > 1]
    if dup:
        problems.append(f"{len(dup)} duplicate keys, e.g. {dup[0]}")
    cells = collections.Counter(tuple(r[k] for k in keyfields) for r in rows)
    missing = [c for c in expected_cells if c not in cells]
    if missing:
        problems.append(f"{len(missing)} expected cells missing, e.g. {missing[0]}")
    extra = [c for c in cells if c not in set(expected_cells)]
    if extra:
        problems.append(f"{len(extra)} cells outside the declared design, e.g. {extra[0]}")
    wrong = [c for c in expected_cells if c in cells and cells[c] != expected_reps]
    if wrong:
        problems.append(f"{len(wrong)} cells without {expected_reps} runs, e.g. {wrong[0]}: "
                        f"{cells[wrong[0]]}")
    for r in rows:
        for k, v in r.items():
            try:
                x = float(v)
            except (TypeError, ValueError):
                continue
            if not math.isfinite(x):
                problems.append(f"non-finite {k} in {name}")
                break
    if problems:
        msg = f"{name}: " + "; ".join(problems)
        if partial:
            print("  WARNING " + msg)
        else:
            raise Incomplete(msg)
    return max(cells.values()) if cells else 0


# ---------------------------------------------------------------------------------------------
#  Main
# ---------------------------------------------------------------------------------------------


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", default="results")
    ap.add_argument("--out", default=None,
                    help="where tables/ and figures/ are written; defaults to the article directory "
                         "above this one when it holds a main.tex, and otherwise to this directory")
    ap.add_argument("--partial", action="store_true",
                    help="warn instead of stopping when a design is incomplete (for iterating on "
                         "the analysis while experiments run; never for the paper)")
    args = ap.parse_args()
    # The published records in this application's own results/ are finished by definition, but a
    # fresh copy (git clone, tar, rsync) gives them a new modification time; the "written in the last
    # two minutes" guard is meant for a re-run still in progress, so it does not apply to them.
    published = os.path.realpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "results"))
    if os.path.realpath(args.results) == published:
        os.environ.setdefault("TT_MIN_AGE", "0")
    if args.out is None:
        args.out = ".." if os.path.exists(os.path.join("..", "main.tex")) else "."
    R, TAB, FIG = args.results, os.path.join(args.out, "tables"), os.path.join(args.out, "figures")
    os.makedirs(TAB, exist_ok=True)
    os.makedirs(FIG, exist_ok=True)
    t.reset_outputs(TAB, FIG)
    facts = {}
    budgets, overshoot = set(), []

    def solver_rows(name):
        path = os.path.join(R, name)
        if not t.complete(path):
            if args.partial:
                return None
            raise Incomplete(f"{name} missing or still being written")
        rows = t.read(path)
        for r in rows:
            budgets.add(int(t.num(r, "work_budget")))
            overshoot.append((t.num(r, "work") - t.num(r, "work_budget")) / t.num(r, "work_budget"))
        return rows

    # --- Ablation ------------------------------------------------------------------------------
    rows = solver_rows("ablation.csv")
    if rows:
        sizes = SIZES
        reps = check_cells("ablation", rows, ["n", "arm"],
                           [(str(n), a) for n in sizes for a in ABL_ARMS], "seed", args.partial,
                           SOLVER_REPS)
        facts["ablRunsPerCell"] = str(reps)
        facts["ablArms"] = str(len(ABL_ARMS) - 1)
        floor = 2.0 / 2 ** reps
        facts["signRankFloor"] = f"{floor:.5f}"
        facts["holmFloorAblation"] = f"{floor * (len(ABL_ARMS) - 1):.4f}"
        by = collections.defaultdict(list)
        for r in rows:
            by[(int(r["n"]), r["arm"])].append(r)
        cell, infeasible = {}, []
        for n in sizes:
            base = by[(n, "full")]
            comps = {a: compare(base, by[(n, a)]) for a in ABL_ARMS[1:] if by.get((n, a))}
            comps = {a: c for a, c in comps.items() if c}
            adj = dict(t.holm([(a, c["p"]) for a, c in comps.items()]))
            for a, c in comps.items():
                c["padj"] = adj[a]
            better, worse = [], []
            for a in ABL_ARMS:
                v = by.get((n, a), [])
                c = comps.get(a)
                feas = sum(1 for r in v if feasible(r))
                fm = t.median([f_of(r) for r in v])
                cell[(n, a)] = [big(fm), ratio_cell(c),
                                p_cell(c["padj"], c["test"]) if c else "--"]
                if feas < len(v):
                    infeasible.append(f"{len(v) - feas} of the \\textsf{{{a}}} runs at $n={n}$")
                T = f"abl{n}{tag(a)}"
                facts[T + "Feas"] = str(feas)
                facts[T + "Obj"] = big(fm)
                if c:
                    facts[T + "P"] = t.fmt_p(c["padj"])
                    facts[T + "WL"] = f"{c['better']}:{c['worse']}"
                    if c["ratio"] == c["ratio"]:
                        facts[T + "Ratio"] = f"{c['ratio']:.2f}"
                        facts[T + "Lo"] = f"{c['lo']:.2f}"
                        facts[T + "Hi"] = f"{c['hi']:.2f}"
                    if c["padj"] < 0.05:
                        (better if c["better"] > c["worse"] else worse).append(a)
            facts[f"ablBetterSig{n}"] = ", ".join(f"\\textsf{{{a}}}" for a in better) or "none"
            facts[f"ablBetterSigCount{n}"] = str(len(better))
            facts[f"ablWorseSig{n}"] = ", ".join(f"\\textsf{{{a}}}" for a in worse) or "none"
            facts[f"ablWorseSigCount{n}"] = str(len(worse))
            facts[f"ablFeasFull{n}"] = str(sum(1 for r in base if feasible(r)))
            # Friedman across all arms, on the lexicographic key (ranks, so no scalarisation).
            seeds = set.intersection(*[{r["seed"] for r in by[(n, a)]} for a in ABL_ARMS])
            blocks = [[lex(by_seed(by[(n, a)])[s]) for a in ABL_ARMS] for s in sorted(seeds, key=int)]
            chi2, df, p, ranks = t.friedman(blocks)
            facts[f"friedmanChi{n}"] = t.fmt(chi2, 1)
            facts[f"friedmanDf{n}"] = str(df)
            facts[f"friedmanP{n}"] = t.fmt_p(p)
        head = (" & " + " & ".join(f"\\multicolumn{{3}}{{c}}{{$n={n}$}}" for n in sizes) + " \\\\\n"
                + "".join(f"\\cmidrule(lr){{{2 + 3 * i}-{4 + 3 * i}}}" for i in range(len(sizes)))
                + "\narm" + " & $f$ & ratio [CI] & $p$" * len(sizes) + " \\\\\n")
        prow = [[f"\\textsf{{{a}}}"] + sum((cell[(n, a)] for n in sizes), []) for a in ABL_ARMS]
        t.write_table(
            os.path.join(TAB, "ablation.tex"),
            "Leave-one-out ablation under cost-normalised credit, at a fixed work budget",
            "tab:ablation", [], prow, align="l" + "rrr" * len(sizes), small=True, header_tex=head,
            note=("%s instances per size, one run each, $W=%s$ bucket recomputations, one worker. "
                  "%s $f$: median objective. "
                  "\\emph{ratio}: median over instances of $f_{\\mathrm{arm}}/f_{\\textsf{full}}$ with "
                  "its distribution-free %.1f\\%% interval $[x_{(2)},x_{(9)}]$, given when every "
                  "pair is feasible. $p$: exact Wilcoxon signed-rank test on $f$, paired on the "
                  "instance, or, marked $\\dagger$, exact sign test on $(\\unp,H,f)$ when a pair "
                  "contains an infeasible run; two-sided, Holm-corrected over the %s comparisons "
                  "of each size. The interval is marginal (not Holm-adjusted, not a test inversion), "
                  "so it can disagree with $p$.") % (reps, thin(max(budgets)),
                                      ("Every run is feasible ($\\unp=0$, $H=0$) except %s."
                                       % ", ".join(infeasible)) if infeasible
                                      else "Every run is feasible ($\\unp=0$, $H=0$).",
                                      100 * (1 - 22 / 1024), facts["ablArms"]))

        # Operator cost and return, in the full portfolio, pooled over both sizes.
        full = [r for r in rows if r["arm"] == "full"]
        totalDraws = sum(t.num(r, f"{o}_draws") for r in full for o in OPS)
        totalGain = sum(t.num(r, f"{o}_gain") for r in full for o in OPS)
        srows, perDraw = [], {}
        for op in OPS:
            draws = sum(t.num(r, f"{op}_draws") for r in full)
            if draws == 0:
                continue
            applied = sum(t.num(r, f"{op}_applied") for r in full)
            improved = sum(t.num(r, f"{op}_improved") for r in full)
            work = sum(t.num(r, f"{op}_work") for r in full)
            gain = sum(t.num(r, f"{op}_gain") for r in full)
            perDraw[op] = work / draws
            srows.append([f"\\textsf{{{op}}}", pct(draws / totalDraws), t.fmt(100 * improved / max(1, applied), 2),
                          t.fmt(work / draws, 1), t.fmt(work / max(1, applied), 1),
                          pct(gain / totalGain)])
            facts[f"op{op}DrawShare"] = pct(draws / totalDraws)
            facts[f"op{op}GainShare"] = pct(gain / totalGain)
            facts[f"op{op}WorkPerDraw"] = t.fmt(work / draws, 1)
            facts[f"op{op}WorkPerApply"] = t.fmt(work / max(1, applied), 0)
        costly = max(perDraw, key=perDraw.get)
        cheap = min(perDraw, key=perDraw.get)
        facts["costliestOperator"] = costly
        facts["cheapestOperator"] = cheap
        facts["costliestWork"] = t.fmt(perDraw[costly], 1)
        facts["cheapestWork"] = t.fmt(perDraw[cheap], 1)
        facts["costRatioObserved"] = t.fmt(perDraw[costly] / perDraw[cheap], 0)
        # Work spent outside the adaptive loop (construction, the escape phases), per size.
        for n in sizes:
            fr = [r for r in full if int(r["n"]) == n]
            tot = sum(t.num(r, "work") for r in fr)
            inloop = sum(t.num(r, f"{o}_work") for r in fr for o in OPS)
            facts[f"outsideLoop{n}"] = t.fmt(100 * (tot - inloop) / tot, 2)
            facts[f"deepMedian{n}"] = t.fmt(t.median([t.num(r, "deep") for r in fr]), 0)
            facts[f"kicksMedian{n}"] = t.fmt(t.median([t.num(r, "kicks") for r in fr]), 0)
            facts[f"freshMedian{n}"] = t.fmt(t.median([t.num(r, "fresh") for r in fr]), 0)
            facts[f"drawsMedian{n}"] = thin(round(t.median([t.num(r, "draws") for r in fr])))
            facts[f"secondsMedian{n}"] = t.fmt(t.median([t.num(r, "elapsed_ms") for r in fr]) / 1000, 1)

    # --- Coverage ------------------------------------------------------------------------------
    rows = solver_rows("covered.csv")
    if rows:
        sizes = SIZES
        reps = check_cells("covered", rows, ["n", "credit", "arm"],
                           [(str(n), c, a) for n in sizes for c in ("application", "work")
                            for a in COV_ARMS], "seed", args.partial, SOLVER_REPS)
        facts["covRuns"] = str(reps)
        by = collections.defaultdict(list)
        for r in rows:
            by[(int(r["n"]), r["credit"], r["arm"])].append(r)
        covcell, worseCells, betterCells = {}, collections.Counter(), collections.Counter()
        covInfeasible = []
        cellsTested = 0
        for n in sizes:
            for credit in ("application", "work"):
                base = by[(n, credit, "base")]
                comps = {a: compare(base, by[(n, credit, a)]) for a in COV_ARMS[1:]}
                adj = dict(t.holm([(a, c["p"]) for a, c in comps.items() if c]))
                cellsTested += 1
                for a in COV_ARMS:
                    v = by[(n, credit, a)]
                    c = comps.get(a)
                    op = a[1:] if a != "base" else None
                    T = f"cov{n}{credit}{tag(a)}"
                    mv = share(v, ["move"])
                    ad = share(v, [op]) if op else float("nan")
                    lns = share(v, ["ruin", "repack"])
                    fm = t.median([f_of(r) for r in v])
                    facts[T + "Obj"] = big(fm)
                    facts[T + "Move"] = pct(mv)
                    facts[T + "Lns"] = pct(lns)
                    facts[T + "Feas"] = str(sum(1 for r in v if feasible(r)))
                    if op:
                        facts[T + "Add"] = pct(ad)
                    pc = "--"
                    if c:
                        pa = adj[a]
                        pc = p_cell(pa, c["test"])
                        facts[T + "P"] = t.fmt_p(pa)
                        facts[T + "WL"] = f"{c['better']}:{c['worse']}"
                        if c["ratio"] == c["ratio"]:
                            facts[T + "Ratio"] = f"{c['ratio']:.2f}"
                            facts[T + "Lo"] = f"{c['lo']:.2f}"
                            facts[T + "Hi"] = f"{c['hi']:.2f}"
                        if pa < 0.05:
                            if c["worse"] > c["better"]:
                                worseCells[op] += 1
                            else:
                                betterCells[op] += 1
                    covcell[(n, credit, a)] = [ratio_cell(c), pc, pct(mv), pct(ad) if op else "--"]
                    for r in v:
                        if not feasible(r):
                            covInfeasible.append(
                                "\\textsf{%s} under %s credit at $n=%d$ (instance %s: $\\unp=%d$, $H=%d$)"
                                % (a, {"application": "per-application", "work": "cost-normalised"}[credit],
                                   n, r["seed"], t.num(r, "unplaced"), t.num(r, "hard")))
        for op in ("winfix", "dayfix"):
            facts[f"{op}WorseCells"] = str(worseCells[op])
            facts[f"{op}BetterCells"] = str(betterCells[op])
            d = sum(t.num(r, f"{op}_draws") for r in rows)
            w = sum(t.num(r, f"{op}_work") for r in rows)
            facts[f"{op}Cost"] = t.fmt(w / d, 1) if d else "--"
        d = sum(t.num(r, "move_draws") for r in rows)
        facts["moveCostCovered"] = t.fmt(sum(t.num(r, "move_work") for r in rows) / d, 1)
        facts["covCells"] = str(cellsTested)
        prow = [[{"application": "appl."}.get(credit, credit) if a == "base" else "", f"\\textsf{{{a}}}"]
                + sum((covcell[(n, credit, a)] for n in sizes), [])
                for credit in ("application", "work") for a in COV_ARMS]
        head = (" & & " + " & ".join(f"\\multicolumn{{4}}{{c}}{{$n={n}$}}" for n in sizes) + " \\\\\n"
                + "".join(f"\\cmidrule(lr){{{3 + 4 * i}-{6 + 4 * i}}}" for i in range(len(sizes)))
                + "\ncredit & arm" + " & ratio [CI] & $p$ & move\\% & added\\%" * len(sizes)
                + " \\\\\n")
        t.write_table(
            os.path.join(TAB, "covered.tex"),
            "Adding \\textsf{winfix} (covered by \\textsf{move}) or \\textsf{dayfix} (aimed, "
            "coverage not established) at a fixed work budget",
            "tab:covered", [], prow, align="ll" + "rrrr" * len(sizes), compact=True, header_tex=head,
            note=("%s instances per size, one run each, $W=%s$, one worker. %s \\emph{ratio}, "
                  "$p$ and $\\dagger$ as in Table~\\ref{tab:ablation}, against \\textsf{base} under "
                  "the same credit rule (\\emph{appl.}: per-application credit); Holm over the two additions of each (size, credit) cell. "
                  "\\emph{move}, \\emph{added}: pooled shares of the draws received by "
                  "\\textsf{move} and by the added operator.")
            % (reps, thin(max(budgets)),
               ("Every run is feasible except one run of " + "; one run of ".join(covInfeasible) + "."
                if covInfeasible else "Every run is feasible.")))

    # --- Credit --------------------------------------------------------------------------------
    rows = solver_rows("credit.csv")
    if rows:
        sizes = SIZES
        reps = check_cells("credit", rows, ["n", "credit"],
                           [(str(n), c) for n in sizes for c in CREDITS], "seed", args.partial,
                           SOLVER_REPS)
        facts["credRuns"] = str(reps)
        by = collections.defaultdict(list)
        for r in rows:
            by[(int(r["n"]), r["credit"])].append(r)
        credcell = {}
        for n in sizes:
            base = by[(n, "application")]
            comps = {g: compare(base, by[(n, g)]) for g in CREDITS if g != "application"}
            adj = dict(t.holm([(g, c["p"]) for g, c in comps.items() if c]))
            meds = {}
            for g in CREDITS:
                v = by[(n, g)]
                c = comps.get(g)
                T = f"cred{n}{tag(g)}"
                fm = t.median([f_of(r) for r in v])
                meds[g] = t.median([lex(r) for r in v]) if False else fm
                lns = share(v, ["ruin", "repack"])
                clip = clip_frac(v) if g != "none" else float("nan")
                facts[T + "Obj"] = big(fm)
                facts[T + "Soft"] = t.fmt(t.median([t.num(r, "soft") for r in v]), 0)
                facts[T + "Lns"] = pct(lns)
                facts[T + "Feas"] = str(sum(1 for r in v if feasible(r)))
                facts[T + "Draws"] = big(t.median([t.num(r, "draws") for r in v]))
                if clip == clip:
                    facts[T + "Clip"] = pct(clip)
                pc = "--"
                if c:
                    pc = p_cell(adj[g], c["test"])
                    facts[T + "P"] = t.fmt_p(adj[g])
                    facts[T + "WL"] = f"{c['better']}:{c['worse']}"
                    if c["ratio"] == c["ratio"]:
                        facts[T + "Ratio"] = f"{c['ratio']:.2f}"
                        facts[T + "Lo"] = f"{c['lo']:.2f}"
                        facts[T + "Hi"] = f"{c['hi']:.2f}"
                credcell[(n, g)] = [ratio_cell(c), pc, pct(lns), pct(clip) if clip == clip else "--"]
            best = min(meds, key=meds.get)
            facts[f"credBest{n}"] = best
            order = sorted(meds, key=meds.get)
            facts[f"credOrder{n}"] = ", ".join(f"\\emph{{{g}}}" for g in order)
        prow = [[g] + sum((credcell[(n, g)] for n in sizes), []) for g in CREDITS]
        head = (" & " + " & ".join(f"\\multicolumn{{4}}{{c}}{{$n={n}$}}" for n in sizes) + " \\\\\n"
                + "".join(f"\\cmidrule(lr){{{2 + 4 * i}-{5 + 4 * i}}}" for i in range(len(sizes)))
                + "\ncredit" + " & ratio [CI] & $p$ & LNS\\% & bound\\%" * len(sizes) + " \\\\\n")
        t.write_table(
            os.path.join(TAB, "credit.tex"),
            "Five credit regimes at a fixed work budget",
            "tab:credit", [], prow, align="l" + "rrrr" * len(sizes), compact=True, sep="4.5pt",
            header_tex=head,
            note=("%s instances per size, one run each, $W=%s$, one worker. \\emph{none}: uniform "
                  "selection; \\emph{application}: reward per draw; \\emph{work}: reward per bucket "
                  "recomputation; \\emph{work-scaled}: reward per %s bucket recomputations; "
                  "\\emph{time}: reward per microsecond; $\\kappa=300$ throughout. \\emph{ratio}, $p$ as in "
                  "Table~\\ref{tab:ablation}, against \\emph{application}, Holm over the four "
                  "comparisons of each size. \\emph{LNS}: share of draws of \\textsf{ruin} and "
                  "\\textsf{repack}; \\emph{bound}: share of weight updates that ended at the floor "
                  "$0.05$ or the ceiling $40$.") % (reps, thin(max(budgets)), "30"))

    allSolver = [r for name in ("ablation.csv", "covered.csv", "credit.csv")
                 for r in (t.read(os.path.join(R, name)) if os.path.exists(os.path.join(R, name)) else [])]
    if allSolver:
        facts["solverRuns"] = str(len(allSolver))
        for k in ("deep", "kicks", "fresh"):
            facts[f"runsWith{k.capitalize()}"] = str(sum(1 for r in allSolver if t.num(r, k) > 0))
        out = [(t.num(r, "work") - sum(t.num(r, f"{o}_work") for o in OPS)) / t.num(r, "work")
               for r in allSolver]
        facts["outsideLoopMaxPct"] = t.fmt(100 * max(out), 2)
    if budgets:
        facts["budgetWork"] = thin(max(budgets))
        facts["budgetWorkMillions"] = t.fmt(max(budgets) / 1e6, 0)
        facts["overshootMaxPct"] = t.fmt(100 * max(overshoot), 3)
        if min(overshoot) < 0:
            raise Incomplete("a solver run stopped before its work budget (time limit?)")

    # --- The generated instances --------------------------------------------------------------
    path = os.path.join(R, "instances.csv")
    if not t.complete(path) and not args.partial:
        raise Incomplete("instances.csv missing or still being written")
    if t.complete(path):
        rows = t.read(path)
        check_cells("instances", rows, ["n"], [(str(n),) for n in SIZES],
                    "seed", args.partial, SOLVER_REPS)
        for k in ("biweekly", "online", "sport", "immovable", "multigroup"):
            facts[f"inst{k.capitalize()}Pct"] = pct(t.mean([t.num(r, k) for r in rows]))
        facts["instRoomDomainMean"] = t.fmt(
            sum(t.num(r, "room_domain_mean") * t.num(r, "room_classes") for r in rows)
            / sum(t.num(r, "room_classes") for r in rows), 1)
        facts["instRoomDomainMin"] = t.fmt(min(t.num(r, "room_domain_min") for r in rows), 0)
        facts["instRoomDomainMax"] = t.fmt(max(t.num(r, "room_domain_max") for r in rows), 0)
        facts["instSlotDomainMean"] = t.fmt(t.mean([t.num(r, "slot_domain_mean") for r in rows]), 0)
        for n in SIZES:
            sub = [r for r in rows if int(r["n"]) == n]
            for k in ("rooms", "groups", "lecturers"):
                facts[f"inst{n}{k.capitalize()}"] = t.fmt(round(t.mean([t.num(r, k) for r in sub])), 0)
            facts[f"inst{n}SportPct"] = pct(t.mean([t.num(r, "sport") for r in sub]))

    # --- Two-armed model -----------------------------------------------------------------------
    path = os.path.join(R, "twoarm.csv")
    if not t.complete(path) and not args.partial:
        raise Incomplete("twoarm.csv missing or still being written")
    if t.complete(path):
        rows = t.read(path)
        kappas = sorted({t.num(r, "kappa") for r in rows}, reverse=True)
        ratios = sorted({int(r["cost_ratio"]) for r in rows})
        reps = check_cells("twoarm", rows, ["kappa", "credit", "cost_ratio"],
                           [(k, c, q) for k in TWOARM_KAPPAS
                            for c in ("none", "application", "work") for q in TWOARM_RATIOS],
                           "rep", args.partial, TWOARM_REPS)
        by = collections.defaultdict(list)
        for r in rows:
            by[(t.num(r, "kappa"), r["credit"], int(r["cost_ratio"]))].append(r)
        facts["twoArmReps"] = str(reps)
        facts["twoArmKappaSolver"] = t.fmt(kappas[0], 0)
        facts["twoArmKappaLow"] = t.fmt(kappas[-1], 0)
        facts["twoArmBudget"] = thin(t.num(rows[0], "work_budget"))
        prow, maxerr, clipm = [], collections.defaultdict(float), collections.defaultdict(float)
        figdat = collections.defaultdict(list)
        for k in kappas:
            K = "Solver" if k == kappas[0] else "Low"
            for q in ratios:
                a, w, u = (by[(k, c, q)] for c in ("application", "work", "none"))
                cells = []
                for v in (a, w):
                    late = [t.num(r, "arm1_share_late") for r in v]
                    pred = t.mean([t.num(r, "predicted_clipped") for r in v])
                    cells += [t.fmt(t.mean(late), 3), t.fmt(pred, 3)]
                for credit, v in (("application", a), ("work", w)):
                    late = [t.num(r, "arm1_share_late") for r in v]
                    err = abs(t.mean(late) - t.mean([t.num(r, "predicted_clipped") for r in v]))
                    maxerr[(K, credit)] = max(maxerr[(K, credit)], err)
                    up = sum(t.num(r, "updates") for r in v)
                    cl = sum(t.num(r, "clipped_updates") for r in v)
                    clipm[(K, credit)] = max(clipm[(K, credit)], cl / up if up else 0)
                ea = [t.num(r, "exclusive") for r in a]
                ew = [t.num(r, "exclusive") for r in w]
                da = [t.num(r, "arm2_draws") for r in a]
                dw = [t.num(r, "arm2_draws") for r in w]
                dEx = 1 - t.mean(ew) / t.mean(ea)
                dDr = 1 - t.mean(dw) / t.mean(da)
                prow.append([t.fmt(k, 0) if q == ratios[0] else "", q] + cells
                            + [t.fmt(t.mean(ea), 1), t.fmt(t.mean(ew), 1), pct(dEx), pct(dDr)])
                T = f"twoArm{K}R{q}"
                facts[T + "ExApp"] = t.fmt(t.mean(ea), 1)
                facts[T + "ExWork"] = t.fmt(t.mean(ew), 1)
                lo, hi = t.mean_ci(ea)
                facts[T + "ExAppLo"], facts[T + "ExAppHi"] = t.fmt(lo, 1), t.fmt(hi, 1)
                lo, hi = t.mean_ci(ew)
                facts[T + "ExWorkLo"], facts[T + "ExWorkHi"] = t.fmt(lo, 1), t.fmt(hi, 1)
                facts[T + "ExLoss"] = pct(dEx)
                facts[T + "DrawLoss"] = pct(dDr)
                facts[T + "ShareApp"] = t.fmt(t.mean([t.num(r, "arm1_share_late") for r in a]), 3)
                facts[T + "ShareWork"] = t.fmt(t.mean([t.num(r, "arm1_share_late") for r in w]), 3)
                facts[T + "ShareAppWhole"] = t.fmt(t.mean([t.num(r, "arm1_share") for r in a]), 3)
                facts[T + "ShareWorkWhole"] = t.fmt(t.mean([t.num(r, "arm1_share") for r in w]), 3)
                # The paired change in exclusive proposals (common seeds across the modes).
                pa, pw = by_seed_rep(a), by_seed_rep(w)
                diffs = [pa[i] - pw[i] for i in sorted(set(pa) & set(pw))]
                lo, hi = t.mean_ci(diffs)
                facts[T + "ExDiff"] = t.fmt(t.mean(diffs), 1)
                facts[T + "ExDiffLo"] = t.fmt(lo, 1)
                facts[T + "ExDiffHi"] = t.fmt(hi, 1)
                if k == kappas[-1]:
                    for credit, v in (("none", u), ("application", a), ("work", w)):
                        late = [t.num(r, "arm1_share_late") for r in v]
                        ex = [t.num(r, "exclusive") for r in v]
                        s_lo, s_hi = t.mean_ci(late)
                        e_lo, e_hi = t.mean_ci(ex)
                        pred = t.mean([t.num(r, "predicted_clipped") for r in v])
                        figdat[credit].append([q, round(t.mean(late), 4), round(t.mean(late) - s_lo, 4),
                                               round(pred, 4), round(t.mean(ex), 2),
                                               round(t.mean(ex) - e_lo, 2)])
            facts[f"twoArmSegmentsMin{K}"] = t.fmt(min(
                t.num(r, "segments") for r in rows if t.num(r, "kappa") == k and r["credit"] != "none"), 0)
        facts["twoArmSegmentsMin"] = t.fmt(min(
            t.num(r, "segments") for r in rows if r["credit"] != "none"), 0)
        for k in kappas:
            K = "Solver" if k == kappas[0] else "Low"
            sub = [r for r in rows if t.num(r, "kappa") == k and r["credit"] != "none"]
            up = sum(t.num(r, "updates") for r in sub)
            facts[f"twoArmClipPooled{K}"] = pct(sum(t.num(r, "clipped_updates") for r in sub) / up)
        for (K, credit), e in maxerr.items():
            facts[f"twoArmErr{K}{credit}"] = t.fmt(e, 3)
        for K in ("Solver", "Low"):
            facts[f"twoArmErr{K}"] = t.fmt(max(maxerr[(K, c)] for c in ("application", "work")), 3)
            facts[f"twoArmClip{K}"] = pct(max(clipm[(K, c)] for c in ("application", "work")))
            for c in ("application", "work"):
                facts[f"twoArmClip{K}{c}"] = pct(clipm[(K, c)])
        nones = [t.num(r, "arm1_share") for r in rows if r["credit"] == "none"]
        facts["twoArmShareNoneMin"] = t.fmt(min(t.mean(v) for kk, v in group_means(rows, "none").items()), 3)
        facts["twoArmShareNoneMax"] = t.fmt(max(t.mean(v) for kk, v in group_means(rows, "none").items()), 3)
        facts["twoArmRatioMax"] = str(ratios[-1])
        head = (" & & \\multicolumn{2}{c}{application} & \\multicolumn{2}{c}{work} & "
                "\\multicolumn{3}{c}{exclusive proposals per run} & arm~2 \\\\\n"
                "\\cmidrule(lr){3-4}\\cmidrule(lr){5-6}\\cmidrule(lr){7-9}\n"
                "$\\kappa$ & $\\omega_2/\\omega_1$ & share & m.-f. & share & m.-f. & appl. & work "
                "& $-$\\% & $-$\\% \\\\\n")
        t.write_table(
            os.path.join(TAB, "twoarm.tex"),
            "The two-armed model: share of the cheap arm and exclusive proposals",
            "tab:twoarm", [], prow, align="rr" + "r" * 8, small=True, header_tex=head,
            note=("Means over %s repetitions of a $%s$-unit work budget; the three credit modes of a "
                  "repetition share its random-number seed. \\emph{share}: share of arm~1 among the "
                  "draws of the second half of the work; \\emph{m.-f.}: mean-field share "
                  "\\eqref{eq:mfshare} with the fixed point clipped to the weight bounds. "
                  "\\emph{exclusive}: proposals, counted with repetition, of states only arm~2 can "
                  "propose; $-$\\%%: reduction under \\emph{work} relative to \\emph{application}, "
                  "of exclusive proposals and of draws of arm~2.")
            % (reps, facts["twoArmBudget"]))
        for credit, d in figdat.items():
            t.write_dat(os.path.join(FIG, f"twoarm-{credit}.dat"),
                        ["ratio", "share", "share_err", "predicted", "exclusive", "exclusive_err"], d)

    # --- Fixed mixture: Theorem 1 --------------------------------------------------------------
    path = os.path.join(R, "fixedq.csv")
    if not t.complete(path) and not args.partial:
        raise Incomplete("fixedq.csv missing or still being written")
    if t.complete(path):
        rows = t.read(path)
        reps = check_cells("fixedq", rows, ["omega2", "q", "work_budget"],
                           [(o, q, w) for o in FIXEDQ_OMEGAS for q in FIXEDQ_QS
                            for w in FIXEDQ_BUDGETS], "rep", args.partial, FIXEDQ_REPS)
        facts["fixedqReps"] = thin(reps)
        by = collections.defaultdict(list)
        for r in rows:
            by[(t.num(r, "omega2"), t.num(r, "q"), t.num(r, "work_budget"))].append(r)
        inside, cells, rateDev, hitDev, hitCover, hitCells = 0, 0, 0.0, 0.0, 0, 0
        fluidErr = 0.0
        hits = {}
        for (o2, q, W), v in sorted(by.items()):
            A = q * 1 + (1 - q) * o2
            n2 = [t.num(r, "n2") for r in v]
            lo, hi = t.mean_ci(n2)
            lower = (1 - q) * W / A
            upper = (1 - q) * (W + max(o2, 1)) / A
            cells += 1
            if hi >= lower and lo <= upper:
                inside += 1
            if W >= 1e5 and q < 1:
                rateDev = max(rateDev, abs(t.mean(n2) / W / ((1 - q) / A) - 1))
            if W < 1e5 and q < 1:
                fluidErr = max(fluidErr, abs(t.mean(n2) / lower - 1))
            hit = [t.num(r, "firsthit_work") for r in v]
            theory = A / ((1 - q) * 0.01)
            hits.setdefault((o2, q), []).extend(hit)
        facts["fixedqHitObs"] = thin(min(len(v) for v in hits.values()))
        for (o2, q), hit in hits.items():
            A = q + (1 - q) * o2
            theory = A / ((1 - q) * 0.01)
            lo, hi = t.mean_ci(hit)
            hitCells += 1
            if lo <= theory <= hi:
                hitCover += 1
            hitDev = max(hitDev, abs(t.mean(hit) / theory - 1))
        facts["fixedqCells"] = str(cells)
        facts["fixedqInside"] = str(inside)
        facts["fixedqRateDevPct"] = t.fmt(100 * rateDev, 2)
        facts["fixedqFluidErrPct"] = t.fmt(100 * fluidErr, 1)
        facts["fixedqHitCells"] = str(hitCells)
        facts["fixedqHitCover"] = str(hitCover)
        facts["fixedqHitDevPct"] = t.fmt(100 * hitDev, 1)
        o2max = max(k[0] for k in hits)
        qmax = max(k[1] for k in hits)
        rho_th = 1 + qmax / (1 - qmax) / o2max
        rho_me = t.mean(hits[(o2max, qmax)]) / t.mean(hits[(o2max, 0.0)])
        facts["fixedqOmega"] = t.fmt(o2max, 0)
        facts["fixedqQ"] = t.fmt(qmax, 2)
        facts["fixedqRhoTheory"] = t.fmt(rho_th, 3)
        facts["fixedqRhoMeasured"] = t.fmt(rho_me, 3)
        # The small-budget example: the fluid value against the measured mean and the bounds.
        Ws = min(k[2] for k in by)
        ex = by[(o2max, 0.5, Ws)]
        A = 0.5 + 0.5 * o2max
        facts["fixedqSmallW"] = t.fmt(Ws, 0)
        facts["fixedqSmallMean"] = f'{t.mean([t.num(r, "n2") for r in ex]):.2f}'
        facts["fixedqSmallLower"] = t.fmt(0.5 * Ws / A, 2)
        facts["fixedqSmallUpper"] = t.fmt(0.5 * (Ws + o2max) / A, 2)

    t.write_macros(os.path.join(TAB, "facts.tex"), facts)
    print(f"wrote {len(facts)} macros")


def by_seed_rep(rows):
    return {r["rep"]: t.num(r, "exclusive") for r in rows}


def group_means(rows, credit):
    out = collections.defaultdict(list)
    for r in rows:
        if r["credit"] == credit:
            out[(r["kappa"], r["cost_ratio"])].append(t.num(r, "arm1_share_late"))
    return out


if __name__ == "__main__":
    try:
        main()
    except Incomplete as e:
        print(f"analysis stopped: {e}")
        sys.exit(1)
