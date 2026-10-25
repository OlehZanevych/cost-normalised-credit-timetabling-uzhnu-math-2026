// Experiments for "Cost-normalised credit can favour cheap operators".
//
//   ablation   leave one operator out, one at a time                       (Table 2)
//   covered    add winfix (covered by move) or dayfix (aimed, coverage not
//              established) under per-application and cost-normalised credit (Table 4)
//   credit     five credit regimes: none, per application, per work, per
//              work in draw equivalents, per microsecond                    (Table 3)
//   twoarm     the allocation effect in a two-armed model without timetabling (Table 1, Fig. 1)
//   fixedq     Theorem 1 at a fixed mixture: draws bought by a work budget, and work to first hit
//   instancestats  the realised composition of the generated instances (Section 7)
//
// Every solver run uses one worker and a fixed budget of work (bucket recomputations).
// `instances` writes the instance files; `selftest` checks the generator, the file format and the
// two evaluators against each other.
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <numeric>
#include <string>
#include <vector>

#include "core/generate.hpp"
#include "core/io.hpp"
#include "core/rng.hpp"
#include "core/search.hpp"
#include "core/selftest.hpp"
#include "core/validate.hpp"

using namespace tt;

namespace {

std::string gResults = "results";
// The part of a solver design to run: by default both sizes and instances 1..seeds.  Splitting a
// design (--sizes, --seedfrom) runs exactly the same records, since each run depends only on its own
// size, instance seed and search seed; the parts are concatenated afterwards.
std::vector<int> gSizes = {800, 3200};
int gSeedFrom = 1;

std::ofstream openCsv(const std::string& name, const std::string& header) {
  // A missing results directory would otherwise make every record vanish silently.
  std::error_code ec;
  std::filesystem::create_directories(gResults, ec);
  const std::string path = gResults + "/" + name;
  std::ofstream f(path);
  if (!f) {
    std::fprintf(stderr, "cannot open %s for writing\n", path.c_str());
    std::exit(1);
  }
  // Unbuffered: one run is one line and runs are seconds apart, so flushing costs nothing and a
  // long experiment can be watched, or salvaged, while it is still running.
  f << std::unitbuf;
  f << header << "\n";
  return f;
}

/// Stop with an error if any record of an experiment failed to reach the disk.
void finishCsv(std::ofstream& f, const char* name) {
  f.flush();
  if (!f) {
    std::fprintf(stderr, "writing %s/%s failed\n", gResults.c_str(), name);
    std::exit(1);
  }
  std::printf("%s done\n", name);
}

GeneratedInstance makeInstance(int n, uint64_t seed) {
  GenOptions go;
  go.nClasses = n;
  go.seed = seed;
  return generateInstance(go);
}

void writeOperatorColumns(std::ofstream& f, const SearchResult& r) {
  for (int i = 0; i < kNumOperators; ++i)
    f << "," << r.ops[i].draws << "," << r.ops[i].applied << "," << r.ops[i].improved << ","
      << r.ops[i].work << "," << r.ops[i].gain << "," << r.ops[i].finalWeight;
}

std::string operatorHeader() {
  std::string h;
  for (int i = 0; i < kNumOperators; ++i) {
    const std::string n = operatorName(i);
    h += "," + n + "_draws," + n + "_applied," + n + "_improved," + n + "_work," + n + "_gain," +
         n + "_weight";
  }
  return h;
}

// =================================================================================================
//
// Every solver experiment runs ONE worker (so a run is a deterministic function of its seeds), with
// elite adoption off, under a fixed budget of WORK -- bucket recomputations, counted by the state
// itself and therefore including construction, the escape phases and every draw that produced no
// candidate.  A run stops at the first check after the budget is reached, so it may overshoot by the
// cost of one draw or one escape phase; the work actually spent is recorded beside the budget.

const char* kRunHeader =
    "n,seed,%s,hard,soft,objective,unplaced,clean,work_budget,work,candidates,draws,elapsed_ms,"
    "deep,kicks,fresh,adoptions,weight_updates,ceiling_hits,floor_hits";

std::string runHeader(const char* key) {
  char buf[512];
  std::snprintf(buf, sizeof(buf), kRunHeader, key);
  return std::string(buf) + operatorHeader();
}

SearchOptions baseOptions(long long workBudget, uint64_t seed) {
  SearchOptions so;
  so.timeLimitMs = 4LL * 3600 * 1000;   // a safety net only: every recorded run ends on its work budget
  so.candidateLimit = 0;
  so.workLimit = workBudget;
  so.threads = 1;
  so.cooperate = false;
  so.seed = seed;
  return so;
}

void writeRun(std::ofstream& f, int n, int seed, const std::string& key, const GeneratedInstance& gi,
              const SearchResult& r, long long workBudget) {
  Score sc = validate(gi.instance, r.best);
  f << n << "," << seed << "," << key << "," << sc.hard() << "," << sc.soft() << ","
    << sc.objective(kDefaultBeta) << "," << sc.unplaced << "," << (sc.clean() ? 1 : 0) << ","
    << workBudget << "," << r.workUnits << "," << r.candidates << "," << r.draws << ","
    << r.elapsedMs << "," << r.deepPhases << "," << r.kicks << "," << r.freshRestarts << ","
    << r.adoptions << "," << r.weightUpdates << "," << r.ceilingHits << "," << r.floorHits;
  writeOperatorColumns(f, r);
  f << "\n";
}

void expAblation(int seeds, long long budget) {
  auto f = openCsv("ablation.csv", runHeader("arm"));
  struct Arm { const char* name; void (*apply)(SearchOptions&); };
  const Arm arms[] = {
      {"full", [](SearchOptions&) {}},
      {"no-move", [](SearchOptions& o) { o.useMove = false; }},
      {"no-swap", [](SearchOptions& o) { o.useSwap = false; }},
      {"no-chain", [](SearchOptions& o) { o.useChain = false; }},
      {"no-kempe", [](SearchOptions& o) { o.useKempe = false; }},
      {"no-ruin", [](SearchOptions& o) { o.useLns = false; }},
      {"no-repack", [](SearchOptions& o) { o.useRepack = false; }},
      {"no-kopt", [](SearchOptions& o) { o.useKopt = false; }},
      // Not an operator removal: ruin keeps running, without its community selector.
      {"no-cluster", [](SearchOptions& o) { o.useClusters = false; }},
      // Only move and swap in the adaptive loop; the escape phases are unchanged.
      {"simple", [](SearchOptions& o) {
         o.useChain = o.useKempe = o.useLns = o.useRepack = o.useKopt = false;
       }},
  };
  for (int n : gSizes) {
    for (int seed = gSeedFrom; seed <= seeds; ++seed) {
      GeneratedInstance gi = makeInstance(n, seed);
      for (const Arm& a : arms) {
        SearchOptions so = baseOptions(budget, 20260912 + seed);
        a.apply(so);
        writeRun(f, n, seed, a.name, gi, solve(gi.instance, so), budget);
      }
    }
  }
  finishCsv(f, "ablation.csv");
}

void expCovered(int seeds, long long budget) {
  auto f = openCsv("covered.csv", runHeader("credit,arm"));
  struct Arm { const char* name; bool dayfix; bool winfix; };
  const Arm arms[] = {{"base", false, false}, {"+winfix", false, true}, {"+dayfix", true, false}};
  const Credit credits[] = {Credit::Application, Credit::Work};
  const char* creditName[] = {"application", "work"};
  for (int n : gSizes) {
    for (int seed = gSeedFrom; seed <= seeds; ++seed) {
      GeneratedInstance gi = makeInstance(n, seed);
      for (size_t ci = 0; ci < 2; ++ci) {
        for (const Arm& a : arms) {
          SearchOptions so = baseOptions(budget, 777 + seed);
          so.useDayFix = a.dayfix;
          so.useWinFix = a.winfix;
          so.credit = credits[ci];
          writeRun(f, n, seed, std::string(creditName[ci]) + "," + a.name, gi,
                   solve(gi.instance, so), budget);
        }
      }
    }
  }
  finishCsv(f, "covered.csv");
}

/// The reference cost that makes the work normaliser dimensionless: 30 buckets, about the mean work
/// of one draw of the default portfolio.  Declared, not fitted.
constexpr double kReferenceCost = 30.0;

void expCredit(int seeds, long long budget) {
  auto f = openCsv("credit.csv", runHeader("credit"));
  struct Regime { const char* name; Credit credit; double unit; };
  const Regime regimes[] = {
      {"none", Credit::None, 1.0},
      {"application", Credit::Application, 1.0},
      {"work", Credit::Work, 1.0},                     // raw bucket count, kappa unchanged
      {"work-scaled", Credit::Work, kReferenceCost},   // bucket count / 30: same scale as draws
      {"time", Credit::Time, 1.0},                     // microseconds, kappa unchanged
  };
  for (int n : gSizes) {
    for (int seed = gSeedFrom; seed <= seeds; ++seed) {
      GeneratedInstance gi = makeInstance(n, seed);
      for (const Regime& g : regimes) {
        SearchOptions so = baseOptions(budget, 5150 + seed);
        so.credit = g.credit;
        so.creditUnit = g.unit;
        writeRun(f, n, seed, g.name, gi, solve(gi.instance, so), budget);
      }
    }
  }
  finishCsv(f, "credit.csv");
}

/// A two-armed model with the timetabling removed.
///
/// Arm 1 costs one unit of work and arm 2 costs omega_2.  Each arm earns the reward 1.5 with a fixed
/// probability per draw (0.12 and 0.06); a draw of arm 2 additionally produces, with probability
/// 0.01, an "exclusive" event -- a proposal of a state arm 1 cannot propose.  The exclusive events
/// do not feed the reward: the model shows how the proxy reward and the cost normaliser allocate a
/// work budget, not which allocation maximises the reward itself.  The weight update is the
/// solver's.  All three credit modes of one repetition share their random-number stream.
void expTwoArm(int reps) {
  auto f = openCsv("twoarm.csv",
                   "kappa,credit,cost_ratio,rep,work_budget,work,segments,arm1_share,"
                   "arm1_share_late,arm2_draws,exclusive,clipped_updates,updates,"
                   "predicted_arm1_share,predicted_clipped");
  const char* names[] = {"none", "application", "work"};
  for (double kappa : {300.0, 60.0}) {
    const double ceiling = 40.0, floorW = 0.05;
    for (int ratio : {2, 5, 10, 20, 40}) {
      for (int mode = 0; mode < 3; ++mode) {
        for (int rep = 1; rep <= reps; ++rep) {
          Rng rng(9001 + rep * 31 + ratio * 7 + static_cast<int>(kappa));
          const double c1 = 1.0, c2 = ratio;
          const double p1 = 0.12, p2 = 0.06, p2Exclusive = 0.010;
          const double workBudget = 400000;
          double w1 = 4.0, w2 = 1.0;
          double work = 0, seg1 = 0, seg2 = 0, sw1 = 0, sw2 = 0, d1 = 0, d2 = 0;
          long long draws1 = 0, draws2 = 0, late1 = 0, late2 = 0, exclusive = 0;
          long long segDraws = 0, segments = 0, clipped = 0, updates = 0;
          while (work < workBudget) {
            const bool pickOne = (mode == 0) ? rng.coin(0.5) : rng.coin(w1 / (w1 + w2));
            const bool late = work >= 0.5 * workBudget;   // second half: after the burn-in
            if (pickOne) {
              ++draws1;
              if (late) ++late1;
              work += c1;
              sw1 += c1;
              d1 += 1;
              if (rng.coin(p1)) seg1 += 1.5;
            } else {
              ++draws2;
              if (late) ++late2;
              work += c2;
              sw2 += c2;
              d2 += 1;
              if (rng.coin(p2)) seg2 += 1.5;
              if (rng.coin(p2Exclusive)) ++exclusive;
            }
            if (++segDraws >= 2000) {
              ++segments;
              if (mode != 0) {
                const double n1 = (mode == 1) ? std::max(1.0, d1) : std::max(1.0, sw1);
                const double n2 = (mode == 1) ? std::max(1.0, d2) : std::max(1.0, sw2);
                w1 = 0.7 * w1 + 0.3 * (1 + kappa * seg1 / n1);
                w2 = 0.7 * w2 + 0.3 * (1 + kappa * seg2 / n2);
                for (double* w : {&w1, &w2}) {
                  ++updates;
                  if (*w >= ceiling || *w <= floorW) ++clipped;
                  *w = std::min(ceiling, std::max(floorW, *w));
                }
              }
              seg1 = seg2 = sw1 = sw2 = d1 = d2 = 0;
              segDraws = 0;
            }
          }
          // The mean-field share (fixed credit rates, deterministic update) and the same share
          // with the weight bounds applied to the fixed point.
          const double rate1 = (mode == 2) ? p1 * 1.5 / c1 : p1 * 1.5;
          const double rate2 = (mode == 2) ? p2 * 1.5 / c2 : p2 * 1.5;
          const double f1 = 1 + kappa * rate1, f2 = 1 + kappa * rate2;
          const double c1w = std::min(ceiling, std::max(floorW, f1));
          const double c2w = std::min(ceiling, std::max(floorW, f2));
          const double pred = (mode == 0) ? 0.5 : f1 / (f1 + f2);
          const double predClipped = (mode == 0) ? 0.5 : c1w / (c1w + c2w);
          f << kappa << "," << names[mode] << "," << ratio << "," << rep << "," << workBudget << ","
            << work << "," << segments << ","
            << double(draws1) / std::max(1LL, draws1 + draws2) << ","
            << double(late1) / std::max(1LL, late1 + late2) << "," << draws2 << "," << exclusive
            << "," << clipped << "," << updates << "," << pred << "," << predClipped << "\n";
        }
      }
    }
  }
  finishCsv(f, "twoarm.csv");
}

/// Theorem 1 checked directly: a FIXED mixture q of a unit-cost operator and an operator of cost
/// omega_2, with no adaptation.  Part (a): the number of expensive draws a work budget W buys, with
/// a draw executed whenever the work so far is below W (so the last draw may overshoot).  Part (b):
/// the work spent until the first exclusive proposal, each expensive draw producing one with
/// probability p independently of the past, with no budget (the run continues until it occurs).
void expFixedQ(int reps) {
  auto f = openCsv("fixedq.csv", "omega2,q,rep,work_budget,n2,work,firsthit_work");
  const double p = 0.01;
  for (int omega2 : {10, 40}) {
    for (double q : {0.0, 0.25, 0.5, 0.75, 0.9}) {
      for (double W : {100.0, 100000.0}) {
        for (int rep = 1; rep <= reps; ++rep) {
          Rng rng(424242 + rep * 101 + omega2 * 7 + static_cast<int>(q * 1000) +
                  static_cast<int>(W));
          double work = 0;
          long long n2 = 0;
          while (work < W) {
            if (rng.coin(q)) {
              work += 1.0;
            } else {
              work += omega2;
              ++n2;
            }
          }
          double hit = 0;
          for (;;) {
            if (rng.coin(q)) {
              hit += 1.0;
            } else {
              hit += omega2;
              if (rng.coin(p)) break;
            }
          }
          f << omega2 << "," << q << "," << rep << "," << W << "," << n2 << "," << work << ","
            << hit << "\n";
        }
      }
    }
  }
  finishCsv(f, "fixedq.csv");
}

/// What the generator actually produced, instance by instance: the realised proportions of the
/// class kinds and the sizes of the domains, for the article's description of the instances (its
/// generator settings are probabilities and sizing factors, not exact proportions).
void expInstanceStats(int seeds) {
  auto f = openCsv("instances.csv",
                   "n,seed,classes,lecturers,groups,rooms,shared_places,buildings,biweekly,online,"
                   "sport,immovable,multigroup,room_classes,room_domain_mean,room_domain_min,"
                   "room_domain_max,slot_domain_mean");
  for (int n : gSizes) {
    for (int seed = gSeedFrom; seed <= seeds; ++seed) {
      GeneratedInstance gi = makeInstance(n, seed);
      const Instance& in = gi.instance;
      int biweekly = 0, online = 0, sport = 0, immovable = 0, multi = 0, roomCls = 0;
      long long roomSum = 0, slotSum = 0;
      int roomMin = 1 << 30, roomMax = 0;
      for (int c = 0; c < in.nClasses(); ++c) {
        const ClassSpec& cs = in.classes[c];
        if (gi.planted[c].parity != Parity::Weekly) ++biweekly;
        if (cs.kind == PlaceKind::Online) ++online;
        if (cs.abstractRoom >= 0) ++sport;
        if (!cs.movable) ++immovable;
        if (cs.groups.size() > 1) ++multi;
        slotSum += static_cast<long long>(cs.slots.size());
        if (!cs.rooms.empty()) {
          ++roomCls;
          const int k = static_cast<int>(cs.rooms.size());
          roomSum += k;
          roomMin = std::min(roomMin, k);
          roomMax = std::max(roomMax, k);
        }
      }
      const double N = std::max(1, in.nClasses());
      f << n << "," << seed << "," << in.nClasses() << "," << in.nLecturers << "," << in.nGroups
        << "," << in.nRooms << "," << in.nAbstract << "," << in.nBuildings << "," << biweekly / N
        << "," << online / N << "," << sport / N << "," << immovable / N << "," << multi / N << ","
        << roomCls << "," << double(roomSum) / std::max(1, roomCls) << "," << roomMin << ","
        << roomMax << "," << double(slotSum) / N << "\n";
    }
  }
  finishCsv(f, "instances.csv");
}

/// Writes the instance files the four experiments above run on, so that the study can be repeated
/// without trusting this program's generator: the files are plain text and any solver can read them.
void writeInstances(const std::string& dir) {
  std::error_code ec;
  std::filesystem::create_directories(dir, ec);
  for (int n : {200, 400, 800, 1600, 3200}) {
    for (int seed = 1; seed <= 10; ++seed) {
      GeneratedInstance gi = makeInstance(n, seed);
      char path[512];
      std::snprintf(path, sizeof(path), "%s/n%05d-s%d.txt", dir.c_str(), n, seed);
      if (!writeInstance(path, gi.instance)) {
        std::fprintf(stderr, "cannot write %s\n", path);
        std::exit(1);
      }
      std::snprintf(path, sizeof(path), "%s/n%05d-s%d.planted.txt", dir.c_str(), n, seed);
      if (!writeSchedule(path, gi.planted)) {
        std::fprintf(stderr, "cannot write %s\n", path);
        std::exit(1);
      }
    }
  }
  std::printf("instances written to %s\n", dir.c_str());
}

}  // namespace

int main(int argc, char** argv) {
  const char* usage =
      "usage: experiments <ablation|covered|credit|twoarm|fixedq|instancestats|instances|selftest>\n"
      "                   [--seeds N] [--seedfrom S] [--sizes 800,3200] [--budget K]\n"
      "                   [--results DIR] [--data DIR]\n"
      "Without --seeds and --budget every subcommand runs the article's design: 10 instances per size\n"
      "and 10 000 000 bucket recomputations per solver run, 50 repetitions per two-armed cell, 2 000\n"
      "per fixed-mixture cell.\n";
  if (argc < 2) {
    std::printf("%s", usage);
    return 2;
  }
  const std::string cmd = argv[1];
  int seeds = -1;                // -1: the article's design for this subcommand
  long long budget = 10000000;   // bucket recomputations per solver run
  std::string data = "data";
  bool resultsGiven = false;
  for (int i = 2; i < argc; ++i) {
    if (!std::strcmp(argv[i], "--seeds") && i + 1 < argc) seeds = std::atoi(argv[++i]);
    else if (!std::strcmp(argv[i], "--budget") && i + 1 < argc) budget = std::atoll(argv[++i]);
    else if (!std::strcmp(argv[i], "--seedfrom") && i + 1 < argc) gSeedFrom = std::atoi(argv[++i]);
    else if (!std::strcmp(argv[i], "--sizes") && i + 1 < argc) {
      gSizes.clear();
      for (const char* p = argv[++i]; *p;) {
        gSizes.push_back(std::atoi(p));
        while (*p && *p != ',') ++p;
        if (*p == ',') ++p;
      }
    }
    else if (!std::strcmp(argv[i], "--results") && i + 1 < argc) {
      gResults = argv[++i];
      resultsGiven = true;
    } else if (!std::strcmp(argv[i], "--data") && i + 1 < argc) data = argv[++i];
    else {
      // A mistyped option must not silently run a different design.
      std::fprintf(stderr, "unknown or incomplete option %s\n%s", argv[i], usage);
      return 2;
    }
  }
  bool sizesOk = !gSizes.empty();
  for (int n : gSizes) sizesOk = sizesOk && n > 0;
  if (seeds == 0 || seeds < -1 || budget <= 0 || gSeedFrom < 1 || !sizesOk) {
    std::fprintf(stderr, "--seeds, --seedfrom, --sizes and --budget must be positive\n");
    return 2;
  }
  auto design = [&](int paper) { return seeds > 0 ? seeds : paper; };
  if (cmd == "ablation") expAblation(design(10), budget);
  else if (cmd == "covered") expCovered(design(10), budget);
  else if (cmd == "credit") expCredit(design(10), budget);
  else if (cmd == "twoarm") expTwoArm(design(50));
  else if (cmd == "fixedq") expFixedQ(design(2000));
  else if (cmd == "instancestats") expInstanceStats(design(10));
  else if (cmd == "instances") writeInstances(data);
  else if (cmd == "selftest") {
    // The self-test writes two scratch files; keep them out of results/ unless a directory is given.
    const std::string dir = resultsGiven
        ? gResults
        : (std::filesystem::temp_directory_path() / "tt-a3-selftest").string();
    return runSelfTest(dir);
  } else {
    std::printf("unknown subcommand %s\n%s", cmd.c_str(), usage);
    return 2;
  }
  return 0;
}
