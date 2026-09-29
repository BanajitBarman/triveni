// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Banajit Barman

/// triveni-pp: inelastic pp events with the configured proton shape.
///
/// Calibrates sigma_gg to PP:sigmaInel, generates Main:numberOfEvents inelastic
/// events and writes
///   <stem>.txt         one row per event: b, T_pp, P_coll, entropy S, eps2-4, string lengths
///   <stem>_maps.txt    entropy density s(x) of the first PP:nMaps events
///   <stem>_quarks.txt  quark and junction positions of those events
///   <stem>.log         settings, calibration and the independent sigma_inel check
/// Rows are ordered by random-number block, so output is independent of Main:threads.

#include "triveni/PPCollision.h"
#include "triveni/Parallel.h"
#include "triveni/Random.h"
#include "triveni/Settings.h"
#include "triveni/Version.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <vector>

namespace
{

constexpr double kPi = 3.14159265358979323846;
constexpr int kBlock = 200;

struct Row {
  double b, overlap, pCollision, entropy, eps2, eps3, eps4, length1, length2;
};

struct Snapshot {
  triveni::Grid entropy;
  triveni::PPEvent event;
};

std::ofstream openOutput(const std::string& path)
{
  std::ofstream out(path);
  if (!out) {
    throw std::runtime_error("cannot write '" + path + "'");
  }
  out.precision(7);
  return out;
}

} // namespace

int main(int argc, char** argv)
{
  try {
    triveni::Settings settings;
    if (!settings.readArguments(argc, argv, std::cout)) {
      return 0;
    }
    const auto start = std::chrono::steady_clock::now();
    const triveni::PPParams params = triveni::PPParams::fromSettings(settings);
    triveni::PPCollider collider(params);
    const std::string stem = settings.word("Main:output");
    const std::uint64_t seed = static_cast<std::uint64_t>(settings.mode("Random:seed"));
    const int nEvents = settings.mode("Main:numberOfEvents");
    const int nThreads = settings.mode("Main:threads");
    const int nMaps = settings.mode("PP:nMaps");

    const triveni::Calibration cal =
      collider.calibrate(triveni::Random::streamSeed(seed, 1), settings.mode("PP:calibrationPairs"), nThreads);

    // Event generation: each block produces a fixed number of accepted events from its own stream.
    const std::uint64_t eventSeed = triveni::Random::streamSeed(seed, 2);
    const int nBlocks = (nEvents + kBlock - 1) / kBlock;
    std::vector<std::vector<Row>> blockRows(nBlocks);
    std::vector<long long> blockTrials(nBlocks, 0);
    std::vector<Snapshot> snapshots;

    triveni::parallelBlocks(nBlocks, nThreads, [&](int block) {
      triveni::Random rng(triveni::Random::streamSeed(eventSeed, static_cast<std::uint64_t>(block)));
      triveni::PPCollider::Workspace ws(params);
      triveni::PPEvent event;
      const int want = std::min(kBlock, nEvents - block * kBlock);
      auto& rows = blockRows[block];
      rows.reserve(want);
      while (static_cast<int>(rows.size()) < want) {
        ++blockTrials[block];
        if (!collider.tryEvent(rng, ws, event)) {
          continue;
        }
        rows.push_back({event.b, event.overlap, event.pCollision, event.entropy, event.ecc.eps[2], event.ecc.eps[3],
                        event.ecc.eps[4], event.p1.stringLength, event.p2.stringLength});
        if (block == 0 && static_cast<int>(snapshots.size()) < nMaps) {
          snapshots.push_back({ws.t1, event});
        }
      }
    });

    long long trials = 0;
    {
      auto out = openOutput(stem + ".txt");
      out << "b/D:overlap/D:pColl/D:S/D:eps2/D:eps3/D:eps4/D:L1/D:L2/D\n";
      for (int block = 0; block < nBlocks; ++block) {
        trials += blockTrials[block];
        for (const Row& r : blockRows[block]) {
          out << r.b << ' ' << r.overlap << ' ' << r.pCollision << ' ' << r.entropy << ' ' << r.eps2 << ' ' << r.eps3
              << ' ' << r.eps4 << ' ' << r.length1 << ' ' << r.length2 << '\n';
        }
      }
    }
    {
      auto maps = openOutput(stem + "_maps.txt");
      auto quarks = openOutput(stem + "_quarks.txt");
      maps << "event/I:x/D:y/D:s/D\n";
      quarks << "event/I:proton/I:item/I:x/D:y/D\n";
      for (std::size_t e = 0; e < snapshots.size(); ++e) {
        const triveni::Grid& g = snapshots[e].entropy;
        for (int iy = 0; iy < g.size(); ++iy) {
          for (int ix = 0; ix < g.size(); ++ix) {
            maps << e << ' ' << g.coord(ix) << ' ' << g.coord(iy) << ' ' << g(ix, iy) << '\n';
          }
        }
        const triveni::PPEvent& ev = snapshots[e].event;
        const triveni::Proton* protons[2] = {&ev.p1, &ev.p2};
        for (int k = 0; k < 2; ++k) {
          const double shift = (k == 0 ? -0.5 : 0.5) * ev.b;
          for (int item = 0; item < 4; ++item) {
            const triveni::Vec3 v = item < 3 ? protons[k]->quarks[item] : protons[k]->junction;
            quarks << e << ' ' << k << ' ' << item << ' ' << v.x + shift << ' ' << v.y << '\n';
          }
        }
      }
    }

    // Independent check: the accepted fraction measures sigma_inel with fresh random numbers.
    const double area = kPi * params.bMax * params.bMax;
    const double fraction = static_cast<double>(nEvents) / static_cast<double>(trials);
    const double sigmaCheck = 10.0 * area * fraction;
    const double sigmaCheckError = 10.0 * area * std::sqrt(fraction * (1.0 - fraction) / static_cast<double>(trials));
    const double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();

    std::ostringstream summary;
    summary << "TRIVENI " << triveni::kVersion << '\n'
            << "triveni-pp: shape " << settings.word("Proton:shape") << ", " << nEvents << " inelastic events from "
            << trials << " trials in " << seconds << " s\n"
            << "  sigma_gg                = " << cal.sigmaGG << " fm^2 (" << cal.pairs << " calibration pairs)\n"
            << "  sigma_inel (calibrated) = " << cal.sigmaInel << " +- " << cal.sigmaInelError << " mb (target "
            << params.sigmaInel << ")\n"
            << "  sigma_inel (check)      = " << sigmaCheck << " +- " << sigmaCheckError << " mb\n"
            << "  <P_coll> for b > 0.9 bMax = " << cal.edgeProbability
            << (cal.edgeProbability > 1e-3 ? "  WARNING: increase PP:bMax" : "") << '\n'
            << "  dln(sigma_inel)/dln(sigma_gg) = " << cal.sensitivity
            << (cal.sensitivity < triveni::Calibration::kBlackDiskSensitivity
                  ? "  WARNING: black-disk regime, sigma_gg ill-determined; the proton is too small for"
                    " this sigma_inel (increase Proton:quarkSpread)"
                  : "")
            << '\n';
    std::cout << summary.str() << "  wrote " << stem << "{.txt,_maps.txt,_quarks.txt,.log}\n";

    auto log = openOutput(stem + ".log");
    log << summary.str() << "\nSettings (* = changed):\n";
    settings.list(log);
    return 0;
  } catch (const std::exception& e) {
    std::cerr << "triveni-pp: " << e.what() << '\n';
    return 1;
  }
}
