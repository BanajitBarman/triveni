// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Banajit Barman

/// triveni-proton: single-proton geometry for the hot-spot, Y and Delta shapes.
///
/// Every sampled quark configuration is built with all three shapes, so the
/// comparison is like for like. Writes
///   <stem>_maps.txt     thickness maps of the first ProtonApp:nMaps configurations
///   <stem>_quarks.txt   quark and junction positions of those configurations
///   <stem>_stats.txt    per configuration: string lengths, junction position, eps2/eps3 per shape
///   <stem>_profile.txt  mean radial profile <T_p>(r) per shape (+ analytic hot-spot curve)
/// All files are whitespace-separated with a ROOT TTree::ReadFile header line.

#include "triveni/Grid.h"
#include "triveni/Parallel.h"
#include "triveni/Proton.h"
#include "triveni/Random.h"
#include "triveni/Settings.h"
#include "triveni/Version.h"
#include "triveni/Steiner.h"

#include <array>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <vector>

namespace
{

constexpr double kPi = 3.14159265358979323846;
constexpr std::array<triveni::ProtonShape, 3> kShapes = {triveni::ProtonShape::HotSpot, triveni::ProtonShape::Y,
                                                     triveni::ProtonShape::Delta};
constexpr int kBlock = 250;

std::ofstream openOutput(const std::string& path)
{
  std::ofstream out(path);
  if (!out) {
    throw std::runtime_error("cannot write '" + path + "'");
  }
  out.precision(7);
  return out;
}

struct StatsRow {
  double steinerLength = 0.0;
  double perimeter = 0.0;
  double junctionOffset = 0.0;
  int atQuark = 0;
  std::array<double, 3> eps2{};
  std::array<double, 3> eps3{};
};

} // namespace

int main(int argc, char** argv)
{
  try {
    triveni::Settings settings;
    if (!settings.readArguments(argc, argv, std::cout)) {
      return 0;
    }
    const auto start = std::chrono::steady_clock::now();
    const triveni::ProtonParams params = triveni::ProtonParams::fromSettings(settings);
    const triveni::ProtonBuilder builder(params);
    const std::string stem = settings.word("Main:output");
    const std::uint64_t seed = static_cast<std::uint64_t>(settings.mode("Random:seed"));
    const int nStats = settings.mode("Main:numberOfEvents");
    const int nMaps = settings.mode("ProtonApp:nMaps");
    const int nThreads = settings.mode("Main:threads");
    const double half = settings.parm("ProtonApp:mapHalfWidth");
    const double step = settings.parm("ProtonApp:mapStep");

    // 1) Thickness maps of a few configurations.
    {
      triveni::Random rng(triveni::Random::streamSeed(seed, 0));
      triveni::Grid grid(half, step);
      std::array<triveni::Grid, 3> maps{grid, grid, grid};
      auto mapsOut = openOutput(stem + "_maps.txt");
      auto quarksOut = openOutput(stem + "_quarks.txt");
      mapsOut << "config/I:x/D:y/D:hotspot/D:Y/D:Delta/D\n";
      quarksOut << "config/I:shape/I:item/I:x/D:y/D:z/D:atQuark/I\n";
      for (int c = 0; c < nMaps; ++c) {
        const auto quarks = builder.sampleQuarks(rng);
        for (std::size_t k = 0; k < kShapes.size(); ++k) {
          const triveni::Proton p = builder.build(quarks, kShapes[k], rng);
          maps[k].clear();
          maps[k].addProton(p, {0.0, 0.0}, false);
          for (int item = 0; item < 4; ++item) {
            const triveni::Vec3 v = item < 3 ? p.quarks[item] : p.junction;
            quarksOut << c << ' ' << k << ' ' << item << ' ' << v.x << ' ' << v.y << ' ' << v.z << ' '
                      << (p.junctionOnQuark >= 0 ? 1 : 0) << '\n';
          }
        }
        for (int iy = 0; iy < grid.size(); ++iy) {
          for (int ix = 0; ix < grid.size(); ++ix) {
            mapsOut << c << ' ' << grid.coord(ix) << ' ' << grid.coord(iy) << ' ' << maps[0](ix, iy) << ' '
                    << maps[1](ix, iy) << ' ' << maps[2](ix, iy) << '\n';
          }
        }
      }
    }

    // 2) Statistics and the mean radial profile, in parallel blocks.
    const triveni::Grid layout(half, step);
    const int nBins = static_cast<int>(half / (2.0 * layout.step()));
    const double binWidth = half / nBins;
    std::vector<StatsRow> rows(nStats);
    const int nBlocks = (nStats + kBlock - 1) / kBlock;
    std::vector<std::array<std::vector<double>, 3>> blockProfiles(nBlocks);
    std::vector<double> binArea(nBins, 0.0); // grid area per radial bin (same for every block)
    for (int iy = 0; iy < layout.size(); ++iy) {
      for (int ix = 0; ix < layout.size(); ++ix) {
        const int bin = static_cast<int>(std::hypot(layout.coord(ix), layout.coord(iy)) / binWidth);
        if (bin < nBins) {
          binArea[bin] += layout.cellArea();
        }
      }
    }

    triveni::parallelBlocks(nBlocks, nThreads, [&](int block) {
      triveni::Random rng(triveni::Random::streamSeed(seed, 1 + static_cast<std::uint64_t>(block)));
      triveni::Grid grid(half, step);
      auto& profile = blockProfiles[block];
      for (auto& p : profile) {
        p.assign(nBins, 0.0);
      }
      const int end = std::min(nStats, (block + 1) * kBlock);
      for (int i = block * kBlock; i < end; ++i) {
        const auto quarks = builder.sampleQuarks(rng);
        const triveni::FermatPoint junction = triveni::fermatPoint(quarks);
        const triveni::Vec3 centroid = (1.0 / 3.0) * (quarks[0] + quarks[1] + quarks[2]);
        StatsRow& row = rows[i];
        row.steinerLength = triveni::legLength(quarks, junction.point);
        row.perimeter =
          triveni::norm(quarks[0] - quarks[1]) + triveni::norm(quarks[1] - quarks[2]) + triveni::norm(quarks[2] - quarks[0]);
        row.junctionOffset = triveni::norm(junction.point - centroid);
        row.atQuark = junction.vertex >= 0 ? 1 : 0;
        for (std::size_t k = 0; k < kShapes.size(); ++k) {
          const triveni::Proton p = builder.build(quarks, kShapes[k], rng);
          grid.clear();
          grid.addProton(p, {0.0, 0.0}, false);
          const triveni::Eccentricities e = triveni::eccentricities(grid);
          row.eps2[k] = e.eps[2];
          row.eps3[k] = e.eps[3];
          for (int iy = 0; iy < grid.size(); ++iy) {
            for (int ix = 0; ix < grid.size(); ++ix) {
              const int bin = static_cast<int>(std::hypot(grid.coord(ix), grid.coord(iy)) / binWidth);
              if (bin < nBins) {
                profile[k][bin] += grid(ix, iy) * grid.cellArea();
              }
            }
          }
        }
      }
    });

    {
      auto out = openOutput(stem + "_stats.txt");
      out << "Lsteiner/D:perimeter/D:junctionOffset/D:atQuark/I:"
             "eps2_hotspot/D:eps3_hotspot/D:eps2_Y/D:eps3_Y/D:eps2_Delta/D:eps3_Delta/D\n";
      for (const StatsRow& r : rows) {
        out << r.steinerLength << ' ' << r.perimeter << ' ' << r.junctionOffset << ' ' << r.atQuark;
        for (std::size_t k = 0; k < kShapes.size(); ++k) {
          out << ' ' << r.eps2[k] << ' ' << r.eps3[k];
        }
        out << '\n';
      }
    }
    {
      // After recentring on the quark centroid the spread per coordinate is sqrt(2/3) * quarkSpread.
      const double width2 = (2.0 / 3.0) * params.quarkSpread * params.quarkSpread + params.quarkWidth * params.quarkWidth;
      auto out = openOutput(stem + "_profile.txt");
      out << "r/D:hotspot/D:Y/D:Delta/D:hotspotAnalytic/D\n";
      for (int bin = 0; bin < nBins; ++bin) {
        const double r = (bin + 0.5) * binWidth;
        out << r;
        for (std::size_t k = 0; k < kShapes.size(); ++k) {
          double sum = 0.0;
          for (const auto& profile : blockProfiles) {
            sum += profile[k][bin];
          }
          out << ' ' << sum / (binArea[bin] * nStats);
        }
        out << ' ' << std::exp(-r * r / (2.0 * width2)) / (2.0 * kPi * width2) << '\n';
      }
    }

    int atQuark = 0;
    double meanSteiner = 0.0;
    double meanPerimeter = 0.0;
    std::array<double, 3> meanEps2{};
    for (const StatsRow& r : rows) {
      atQuark += r.atQuark;
      meanSteiner += r.steinerLength / nStats;
      meanPerimeter += r.perimeter / nStats;
      for (std::size_t k = 0; k < kShapes.size(); ++k) {
        meanEps2[k] += r.eps2[k] / nStats;
      }
    }
    const double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
    std::cout << "TRIVENI " << triveni::kVersion << '\n'
              << "triveni-proton: " << nStats << " configurations in " << seconds << " s\n"
              << "  <L_Y> = " << meanSteiner << " fm, <L_Delta> = " << meanPerimeter << " fm\n"
              << "  junction on a quark (angle >= 120 deg): " << 100.0 * atQuark / nStats << " %\n"
              << "  <eps2> single proton: hotspot " << meanEps2[0] << ", Y " << meanEps2[1] << ", Delta "
              << meanEps2[2] << '\n'
              << "  wrote " << stem << "_{maps,quarks,stats,profile}.txt\n";
    return 0;
  } catch (const std::exception& e) {
    std::cerr << "triveni-proton: " << e.what() << '\n';
    return 1;
  }
}
