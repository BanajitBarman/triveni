// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Banajit Barman

/// triveni-exclusive: coherent J/psi production, e + target -> e + J/psi + target.
///
/// The proton dipole model is fitted to eSTARlight's gamma p -> J/psi p
/// parametrisation, so that comparisons with eSTARlight test the nuclear and Q^2
/// physics rather than the input. Writes
///   <stem>_fit.txt  sigma(gamma p -> J/psi p) at Q^2 = 0: target vs model
///   <stem>_q2.txt   e + target cross section in each Q^2 bin (W in [wMin, wMax])
///   <stem>_t.txt    flux-weighted dsigma/d|t| in each Q^2 bin
///   <stem>.log      settings, fitted parameters, totals

#include "triveni/Exclusive.h"
#include "triveni/Settings.h"
#include "triveni/Version.h"

#include <chrono>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace
{

std::vector<double> parseList(const std::string& text)
{
  std::vector<double> out;
  std::stringstream ss(text);
  std::string item;
  while (std::getline(ss, item, ',')) {
    out.push_back(std::stod(item));
  }
  return out;
}

std::ofstream openOutput(const std::string& path)
{
  std::ofstream out(path);
  if (!out) {
    throw std::runtime_error("cannot write '" + path + "'");
  }
  out.precision(7);
  return out;
}

/// eSTARlight's sigma(gamma p -> J/psi p) in nb (photonNucleusCrossSection::sigmagp).
double estarlightJpsiProton(double W)
{
  const double threshold = 1.0 - std::pow(triveni::constants::jpsiMass + triveni::constants::protonMass, 2) / (W * W);
  return threshold * threshold * 4.06 * std::pow(W, 0.65);
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
    const std::string stem = settings.word("Main:output");
    const int threads = settings.mode("Main:threads");
    const std::string targetName = settings.word("Exclusive:target");
    const double gammaE = settings.parm("Exclusive:electronGamma");
    const double gammaT = settings.parm("Exclusive:targetGamma");
    const double wMin = settings.parm("Exclusive:wMin");
    const double wMax = settings.parm("Exclusive:wMax");
    const auto q2Edges = parseList(settings.word("Exclusive:q2Edges"));
    const auto fitW = parseList(settings.word("Exclusive:fitW"));
    const double protonB = settings.parm("Exclusive:protonB");

    // Electron energy in the target rest frame, as eSTARlight computes it.
    const double electronEnergy = triveni::constants::electronMass * std::cosh(std::acosh(gammaE) + std::acosh(gammaT));

    const auto jpsi = triveni::VectorMesonWavefunction::jpsi();
    const triveni::DipoleParams dipole = triveni::fitDipoleToPhotoproduction(jpsi, protonB, estarlightJpsiProton, fitW);

    {
      triveni::CoherentVectorMeson proton(jpsi, dipole);
      auto out = openOutput(stem + "_fit.txt");
      out << "W/D:estarlight/D:triveni/D\n";
      for (double W = 20.0; W <= 250.0; W *= 1.1) {
        out << W << ' ' << estarlightJpsiProton(W) << ' ' << proton.sigma(W, 0.0).T << '\n';
      }
    }

    std::optional<triveni::WoodsSaxon> nucleus;
    if (targetName != "p") {
      nucleus = triveni::WoodsSaxon::fromName(targetName);
    }
    const triveni::CoherentVectorMeson model(jpsi, dipole, nucleus);

    auto q2Out = openOutput(stem + "_q2.txt");
    auto tOut = openOutput(stem + "_t.txt");
    q2Out << "q2lo/D:q2hi/D:sigma/D\n";
    tOut << "bin/I:absT/D:dsdt/D\n";
    double total = 0.0;
    for (std::size_t i = 0; i + 1 < q2Edges.size(); ++i) {
      const auto res = triveni::electronTargetCrossSection(model, electronEnergy, wMin, wMax, q2Edges[i], q2Edges[i + 1],
                                                           threads);
      total += res.sigma;
      q2Out << q2Edges[i] << ' ' << q2Edges[i + 1] << ' ' << res.sigma << '\n';
      for (std::size_t k = 0; k < res.absT.size(); ++k) {
        tOut << i << ' ' << res.absT[k] << ' ' << res.dsdt[k] << '\n';
      }
      std::cout << "  Q2 in [" << q2Edges[i] << ", " << q2Edges[i + 1] << "] GeV^2: sigma = " << res.sigma << " nb\n";
    }

    const double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
    std::ostringstream summary;
    summary << "TRIVENI " << triveni::kVersion << '\n'
            << "triveni-exclusive: coherent J/psi on " << targetName << ", E_e(target frame) = " << electronEnergy
            << " GeV, W in [" << wMin << ", " << std::min(wMax, std::sqrt(2.0 * 0.938272046 * electronEnergy)) << "] GeV\n"
            << "  fitted dipole: Omega_0 = " << dipole.omega0 << ", lambda = " << dipole.lambda << ", B_p = " << dipole.protonB
            << " GeV^-2\n"
            << "  total over Q2 in [" << q2Edges.front() << ", " << q2Edges.back() << "]: " << total << " nb (" << seconds
            << " s)\n";
    std::cout << summary.str();
    auto log = openOutput(stem + ".log");
    log << summary.str() << "\nSettings (* = changed):\n";
    settings.list(log);
    return 0;
  } catch (const std::exception& e) {
    std::cerr << "triveni-exclusive: " << e.what() << '\n';
    return 1;
  }
}
