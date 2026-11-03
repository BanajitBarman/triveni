// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Banajit Barman

#include "Check.h"

#include "triveni/Bessel.h"
#include "triveni/Exclusive.h"
#include "triveni/Nucleus.h"
#include "triveni/VectorMeson.h"

#include <cmath>

int main()
{
  constexpr double pi = 3.14159265358979323846;

  // Bessel functions against reference values.
  check::near(triveni::besselJ0(1.0), 0.7651976866, 1e-7, "J0(1)");
  check::near(triveni::besselJ0(10.0), -0.2459357645, 1e-7, "J0(10)");
  check::near(triveni::besselK0(1.0) / 0.4210244382, 1.0, 3e-7, "K0(1)");
  check::near(triveni::besselK1(1.0) / 0.6019072302, 1.0, 3e-7, "K1(1)");
  check::near(triveni::besselK0(0.1) / 2.4270690247, 1.0, 3e-7, "K0(0.1)");
  check::near(triveni::besselK1(0.1) / 9.8538447809, 1.0, 3e-7, "K1(0.1)");
  check::near(triveni::besselK0(5.0) / 0.0036910983, 1.0, 3e-7, "K0(5)");
  check::near(triveni::besselK1(5.0) / 0.0040446134, 1.0, 3e-7, "K1(5)");

  // J/psi boosted Gaussian: the transverse wavefunction is normalised.
  const auto jpsi = triveni::VectorMesonWavefunction::jpsi();
  check::near(triveni::PhotonVectorMesonOverlap(jpsi).normalisationT(), 1.0, 5e-3, "J/psi normalisation (KMW parameters)");

  // Woods–Saxon gold: normalised to A, rho0 equal to eSTARlight's.
  const triveni::Nucleus au(triveni::WoodsSaxon::gold());
  double integral = 0.0;
  const int n = 20000;
  for (int i = 0; i < n; ++i) {
    const double b = (i + 0.5) * au.maxImpactParameter() / n;
    integral += 2.0 * pi * b * au.thickness(b) * au.maxImpactParameter() / n;
  }
  check::near(integral, 197.0, 1e-3, "Au: int T_A d^2b = A");
  check::near(au.rho0(), 0.169551, 1e-6, "Au: rho0 equals eSTARlight's");

  // Proton dipole cross section: dilute limit 2 r^2 Omega, and a numeric b integral elsewhere.
  triveni::DipoleParams dp{0.3, 0.2, 4.0};
  triveni::CoherentVectorMeson proton(jpsi, dp);
  const double x = 1e-3;
  const double omega = dp.omega0 * std::pow(x, -dp.lambda);
  check::near(proton.protonDipoleCrossSection(1e-3, x) / (2e-6 * omega), 1.0, 1e-6, "sigma_p dilute limit");
  for (double r : {0.5, 3.0, 10.0}) {
    double numeric = 0.0;
    const int nb = 40000;
    const double bMax = 40.0;
    for (int i = 0; i < nb; ++i) {
      const double b = (i + 0.5) * bMax / nb;
      const double tp = std::exp(-b * b / (2.0 * dp.protonB)) / (2.0 * pi * dp.protonB);
      numeric += 2.0 * 2.0 * pi * b * (1.0 - std::exp(-r * r * omega * tp)) * bMax / nb;
    }
    check::near(proton.protonDipoleCrossSection(r, x) / numeric, 1.0, 1e-5, "sigma_p closed form at r = " + std::to_string(r));
  }

  // Dilute limit: A_Au(Delta) / A_p(Delta) = F_A(Delta) / F_p(Delta) exactly, so at t = 0
  // dsigma/dt(Au) / dsigma/dt(p) = A^2.
  {
    triveni::DipoleParams dilute{1e-7, 0.2, 4.0};
    triveni::CoherentVectorMeson p(jpsi, dilute);
    triveni::CoherentVectorMeson a(jpsi, dilute, triveni::WoodsSaxon::gold());
    std::vector<double> pT, pL, aT, aL;
    p.dsdtOnGrid(90.0, 2.0, pT, pL);
    a.dsdtOnGrid(90.0, 2.0, aT, aL);
    check::near(aT.front() / pT.front() / (197.0 * 197.0), 1.0, 2e-3, "dilute limit: T, t = 0 ratio is A^2");
    check::near(aL.front() / pL.front() / (197.0 * 197.0), 1.0, 2e-3, "dilute limit: L, t = 0 ratio is A^2");
  }

  // Photon flux: positive inside the kinematic range, zero outside.
  check::expect(triveni::electronPhotonFlux(100.0, 1.0, 4220.0) > 0.0, "flux positive");
  check::expect(triveni::electronPhotonFlux(5000.0, 1.0, 4220.0) == 0.0, "flux zero for E > E_e");
  check::near(triveni::photonPolarisation(1e-9, 4220.0), 1.0, 1e-9, "epsilon -> 1 at small y");

  return check::result("test_exclusive");
}
