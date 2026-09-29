// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Banajit Barman

#include "Check.h"

#include "triveni/Grid.h"
#include "triveni/PPCollision.h"
#include "triveni/Proton.h"
#include "triveni/Random.h"

#include <cmath>

int main()
{
  check::near(triveni::reducedThickness(2.0, 8.0, 0.0), 4.0, 1e-12, "p = 0 is the geometric mean");
  check::near(triveni::reducedThickness(2.0, 8.0, 1.0), 5.0, 1e-12, "p = 1 is the arithmetic mean");
  check::near(triveni::reducedThickness(2.0, 8.0, -1.0), 3.2, 1e-12, "p = -1 is the harmonic mean");
  check::near(triveni::reducedThickness(0.0, 8.0, -1.0), 0.0, 1e-12, "p < 0 vanishes if one side is empty");

  // Two equal Gaussians at +-d: eps2 = d^2 / (d^2 + 2 w^2), eps3 = 0.
  // Sources are cut 5 widths from their centre (tail ~ e^-12.5), hence the 1e-5 tolerance.
  {
    const double d = 0.6;
    const double w = 0.3;
    triveni::Proton p;
    p.hotSpots.push_back({{-d, 0.0}, w, 0.5, 1.0});
    p.hotSpots.push_back({{+d, 0.0}, w, 0.5, 1.0});
    triveni::Grid grid(3.0, 0.02);
    grid.addProton(p, {0.0, 0.0}, false);
    const triveni::Eccentricities e = triveni::eccentricities(grid);
    check::near(e.eps[2], d * d / (d * d + 2.0 * w * w), 1e-5, "two-blob eps2");
    check::near(e.eps[3], 0.0, 1e-5, "two-blob eps3");
  }

  // sigma_gg calibration, then an independent check of sigma_inel from the accepted fraction.
  for (auto shape : {triveni::ProtonShape::HotSpot, triveni::ProtonShape::Y}) {
    triveni::PPParams params;
    params.proton.shape = shape;
    params.gridStep = 0.08;
    triveni::PPCollider collider(params);
    const triveni::Calibration cal = collider.calibrate(123, 4000, 4);
    const std::string label = triveni::toString(shape);
    check::near(cal.sigmaInel, params.sigmaInel, 1e-6 * params.sigmaInel, label + ": calibration hits the target");
    check::expect(cal.edgeProbability < 1e-3, label + ": P(b) negligible at bMax");

    triveni::Random rng(456);
    triveni::PPCollider::Workspace ws(params);
    triveni::PPEvent event;
    constexpr int trials = 6000;
    int accepted = 0;
    for (int i = 0; i < trials; ++i) {
      accepted += collider.tryEvent(rng, ws, event) ? 1 : 0;
    }
    const double area = 3.14159265358979323846 * params.bMax * params.bMax;
    const double f = static_cast<double>(accepted) / trials;
    const double sigma = 10.0 * area * f;
    const double error = std::hypot(10.0 * area * std::sqrt(f * (1.0 - f) / trials), cal.sigmaInelError);
    check::near(sigma, params.sigmaInel, 4.0 * error, label + ": accepted fraction reproduces sigma_inel");
  }

  return check::result("test_pp");
}
