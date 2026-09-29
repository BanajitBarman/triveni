// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Banajit Barman

/// \file PPCollision.h
/// \brief Proton–proton collisions: collision probability, entropy deposition and
///        eccentricities (TRENTo-style), for any proton shape.
///
/// A pair at impact parameter b collides with probability
///   P(b) = 1 - exp(-sigma_gg T_pp(b)),   T_pp(b) = int d^2x T_1(x + b/2) T_2(x - b/2),
/// computed event by event from the sampled configurations. sigma_gg is fixed so
/// that int d^2b <P(b)> equals the measured inelastic cross section.
/// The entropy density is s = ((T1^p + T2^p)/2)^(1/p) with (optionally) fluctuated thicknesses.

#pragma once

#include "triveni/Grid.h"
#include "triveni/Proton.h"
#include "triveni/Random.h"

#include <cstdint>

namespace triveni
{

class Settings;

struct PPParams {
  ProtonParams proton;
  double sigmaInel = 42.0; ///< mb (sqrt(s) = 200 GeV)
  double bMax = 3.5;       ///< fm
  double reducedP = 0.0;
  double gridHalfWidth = 3.5;
  double gridStep = 0.05;

  static PPParams fromSettings(const Settings& settings);
};

struct PPEvent {
  double b = 0.0;          ///< Impact parameter (fm).
  double overlap = 0.0;    ///< T_pp(b) for this pair (fm^-2).
  double pCollision = 0.0; ///< 1 - exp(-sigma_gg T_pp).
  double entropy = 0.0;    ///< int s d^2x (arbitrary units, proportional to multiplicity).
  Eccentricities ecc;
  Proton p1;
  Proton p2;
};

struct Calibration {
  double sigmaGG = 0.0;         ///< fm^2
  double sigmaInel = 0.0;       ///< mb, reproduced by the calibration sample
  double sigmaInelError = 0.0;  ///< mb, statistical
  double edgeProbability = 0.0; ///< Mean P(b) for b > 0.9 bMax; must be ~0 or bMax is too small.
  /// d ln(sigma_inel) / d ln(sigma_gg). Small values (black-disk regime) mean sigma_inel hardly
  /// depends on sigma_gg, so sigma_gg is poorly determined: the proton is too small for sigmaInel.
  double sensitivity = 0.0;
  int pairs = 0;

  static constexpr double kBlackDiskSensitivity = 0.1;
};

/// s(T1, T2) of TRENTo; p = 0 is the geometric mean.
double reducedThickness(double t1, double t2, double p);

class PPCollider
{
 public:
  explicit PPCollider(const PPParams& params);

  const PPParams& params() const { return mParams; }

  /// Scratch grids; use one per thread.
  struct Workspace {
    explicit Workspace(const PPParams& params);
    Grid t1;
    Grid t2;
  };

  /// Fixes sigma_gg from nPairs sampled pairs and stores it.
  Calibration calibrate(std::uint64_t seed, int nPairs, int nThreads);

  void setSigmaGG(double sigmaGG) { mSigmaGG = sigmaGG; }
  double sigmaGG() const { return mSigmaGG; }

  /// One pp encounter at a random impact parameter (uniform in the disk b < bMax).
  /// Returns true if it is inelastic; `event` is then filled and ws.t1 holds s(x).
  bool tryEvent(Random& rng, Workspace& ws, PPEvent& event) const;

 private:
  double sampleOverlap(Random& rng, Workspace& ws, double b, Proton& p1, Proton& p2) const;

  PPParams mParams;
  ProtonBuilder mBuilder;
  double mSigmaGG = -1.0;
};

} // namespace triveni
