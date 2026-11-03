// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Banajit Barman

/// \file Nucleus.h
/// \brief Woods–Saxon nucleus: density normalised to A and its thickness T_A(b).
///
/// rho(r) = rho0 / (1 + exp((r - R)/a)), with rho0 fixed by int rho d^3r = A.
/// The presets use the same R and a as eSTARlight, so benchmarks against it
/// compare physics, not nuclear inputs.

#pragma once

#include <string>
#include <vector>

namespace triveni
{

struct WoodsSaxon {
  int A = 197;
  double radius = 6.38;      ///< fm
  double diffuseness = 0.53; ///< fm

  static WoodsSaxon gold() { return {197, 6.38, 0.53}; }
  static WoodsSaxon lead() { return {208, 6.624, 0.53}; }
  static WoodsSaxon copper() { return {63, 4.214, 0.53}; }
  /// "Au", "Pb" or "Cu".
  static WoodsSaxon fromName(const std::string& name);
};

class Nucleus
{
 public:
  explicit Nucleus(const WoodsSaxon& ws);

  const WoodsSaxon& params() const { return mWs; }
  int massNumber() const { return mWs.A; }

  double rho0() const { return mRho0; }       ///< fm^-3
  double density(double r) const;             ///< fm^-3, r in fm
  double thickness(double b) const;           ///< fm^-2, b in fm; int T d^2b = A
  double maxImpactParameter() const { return mBMax; } ///< fm, beyond which T is negligible

 private:
  WoodsSaxon mWs;
  double mRho0 = 0.0;
  double mBMax = 0.0;
  double mStep = 0.0;
  std::vector<double> mThickness; ///< tabulated on b = i * mStep
};

} // namespace triveni
