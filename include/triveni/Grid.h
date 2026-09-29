// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Banajit Barman

/// \file Grid.h
/// \brief Square transverse grid, proton deposition and eccentricities.

#pragma once

#include "triveni/Vec.h"

#include <array>
#include <vector>

namespace triveni
{

struct Proton;

class Grid
{
 public:
  /// Square [-halfWidth, halfWidth]^2; the step is adjusted so the cells tile it exactly.
  Grid(double halfWidth, double step);

  int size() const { return mN; }
  double step() const { return mStep; }
  double halfWidth() const { return mHalf; }
  double cellArea() const { return mStep * mStep; }

  /// Coordinate of the centre of cell i.
  double coord(int i) const { return -mHalf + (i + 0.5) * mStep; }

  double& operator()(int ix, int iy) { return mValues[static_cast<std::size_t>(iy) * mN + ix]; }
  double operator()(int ix, int iy) const { return mValues[static_cast<std::size_t>(iy) * mN + ix]; }

  std::vector<double>& values() { return mValues; }
  const std::vector<double>& values() const { return mValues; }

  void clear();

  /// Sum of the cell values times the cell area.
  double integral() const;

  /// Adds scale * T_p(x - offset). Each source is evaluated only within 5 widths of itself.
  void addProton(const Proton& proton, Vec2 offset, bool fluctuated, double scale = 1.0);

 private:
  int cellIndex(double x) const;

  int mN = 0;
  double mHalf = 0.0;
  double mStep = 0.0;
  std::vector<double> mValues;
};

struct Eccentricities {
  Vec2 centroid;
  /// eps[n] = |sum s r^n e^{i n phi}| / sum s r^n for n = 1..4. eps[0] is unused and
  /// eps[1] vanishes about the centroid by construction.
  std::array<double, 5> eps{};
};

/// Eccentricities of a non-negative density, about its own centroid.
Eccentricities eccentricities(const Grid& density);

} // namespace triveni
