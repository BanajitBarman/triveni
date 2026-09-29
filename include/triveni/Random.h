// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Banajit Barman

/// \file Random.h
/// \brief Reproducible random numbers.
///
/// The engine is std::mt19937_64 (fully specified by the standard), but the
/// distributions are implemented here: std::normal_distribution and friends use
/// different algorithms in libc++ and libstdc++, so the same seed would give
/// different events on macOS and Linux.

#pragma once

#include <cstdint>
#include <random>

namespace triveni
{

class Random
{
 public:
  explicit Random(std::uint64_t seed = 1);

  void seed(std::uint64_t seed);

  /// Uniform in the open interval (0, 1).
  double uniform();

  /// Standard normal (Marsaglia polar method).
  double normal();

  /// Gamma(k, 1/k): mean 1, variance 1/k (Marsaglia–Tsang). Requires k > 0.
  double gamma(double k);

  /// Seed of an independent stream, e.g. one per block of events, so that
  /// results do not depend on the number of threads.
  static std::uint64_t streamSeed(std::uint64_t seed, std::uint64_t stream);

 private:
  double gammaUnitScale(double shape);

  std::mt19937_64 mEngine;
  bool mHasSpare = false;
  double mSpare = 0.0;
};

} // namespace triveni
