// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Banajit Barman

#include "triveni/Random.h"

#include <array>
#include <cmath>
#include <stdexcept>

namespace triveni
{

namespace
{

std::uint64_t splitmix64(std::uint64_t& state)
{
  std::uint64_t z = (state += 0x9E3779B97F4A7C15ULL);
  z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
  z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
  return z ^ (z >> 31);
}

} // namespace

Random::Random(std::uint64_t seed) { this->seed(seed); }

void Random::seed(std::uint64_t seed)
{
  // Spread the 64-bit seed over the full Mersenne-Twister state.
  std::uint64_t state = seed;
  std::array<std::uint32_t, 8> words{};
  for (auto& word : words) {
    word = static_cast<std::uint32_t>(splitmix64(state) >> 32);
  }
  std::seed_seq sequence(words.begin(), words.end());
  mEngine.seed(sequence);
  mHasSpare = false;
}

double Random::uniform()
{
  // 53 random bits, shifted by half a step so that 0 and 1 are never returned.
  return (static_cast<double>(mEngine() >> 11) + 0.5) * 0x1.0p-53;
}

double Random::normal()
{
  if (mHasSpare) {
    mHasSpare = false;
    return mSpare;
  }
  double u = 0.0;
  double v = 0.0;
  double s = 0.0;
  do {
    u = 2.0 * uniform() - 1.0;
    v = 2.0 * uniform() - 1.0;
    s = u * u + v * v;
  } while (s >= 1.0 || s == 0.0);
  const double factor = std::sqrt(-2.0 * std::log(s) / s);
  mSpare = v * factor;
  mHasSpare = true;
  return u * factor;
}

double Random::gamma(double k)
{
  if (!(k > 0.0)) {
    throw std::invalid_argument("Random::gamma: shape must be positive");
  }
  if (k < 1.0) {
    // Boost: Gamma(k) = Gamma(k + 1) * U^(1/k).
    return gammaUnitScale(k + 1.0) * std::pow(uniform(), 1.0 / k) / k;
  }
  return gammaUnitScale(k) / k;
}

double Random::gammaUnitScale(double shape)
{
  // Marsaglia & Tsang, ACM TOMS 26 (2000) 363, valid for shape >= 1.
  const double d = shape - 1.0 / 3.0;
  const double c = 1.0 / std::sqrt(9.0 * d);
  while (true) {
    const double x = normal();
    double v = 1.0 + c * x;
    if (v <= 0.0) {
      continue;
    }
    v = v * v * v;
    const double u = uniform();
    const double x2 = x * x;
    if (u < 1.0 - 0.0331 * x2 * x2) {
      return d * v;
    }
    if (std::log(u) < 0.5 * x2 + d * (1.0 - v + std::log(v))) {
      return d * v;
    }
  }
}

std::uint64_t Random::streamSeed(std::uint64_t seed, std::uint64_t stream)
{
  std::uint64_t state = seed ^ (stream * 0xD1B54A32D192ED03ULL);
  splitmix64(state);
  return splitmix64(state);
}

} // namespace triveni
