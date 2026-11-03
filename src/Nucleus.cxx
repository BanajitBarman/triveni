// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Banajit Barman

#include "triveni/Nucleus.h"

#include <cmath>
#include <stdexcept>

namespace triveni
{

namespace
{

constexpr double kPi = 3.14159265358979323846;

/// Simpson's rule on n (even) intervals.
template <class F>
double simpson(F f, double a, double b, int n)
{
  const double h = (b - a) / n;
  double sum = f(a) + f(b);
  for (int i = 1; i < n; ++i) {
    sum += (i % 2 ? 4.0 : 2.0) * f(a + i * h);
  }
  return sum * h / 3.0;
}

} // namespace

WoodsSaxon WoodsSaxon::fromName(const std::string& name)
{
  if (name == "Au") {
    return gold();
  }
  if (name == "Pb") {
    return lead();
  }
  if (name == "Cu") {
    return copper();
  }
  throw std::invalid_argument("Unknown nucleus '" + name + "' (Au | Pb | Cu)");
}

Nucleus::Nucleus(const WoodsSaxon& ws) : mWs(ws)
{
  if (ws.A < 2 || !(ws.radius > 0.0) || !(ws.diffuseness > 0.0)) {
    throw std::invalid_argument("Nucleus: need A >= 2 and positive R, a");
  }
  // Beyond R + 25a the density is below e^-25 of its central value.
  mBMax = ws.radius + 25.0 * ws.diffuseness;
  const auto shape = [&](double r) { return 1.0 / (1.0 + std::exp((r - ws.radius) / ws.diffuseness)); };
  const double volume = simpson([&](double r) { return 4.0 * kPi * r * r * shape(r); }, 0.0, mBMax, 4000);
  mRho0 = ws.A / volume;

  mStep = 0.005;
  const int n = static_cast<int>(mBMax / mStep) + 2;
  mThickness.resize(n);
  for (int i = 0; i < n; ++i) {
    const double b = i * mStep;
    const double zMax = std::sqrt(std::max(0.0, mBMax * mBMax - b * b));
    mThickness[i] = zMax > 0.0 ? 2.0 * simpson([&](double z) { return density(std::sqrt(b * b + z * z)); }, 0.0, zMax, 800)
                               : 0.0;
  }
}

double Nucleus::density(double r) const
{
  return mRho0 / (1.0 + std::exp((r - mWs.radius) / mWs.diffuseness));
}

double Nucleus::thickness(double b) const
{
  const double u = std::fabs(b) / mStep;
  const auto i = static_cast<std::size_t>(u);
  if (i + 1 >= mThickness.size()) {
    return 0.0;
  }
  const double f = u - i;
  return (1.0 - f) * mThickness[i] + f * mThickness[i + 1];
}

} // namespace triveni
