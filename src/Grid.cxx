// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Banajit Barman

#include "triveni/Grid.h"

#include "triveni/Proton.h"

#include <algorithm>
#include <cmath>
#include <complex>
#include <stdexcept>

namespace triveni
{

namespace
{

constexpr double kReach = 5.0; // sources are evaluated within kReach widths

} // namespace

Grid::Grid(double halfWidth, double step)
{
  if (!(halfWidth > 0.0) || !(step > 0.0) || step > halfWidth) {
    throw std::invalid_argument("Grid: need 0 < step <= halfWidth");
  }
  mN = static_cast<int>(std::lround(2.0 * halfWidth / step));
  mHalf = halfWidth;
  mStep = 2.0 * halfWidth / mN;
  mValues.assign(static_cast<std::size_t>(mN) * mN, 0.0);
}

void Grid::clear() { std::fill(mValues.begin(), mValues.end(), 0.0); }

double Grid::integral() const
{
  double sum = 0.0;
  for (double v : mValues) {
    sum += v;
  }
  return sum * cellArea();
}

int Grid::cellIndex(double x) const
{
  return std::clamp(static_cast<int>(std::floor((x + mHalf) / mStep)), 0, mN - 1);
}

void Grid::addProton(const Proton& p, Vec2 offset, bool fluctuated, double scale)
{
  for (const HotSpot& h : p.hotSpots) {
    const double w = scale * h.weight * (fluctuated ? h.fluctuation : 1.0);
    const Vec2 c = h.centre + offset;
    const double reach = kReach * h.width;
    const int x0 = cellIndex(c.x - reach), x1 = cellIndex(c.x + reach);
    const int y0 = cellIndex(c.y - reach), y1 = cellIndex(c.y + reach);
    for (int iy = y0; iy <= y1; ++iy) {
      for (int ix = x0; ix <= x1; ++ix) {
        (*this)(ix, iy) += w * gaussian2D(Vec2{coord(ix), coord(iy)} - c, h.width);
      }
    }
  }
  for (const FluxTube& f : p.tubes) {
    const double w = scale * f.weight * (fluctuated ? f.fluctuation : 1.0);
    const Vec2 a = f.a + offset;
    const Vec2 b = f.b + offset;
    const double reach = kReach * f.width;
    const int x0 = cellIndex(std::min(a.x, b.x) - reach), x1 = cellIndex(std::max(a.x, b.x) + reach);
    const int y0 = cellIndex(std::min(a.y, b.y) - reach), y1 = cellIndex(std::max(a.y, b.y) + reach);
    for (int iy = y0; iy <= y1; ++iy) {
      for (int ix = x0; ix <= x1; ++ix) {
        (*this)(ix, iy) += w * smearedSegment(Vec2{coord(ix), coord(iy)}, a, b, f.width);
      }
    }
  }
}

Eccentricities eccentricities(const Grid& g)
{
  Eccentricities result;
  const int n = g.size();
  double total = 0.0;
  double sx = 0.0;
  double sy = 0.0;
  for (int iy = 0; iy < n; ++iy) {
    for (int ix = 0; ix < n; ++ix) {
      const double v = g(ix, iy);
      total += v;
      sx += v * g.coord(ix);
      sy += v * g.coord(iy);
    }
  }
  if (!(total > 0.0)) {
    return result;
  }
  result.centroid = {sx / total, sy / total};

  std::array<std::complex<double>, 5> numerator{};
  std::array<double, 5> denominator{};
  for (int iy = 0; iy < n; ++iy) {
    const double dy = g.coord(iy) - result.centroid.y;
    for (int ix = 0; ix < n; ++ix) {
      const double v = g(ix, iy);
      if (v <= 0.0) {
        continue;
      }
      const std::complex<double> z(g.coord(ix) - result.centroid.x, dy);
      const double r = std::abs(z);
      std::complex<double> zn = 1.0;
      double rn = 1.0;
      for (int k = 1; k <= 4; ++k) {
        zn *= z;
        rn *= r;
        numerator[k] += v * zn;
        denominator[k] += v * rn;
      }
    }
  }
  for (int k = 1; k <= 4; ++k) {
    result.eps[k] = denominator[k] > 0.0 ? std::abs(numerator[k]) / denominator[k] : 0.0;
  }
  return result;
}

} // namespace triveni
