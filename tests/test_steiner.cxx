// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Banajit Barman

#include "Check.h"

#include "triveni/Random.h"
#include "triveni/Steiner.h"

#include <algorithm>
#include <array>
#include <cmath>

namespace
{

constexpr double kPi = 3.14159265358979323846;

triveni::Vec3 unit(triveni::Vec3 v) { return (1.0 / triveni::norm(v)) * v; }

double angleAt(const std::array<triveni::Vec3, 3>& q, int i)
{
  const triveni::Vec3 u = q[(i + 1) % 3] - q[i];
  const triveni::Vec3 v = q[(i + 2) % 3] - q[i];
  return std::acos(std::clamp(triveni::dot(u, v) / (triveni::norm(u) * triveni::norm(v)), -1.0, 1.0));
}

/// Weiszfeld iteration: an independent minimiser of the summed distance.
triveni::Vec3 weiszfeld(const std::array<triveni::Vec3, 3>& q)
{
  triveni::Vec3 p = (1.0 / 3.0) * (q[0] + q[1] + q[2]);
  // An unconverged result is only longer, so the "never longer" check below stays valid.
  for (int iter = 0; iter < 5000; ++iter) {
    triveni::Vec3 num;
    double den = 0.0;
    for (const auto& v : q) {
      const double d = triveni::norm(v - p);
      if (d < 1e-14) {
        return v;
      }
      num += (1.0 / d) * v;
      den += 1.0 / d;
    }
    const triveni::Vec3 next = (1.0 / den) * num;
    const double move = triveni::norm(next - p);
    p = next;
    if (move < 1e-15) {
      break;
    }
  }
  return p;
}

} // namespace

int main()
{
  {
    const std::array<triveni::Vec3, 3> q = {triveni::Vec3{1, 0, 0}, triveni::Vec3{-0.5, std::sqrt(3.0) / 2, 0},
                                        triveni::Vec3{-0.5, -std::sqrt(3.0) / 2, 0}};
    const triveni::FermatPoint f = triveni::fermatPoint(q);
    check::expect(f.vertex == -1, "equilateral: junction is interior");
    check::near(triveni::norm(f.point), 0.0, 1e-12, "equilateral: junction at the centre");
    check::near(triveni::legLength(q, f.point), 3.0, 1e-12, "equilateral: total length 3");
  }
  {
    const std::array<triveni::Vec3, 3> q = {triveni::Vec3{0, 0, 0}, triveni::Vec3{1, 0, 0}, triveni::Vec3{-0.5, 0.1, 0}};
    const triveni::FermatPoint f = triveni::fermatPoint(q);
    check::expect(f.vertex == 0, "169 deg angle: junction sits on that quark");
  }
  {
    const std::array<triveni::Vec3, 3> q = {triveni::Vec3{0, 0, 0}, triveni::Vec3{1, 0, 0}, triveni::Vec3{3, 0, 0}};
    check::expect(triveni::fermatPoint(q).vertex == 1, "collinear: junction on the middle quark");
  }
  {
    const std::array<triveni::Vec3, 3> q = {triveni::Vec3{0.2, 0.1, 0}, triveni::Vec3{0.2, 0.1, 0}, triveni::Vec3{1, 1, 1}};
    const triveni::FermatPoint f = triveni::fermatPoint(q);
    check::near(triveni::norm(f.point - q[0]), 0.0, 1e-15, "coincident quarks: junction on the pair");
  }

  // Random 3D triangles: interior junctions have 120 deg legs, vertex junctions an angle >= 120 deg,
  // and the result agrees with an independent minimiser.
  triveni::Random rng(7);
  int interior = 0;
  int worseThanWeiszfeld = 0;
  double maxLegImbalance = 0.0;
  for (int i = 0; i < 20000; ++i) {
    std::array<triveni::Vec3, 3> q{};
    for (auto& v : q) {
      v = {rng.normal(), rng.normal(), rng.normal()};
    }
    const triveni::FermatPoint f = triveni::fermatPoint(q);
    if (f.vertex < 0) {
      ++interior;
      const triveni::Vec3 sum = unit(q[0] - f.point) + unit(q[1] - f.point) + unit(q[2] - f.point);
      maxLegImbalance = std::max(maxLegImbalance, triveni::norm(sum));
    } else {
      check::expect(angleAt(q, f.vertex) >= 2.0 * kPi / 3.0 - 1e-12, "vertex junction has an angle >= 120 deg");
    }
    if (triveni::legLength(q, f.point) > triveni::legLength(q, weiszfeld(q)) + 1e-9) {
      ++worseThanWeiszfeld;
    }
  }
  check::near(maxLegImbalance, 0.0, 1e-9, "interior junction: unit legs sum to zero (120 deg)");
  check::expect(worseThanWeiszfeld == 0, "never longer than the Weiszfeld minimum");
  check::expect(interior > 10000, "most random triangles have an interior junction");

  return check::result("test_steiner");
}
