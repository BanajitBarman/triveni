// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Banajit Barman

#include "triveni/Steiner.h"

#include <algorithm>
#include <cmath>

namespace triveni
{

FermatPoint fermatPoint(const std::array<Vec3, 3>& q)
{
  // side[i] is the side opposite quark i.
  const double side[3] = {norm(q[1] - q[2]), norm(q[2] - q[0]), norm(q[0] - q[1])};
  const double scale = std::max({side[0], side[1], side[2]});
  if (scale == 0.0) {
    return {q[0], 0};
  }

  // Two coincident quarks: the doubled point wins (2|P - A| + |P - C| is minimal at A).
  const double tiny = 1e-12 * scale;
  if (side[2] <= tiny || side[1] <= tiny) {
    return {q[0], 0};
  }
  if (side[0] <= tiny) {
    return {q[1], 1};
  }

  double angle[3];
  for (int i = 0; i < 3; ++i) {
    const Vec3 u = q[(i + 1) % 3] - q[i];
    const Vec3 v = q[(i + 2) % 3] - q[i];
    const double c = dot(u, v) / (norm(u) * norm(v));
    angle[i] = std::acos(std::clamp(c, -1.0, 1.0));
  }

  constexpr double pi = 3.14159265358979323846;
  for (int i = 0; i < 3; ++i) {
    if (angle[i] >= 2.0 * pi / 3.0) {
      return {q[i], i};
    }
  }

  // First isogonic centre: barycentric coordinates a_i / sin(A_i + 60 deg).
  double weight[3];
  for (int i = 0; i < 3; ++i) {
    weight[i] = side[i] / std::sin(angle[i] + pi / 3.0);
  }
  const double sum = weight[0] + weight[1] + weight[2];
  const Vec3 point = (weight[0] / sum) * q[0] + (weight[1] / sum) * q[1] + (weight[2] / sum) * q[2];
  return {point, -1};
}

double legLength(const std::array<Vec3, 3>& q, Vec3 point)
{
  return norm(q[0] - point) + norm(q[1] - point) + norm(q[2] - point);
}

} // namespace triveni
