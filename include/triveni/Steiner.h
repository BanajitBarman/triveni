// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Banajit Barman

/// \file Steiner.h
/// \brief Position of the baryon junction for three valence quarks.
///
/// The Y-shaped flux tubes minimise their total length, so the junction sits at
/// the Fermat–Torricelli (Steiner) point of the quark triangle: the three legs
/// meet at 120 degrees. If one angle of the triangle is >= 120 degrees, the
/// minimum is on that quark and the Y degenerates into two legs.

#pragma once

#include "triveni/Vec.h"

#include <array>

namespace triveni
{

struct FermatPoint {
  Vec3 point;
  int vertex = -1; ///< Index of the quark the junction sits on, or -1 if it is interior.
};

/// Point minimising |P - q0| + |P - q1| + |P - q2|. Works for any triangle in 3D,
/// including collinear and coincident quarks.
FermatPoint fermatPoint(const std::array<Vec3, 3>& quarks);

/// Total length of the three legs from the quarks to `point`.
double legLength(const std::array<Vec3, 3>& quarks, Vec3 point);

} // namespace triveni
