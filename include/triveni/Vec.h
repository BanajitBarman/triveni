// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Banajit Barman

/// \file Vec.h
/// \brief Minimal 2D (transverse) and 3D vectors. Lengths are in fm.

#pragma once

#include <cmath>

namespace triveni
{

struct Vec2 {
  double x = 0.0;
  double y = 0.0;
};

inline Vec2 operator+(Vec2 a, Vec2 b) { return {a.x + b.x, a.y + b.y}; }
inline Vec2 operator-(Vec2 a, Vec2 b) { return {a.x - b.x, a.y - b.y}; }
inline Vec2 operator*(double s, Vec2 a) { return {s * a.x, s * a.y}; }
inline double dot(Vec2 a, Vec2 b) { return a.x * b.x + a.y * b.y; }
inline double norm(Vec2 a) { return std::hypot(a.x, a.y); }

struct Vec3 {
  double x = 0.0;
  double y = 0.0;
  double z = 0.0;
};

inline Vec3 operator+(Vec3 a, Vec3 b) { return {a.x + b.x, a.y + b.y, a.z + b.z}; }
inline Vec3 operator-(Vec3 a, Vec3 b) { return {a.x - b.x, a.y - b.y, a.z - b.z}; }
inline Vec3 operator*(double s, Vec3 a) { return {s * a.x, s * a.y, s * a.z}; }
inline Vec3& operator+=(Vec3& a, Vec3 b) { return a = a + b; }
inline double dot(Vec3 a, Vec3 b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
inline double norm(Vec3 a) { return std::sqrt(dot(a, a)); }

/// Projection on the transverse (x, y) plane; the beam is along z.
inline Vec2 transverse(Vec3 a) { return {a.x, a.y}; }

} // namespace triveni
