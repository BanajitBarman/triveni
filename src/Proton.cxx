// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Banajit Barman

#include "triveni/Proton.h"

#include "triveni/Settings.h"
#include "triveni/Steiner.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <stdexcept>
#include <utility>

namespace triveni
{

namespace
{

constexpr double kPi = 3.14159265358979323846;

} // namespace

ProtonShape shapeFromString(const std::string& name)
{
  std::string key = name;
  std::transform(key.begin(), key.end(), key.begin(), [](unsigned char c) { return std::tolower(c); });
  if (key == "hotspot") {
    return ProtonShape::HotSpot;
  }
  if (key == "y") {
    return ProtonShape::Y;
  }
  if (key == "delta") {
    return ProtonShape::Delta;
  }
  if (key == "ydelta") {
    return ProtonShape::YDelta;
  }
  throw std::invalid_argument("Unknown proton shape '" + name + "' (hotspot | Y | Delta | YDelta)");
}

std::string toString(ProtonShape shape)
{
  switch (shape) {
    case ProtonShape::HotSpot:
      return "hotspot";
    case ProtonShape::Y:
      return "Y";
    case ProtonShape::Delta:
      return "Delta";
    case ProtonShape::YDelta:
      return "YDelta";
  }
  return "?";
}

ProtonParams ProtonParams::fromSettings(const Settings& s)
{
  ProtonParams p;
  p.shape = shapeFromString(s.word("Proton:shape"));
  p.quarkSpread = s.parm("Proton:quarkSpread");
  p.quarkWidth = s.parm("Proton:quarkWidth");
  p.stringWidth = s.parm("Proton:stringWidth");
  p.stringFraction = s.parm("Proton:stringFraction");
  p.dYDelta = s.parm("Proton:dYDelta");
  p.fluctuationK = s.parm("Proton:fluctuationK");
  if (p.quarkSpread < 0.0 || p.quarkWidth <= 0.0 || p.stringWidth <= 0.0) {
    throw std::invalid_argument("Proton: widths must be positive");
  }
  if (p.stringFraction < 0.0 || p.stringFraction > 1.0) {
    throw std::invalid_argument("Proton:stringFraction must be in [0, 1]");
  }
  return p;
}

double gaussian2D(Vec2 d, double width)
{
  const double w2 = width * width;
  return std::exp(-dot(d, d) / (2.0 * w2)) / (2.0 * kPi * w2);
}

double smearedSegment(Vec2 x, Vec2 a, Vec2 b, double width)
{
  const Vec2 ab = b - a;
  const double length = norm(ab);
  if (length < 1e-6 * width) {
    return gaussian2D(x - (0.5 * (a + b)), width);
  }
  const Vec2 u = (1.0 / length) * ab;
  const Vec2 d = x - a;
  const double t = dot(d, u);                           // along the segment
  const double perp2 = std::max(0.0, dot(d, d) - t * t); // squared distance from its line
  const double s = std::sqrt(2.0) * width;
  // erf(t/s) - erf((t-L)/s), written with erfc to avoid cancellation far from the segment.
  const double along = (t > 0.5 * length) ? std::erfc((t - length) / s) - std::erfc(t / s)
                                          : std::erfc(-t / s) - std::erfc((length - t) / s);
  return std::exp(-perp2 / (2.0 * width * width)) / (std::sqrt(2.0 * kPi) * width) * along / (2.0 * length);
}

double Proton::thickness(Vec2 x, bool fluctuated) const
{
  double t = 0.0;
  for (const HotSpot& h : hotSpots) {
    t += h.weight * (fluctuated ? h.fluctuation : 1.0) * gaussian2D(x - h.centre, h.width);
  }
  for (const FluxTube& f : tubes) {
    t += f.weight * (fluctuated ? f.fluctuation : 1.0) * smearedSegment(x, f.a, f.b, f.width);
  }
  return t;
}

ProtonBuilder::ProtonBuilder(const ProtonParams& params) : mParams(params) {}

std::array<Vec3, 3> ProtonBuilder::sampleQuarks(Random& rng) const
{
  std::array<Vec3, 3> q{};
  for (Vec3& v : q) {
    v = {mParams.quarkSpread * rng.normal(), mParams.quarkSpread * rng.normal(), mParams.quarkSpread * rng.normal()};
  }
  return q;
}

Proton ProtonBuilder::build(const std::array<Vec3, 3>& quarks, ProtonShape shape, Random& rng) const
{
  Proton p;
  p.quarks = quarks;
  const FermatPoint junction = fermatPoint(quarks);
  p.junction = junction.point;
  p.junctionOnQuark = junction.vertex;

  if (shape == ProtonShape::YDelta) {
    const double meanSeparation =
      (norm(quarks[0] - quarks[1]) + norm(quarks[1] - quarks[2]) + norm(quarks[2] - quarks[0])) / 3.0;
    shape = (meanSeparation < mParams.dYDelta) ? ProtonShape::Delta : ProtonShape::Y;
  }

  // Flux tubes in 3D.
  std::vector<std::pair<Vec3, Vec3>> segments;
  if (shape == ProtonShape::Y) {
    p.topology = Topology::Y;
    for (int i = 0; i < 3; ++i) {
      if (i != junction.vertex) {
        segments.emplace_back(quarks[i], junction.point);
      }
    }
  } else if (shape == ProtonShape::Delta) {
    p.topology = Topology::Delta;
    for (int i = 0; i < 3; ++i) {
      segments.emplace_back(quarks[i], quarks[(i + 1) % 3]);
    }
  }
  double totalLength = 0.0;
  for (const auto& [a, b] : segments) {
    totalLength += norm(b - a);
  }
  p.stringLength = totalLength;
  const double fs = (totalLength > 1e-9) ? mParams.stringFraction : 0.0;

  // Recentre on the centre of mass of the unfluctuated density.
  Vec3 centre = ((1.0 - fs) / 3.0) * (quarks[0] + quarks[1] + quarks[2]);
  if (fs > 0.0) {
    for (const auto& [a, b] : segments) {
      centre += (fs * norm(b - a) / totalLength) * (0.5 * (a + b));
    }
  }
  for (Vec3& q : p.quarks) {
    q = q - centre;
  }
  p.junction = p.junction - centre;
  for (auto& [a, b] : segments) {
    a = a - centre;
    b = b - centre;
  }

  const bool fluctuate = mParams.fluctuationK > 0.0;
  auto draw = [&]() { return fluctuate ? rng.gamma(mParams.fluctuationK) : 1.0; };

  for (const Vec3& q : p.quarks) {
    p.hotSpots.push_back({transverse(q), mParams.quarkWidth, (1.0 - fs) / 3.0, draw()});
  }
  if (fs > 0.0) {
    for (const auto& [a, b] : segments) {
      const double length = norm(b - a);
      if (length <= 0.0) {
        continue;
      }
      p.tubes.push_back({transverse(a), transverse(b), mParams.stringWidth, fs * length / totalLength, draw(), length});
    }
  }
  return p;
}

} // namespace triveni
