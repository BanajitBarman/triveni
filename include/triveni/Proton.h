// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Banajit Barman

/// \file Proton.h
/// \brief Transverse structure of a proton: three valence-quark hot spots joined
///        by gluonic flux tubes in a Y (junction) or Delta topology.
///
/// Geometry is built in 3D and projected on the transverse plane. A flux tube is
/// the segment between its end points convolved with a 3D Gaussian of width
/// `stringWidth`; its projection along the beam is then exactly the projected
/// segment convolved with a 2D Gaussian, which has a closed form (erf).
///
/// The thickness T_p(x) (fm^-2) is normalised to 1:
///   T_p = (1 - f_s)/3 * sum_quarks G(x - q_i)  +  f_s * sum_tubes (L_j / L) * S_j(x),
/// with f_s = stringFraction, L_j the 3D length of tube j and L the total.
/// Tube weights follow the string energy sigma * L_j.

#pragma once

#include "triveni/Random.h"
#include "triveni/Vec.h"

#include <array>
#include <string>
#include <vector>

namespace triveni
{

class Settings;

enum class ProtonShape { HotSpot, Y, Delta, YDelta };

ProtonShape shapeFromString(const std::string& name);
std::string toString(ProtonShape shape);

enum class Topology { None, Y, Delta };

struct ProtonParams {
  ProtonShape shape = ProtonShape::Y;
  double quarkSpread = 0.36;   ///< fm, sigma of each coordinate of the quark positions
  double quarkWidth = 0.17;    ///< fm, width of the hot spot around each quark
  double stringWidth = 0.17;   ///< fm, width of a flux tube
  double stringFraction = 0.5; ///< share of T_p carried by the flux tubes
  double dYDelta = 0.5;        ///< fm, YDelta switch: mean quark separation below which Delta is used
  double fluctuationK = 0.0;   ///< Gamma(k, 1/k) weight fluctuations; <= 0 disables

  static ProtonParams fromSettings(const Settings& settings);
};

/// Gaussian gluon cloud around a valence quark (transverse plane).
struct HotSpot {
  Vec2 centre;
  double width = 0.0;
  double weight = 0.0;      ///< Share of the normalised thickness.
  double fluctuation = 1.0; ///< Gamma-distributed factor, mean 1.
};

/// Projected flux tube between two transverse points.
struct FluxTube {
  Vec2 a;
  Vec2 b;
  double width = 0.0;
  double weight = 0.0;
  double fluctuation = 1.0;
  double length = 0.0; ///< 3D length (fm), before projection.
};

struct Proton {
  std::array<Vec3, 3> quarks{}; ///< Valence-quark positions (fm), centred on the proton's centre of mass.
  Vec3 junction;                ///< Fermat point of the quarks (computed for every shape).
  int junctionOnQuark = -1;     ///< Quark the junction sits on (angle >= 120 deg), or -1.
  Topology topology = Topology::None;
  double stringLength = 0.0;    ///< Total 3D flux-tube length (fm); string energy = sigma * stringLength.
  std::vector<HotSpot> hotSpots;
  std::vector<FluxTube> tubes;

  /// Thickness at transverse point x (fm^-2). Without fluctuations it integrates to 1.
  double thickness(Vec2 x, bool fluctuated = false) const;
};

class ProtonBuilder
{
 public:
  explicit ProtonBuilder(const ProtonParams& params);

  const ProtonParams& params() const { return mParams; }

  /// Three quark positions, each coordinate Gaussian with sigma = quarkSpread.
  std::array<Vec3, 3> sampleQuarks(Random& rng) const;

  /// Builds a proton from given quark positions. `rng` is only used for weight fluctuations,
  /// so the same quarks can be built with several shapes for a like-for-like comparison.
  Proton build(const std::array<Vec3, 3>& quarks, ProtonShape shape, Random& rng) const;
  Proton build(const std::array<Vec3, 3>& quarks, Random& rng) const { return build(quarks, mParams.shape, rng); }

  Proton generate(Random& rng) const { return build(sampleQuarks(rng), rng); }

 private:
  ProtonParams mParams;
};

/// Normalised 2D Gaussian of width w evaluated at offset d.
double gaussian2D(Vec2 d, double width);

/// Segment a-b convolved with a normalised 2D Gaussian of width w; integrates to 1.
double smearedSegment(Vec2 x, Vec2 a, Vec2 b, double width);

} // namespace triveni
