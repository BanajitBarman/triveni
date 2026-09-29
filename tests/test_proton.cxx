// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Banajit Barman

#include "Check.h"

#include "triveni/Grid.h"
#include "triveni/Proton.h"
#include "triveni/Random.h"
#include "triveni/Steiner.h"

#include <array>
#include <cmath>

namespace
{

/// Direct quadrature of the segment convolution, as an independent reference.
double segmentByQuadrature(triveni::Vec2 x, triveni::Vec2 a, triveni::Vec2 b, double w)
{
  constexpr int n = 4000;
  double sum = 0.0;
  for (int i = 0; i < n; ++i) {
    const double s = (i + 0.5) / n;
    sum += triveni::gaussian2D(x - (a + s * (b - a)), w);
  }
  return sum / n;
}

} // namespace

int main()
{
  // Flux-tube profile: normalisation, agreement with quadrature, zero-length limit.
  {
    const triveni::Vec2 a{-0.3, 0.1};
    const triveni::Vec2 b{0.5, 0.4};
    const double w = 0.17;
    triveni::Grid grid(2.0, 0.01);
    double sum = 0.0;
    double worst = 0.0;
    for (int iy = 0; iy < grid.size(); ++iy) {
      for (int ix = 0; ix < grid.size(); ++ix) {
        const triveni::Vec2 x{grid.coord(ix), grid.coord(iy)};
        const double v = triveni::smearedSegment(x, a, b, w);
        sum += v * grid.cellArea();
        if ((ix % 37 == 0) && (iy % 41 == 0)) {
          const double ref = segmentByQuadrature(x, a, b, w);
          worst = std::max(worst, std::abs(v - ref) / (ref + 1e-3));
        }
      }
    }
    check::near(sum, 1.0, 1e-6, "flux tube integrates to 1");
    check::near(worst, 0.0, 1e-6, "flux tube matches direct quadrature");
    const triveni::Vec2 x{0.05, -0.02};
    check::near(triveni::smearedSegment(x, a, a + triveni::Vec2{1e-9, 0}, w), triveni::gaussian2D(x - a, w), 1e-6,
                "zero-length tube is a Gaussian");
    check::near(triveni::smearedSegment(x, a, a + triveni::Vec2{1e-4, 0}, w), triveni::gaussian2D(x - a, w), 1e-2,
                "short tube is close to a Gaussian");
  }

  triveni::ProtonParams params;
  params.fluctuationK = 0.0;

  // Every shape is normalised and centred.
  for (auto shape : {triveni::ProtonShape::HotSpot, triveni::ProtonShape::Y, triveni::ProtonShape::Delta, triveni::ProtonShape::YDelta}) {
    params.shape = shape;
    const triveni::ProtonBuilder builder(params);
    triveni::Random rng(11);
    triveni::Grid grid(3.0, 0.03);
    double worstNorm = 0.0;
    double worstCentre = 0.0;
    for (int i = 0; i < 200; ++i) {
      const triveni::Proton p = builder.generate(rng);
      grid.clear();
      grid.addProton(p, {0.0, 0.0}, false);
      worstNorm = std::max(worstNorm, std::abs(grid.integral() - 1.0));
      worstCentre = std::max(worstCentre, triveni::norm(triveni::eccentricities(grid).centroid));
    }
    check::near(worstNorm, 0.0, 1e-5, triveni::toString(shape) + ": T_p integrates to 1");
    check::near(worstCentre, 0.0, 1e-6, triveni::toString(shape) + ": T_p is centred");
  }

  // Hot spots: <x^2> of T_p averaged over protons is (2/3) quarkSpread^2 + quarkWidth^2.
  {
    params.shape = triveni::ProtonShape::HotSpot;
    const triveni::ProtonBuilder builder(params);
    triveni::Random rng(12);
    constexpr int n = 200000;
    double sum = 0.0;
    for (int i = 0; i < n; ++i) {
      const triveni::Proton p = builder.generate(rng);
      for (const auto& h : p.hotSpots) {
        sum += h.weight * (h.centre.x * h.centre.x + h.width * h.width);
      }
    }
    const double expected = (2.0 / 3.0) * params.quarkSpread * params.quarkSpread + params.quarkWidth * params.quarkWidth;
    check::near(sum / n, expected, 0.01 * expected, "hot-spot mean square radius");
  }

  // Topology details.
  {
    triveni::Random rng(13);
    const std::array<triveni::Vec3, 3> q = {triveni::Vec3{0.4, 0, 0.1}, triveni::Vec3{-0.2, 0.35, -0.1}, triveni::Vec3{-0.2, -0.3, 0}};
    const double perimeter = triveni::norm(q[0] - q[1]) + triveni::norm(q[1] - q[2]) + triveni::norm(q[2] - q[0]);

    params.shape = triveni::ProtonShape::Y;
    triveni::ProtonBuilder builder(params);
    const triveni::Proton y = builder.build(q, rng);
    check::expect(y.topology == triveni::Topology::Y && y.tubes.size() == 3, "Y: three legs");
    check::near(y.stringLength, triveni::legLength(q, triveni::fermatPoint(q).point), 1e-12, "Y: string length");

    const triveni::Proton d = builder.build(q, triveni::ProtonShape::Delta, rng);
    check::expect(d.topology == triveni::Topology::Delta && d.tubes.size() == 3, "Delta: three sides");
    check::near(d.stringLength, perimeter, 1e-12, "Delta: string length is the perimeter");

    params.dYDelta = 10.0;
    check::expect(triveni::ProtonBuilder(params).build(q, triveni::ProtonShape::YDelta, rng).topology == triveni::Topology::Delta,
                  "YDelta: small proton is Delta");
    params.dYDelta = 0.0;
    check::expect(triveni::ProtonBuilder(params).build(q, triveni::ProtonShape::YDelta, rng).topology == triveni::Topology::Y,
                  "YDelta: large proton is Y");

    params.stringFraction = 0.0;
    check::expect(triveni::ProtonBuilder(params).build(q, triveni::ProtonShape::Y, rng).tubes.empty(),
                  "stringFraction 0: no tubes");
    params = triveni::ProtonParams{};
  }

  // Weight fluctuations have mean 1.
  {
    params.shape = triveni::ProtonShape::Y;
    params.fluctuationK = 2.0;
    const triveni::ProtonBuilder builder(params);
    triveni::Random rng(14);
    double sum = 0.0;
    int count = 0;
    for (int i = 0; i < 50000; ++i) {
      const triveni::Proton p = builder.generate(rng);
      for (const auto& h : p.hotSpots) {
        sum += h.fluctuation;
        ++count;
      }
      for (const auto& t : p.tubes) {
        sum += t.fluctuation;
        ++count;
      }
    }
    check::near(sum / count, 1.0, 5.0 * std::sqrt(0.5 / count), "fluctuation mean is 1");
  }

  return check::result("test_proton");
}
