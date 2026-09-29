// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Banajit Barman

#include "triveni/PPCollision.h"

#include "triveni/Parallel.h"
#include "triveni/Settings.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>

namespace triveni
{

namespace
{

constexpr double kPi = 3.14159265358979323846;
constexpr double kMbPerFm2 = 10.0;
constexpr int kCalibrationBlock = 250;

} // namespace

PPParams PPParams::fromSettings(const Settings& s)
{
  PPParams p;
  p.proton = ProtonParams::fromSettings(s);
  p.sigmaInel = s.parm("PP:sigmaInel");
  p.bMax = s.parm("PP:bMax");
  p.reducedP = s.parm("PP:reducedP");
  p.gridHalfWidth = s.parm("Grid:halfWidth");
  p.gridStep = s.parm("Grid:step");
  if (!(p.sigmaInel > 0.0) || !(p.bMax > 0.0)) {
    throw std::invalid_argument("PP: sigmaInel and bMax must be positive");
  }
  if (p.gridHalfWidth < 0.5 * p.bMax + 1.5) {
    throw std::invalid_argument("PP: Grid:halfWidth must be at least bMax/2 + 1.5 fm to contain both protons");
  }
  return p;
}

double reducedThickness(double t1, double t2, double p)
{
  if (std::abs(p) < 1e-12) {
    return std::sqrt(t1 * t2);
  }
  if (t1 <= 0.0 || t2 <= 0.0) {
    return p > 0.0 ? std::pow(0.5 * (std::pow(std::max(t1, 0.0), p) + std::pow(std::max(t2, 0.0), p)), 1.0 / p) : 0.0;
  }
  return std::pow(0.5 * (std::pow(t1, p) + std::pow(t2, p)), 1.0 / p);
}

PPCollider::Workspace::Workspace(const PPParams& params)
  : t1(params.gridHalfWidth, params.gridStep), t2(params.gridHalfWidth, params.gridStep)
{
}

PPCollider::PPCollider(const PPParams& params) : mParams(params), mBuilder(params.proton) {}

double PPCollider::sampleOverlap(Random& rng, Workspace& ws, double b, Proton& p1, Proton& p2) const
{
  p1 = mBuilder.generate(rng);
  p2 = mBuilder.generate(rng);
  ws.t1.clear();
  ws.t2.clear();
  ws.t1.addProton(p1, {-0.5 * b, 0.0}, false);
  ws.t2.addProton(p2, {+0.5 * b, 0.0}, false);
  const auto& v1 = ws.t1.values();
  const auto& v2 = ws.t2.values();
  double sum = 0.0;
  for (std::size_t i = 0; i < v1.size(); ++i) {
    sum += v1[i] * v2[i];
  }
  return sum * ws.t1.cellArea();
}

Calibration PPCollider::calibrate(std::uint64_t seed, int nPairs, int nThreads)
{
  if (nPairs < 100) {
    throw std::invalid_argument("PPCollider::calibrate: need at least 100 pairs");
  }
  std::vector<double> impact(nPairs);
  std::vector<double> overlap(nPairs);
  const int nBlocks = (nPairs + kCalibrationBlock - 1) / kCalibrationBlock;
  parallelBlocks(nBlocks, nThreads, [&](int block) {
    Random rng(Random::streamSeed(seed, static_cast<std::uint64_t>(block)));
    Workspace ws(mParams);
    Proton p1;
    Proton p2;
    const int end = std::min(nPairs, (block + 1) * kCalibrationBlock);
    for (int i = block * kCalibrationBlock; i < end; ++i) {
      impact[i] = mParams.bMax * std::sqrt(rng.uniform());
      overlap[i] = sampleOverlap(rng, ws, impact[i], p1, p2);
    }
  });

  const double area = kPi * mParams.bMax * mParams.bMax; // fm^2
  const double target = mParams.sigmaInel / kMbPerFm2;   // fm^2
  auto sigmaInelOf = [&](double sigmaGG) {
    double sum = 0.0;
    for (double t : overlap) {
      sum -= std::expm1(-sigmaGG * t);
    }
    return area * sum / nPairs;
  };

  double lo = 1e-6;
  double hi = 1e6;
  if (sigmaInelOf(hi) <= target) {
    throw std::runtime_error("PPCollider::calibrate: sigmaInel cannot be reached; increase PP:bMax");
  }
  for (int iter = 0; iter < 200 && hi / lo > 1.0 + 1e-13; ++iter) {
    const double mid = std::sqrt(lo * hi);
    (sigmaInelOf(mid) < target ? lo : hi) = mid;
  }
  mSigmaGG = std::sqrt(lo * hi);

  Calibration c;
  c.sigmaGG = mSigmaGG;
  c.pairs = nPairs;
  double sumP = 0.0;
  double sumP2 = 0.0;
  double sumSlope = 0.0; // sum of sigma_gg * dP/dsigma_gg
  double edgeSum = 0.0;
  int edgeCount = 0;
  for (int i = 0; i < nPairs; ++i) {
    const double x = mSigmaGG * overlap[i];
    const double p = -std::expm1(-x);
    sumP += p;
    sumP2 += p * p;
    sumSlope += x * std::exp(-x);
    if (impact[i] > 0.9 * mParams.bMax) {
      edgeSum += p;
      ++edgeCount;
    }
  }
  const double mean = sumP / nPairs;
  const double variance = std::max(0.0, sumP2 / nPairs - mean * mean);
  c.sigmaInel = area * mean * kMbPerFm2;
  c.sigmaInelError = area * std::sqrt(variance / nPairs) * kMbPerFm2;
  c.edgeProbability = edgeCount > 0 ? edgeSum / edgeCount : 0.0;
  c.sensitivity = sumP > 0.0 ? sumSlope / sumP : 0.0;
  return c;
}

bool PPCollider::tryEvent(Random& rng, Workspace& ws, PPEvent& event) const
{
  if (!(mSigmaGG > 0.0)) {
    throw std::logic_error("PPCollider::tryEvent: call calibrate() or setSigmaGG() first");
  }
  const double b = mParams.bMax * std::sqrt(rng.uniform());
  const double overlap = sampleOverlap(rng, ws, b, event.p1, event.p2);
  const double pCollision = -std::expm1(-mSigmaGG * overlap);
  if (rng.uniform() >= pCollision) {
    return false;
  }

  if (mParams.proton.fluctuationK > 0.0) {
    ws.t1.clear();
    ws.t2.clear();
    ws.t1.addProton(event.p1, {-0.5 * b, 0.0}, true);
    ws.t2.addProton(event.p2, {+0.5 * b, 0.0}, true);
  }
  auto& s = ws.t1.values();
  const auto& t2 = ws.t2.values();
  for (std::size_t i = 0; i < s.size(); ++i) {
    s[i] = reducedThickness(s[i], t2[i], mParams.reducedP);
  }

  event.b = b;
  event.overlap = overlap;
  event.pCollision = pCollision;
  event.entropy = ws.t1.integral();
  event.ecc = eccentricities(ws.t1);
  return true;
}

} // namespace triveni
