// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Banajit Barman

#include "triveni/Exclusive.h"

#include "triveni/Bessel.h"
#include "triveni/Parallel.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace triveni
{

namespace
{

constexpr double kPi = 3.14159265358979323846;
constexpr double kEulerGamma = 0.57721566490153286;

/// Simpson weights for n points (n odd) spaced by h.
std::vector<double> simpsonWeights(int n, double h)
{
  if (n < 3 || n % 2 == 0) {
    throw std::invalid_argument("simpsonWeights: need an odd number of points >= 3");
  }
  std::vector<double> w(n);
  for (int i = 0; i < n; ++i) {
    w[i] = h / 3.0 * ((i == 0 || i == n - 1) ? 1.0 : (i % 2 ? 4.0 : 2.0));
  }
  return w;
}

/// Gauss–Legendre nodes and weights on [0, 1].
void gaussLegendre(int n, std::vector<double>& x, std::vector<double>& w)
{
  x.resize(n);
  w.resize(n);
  for (int i = 0; i < n; ++i) {
    double z = std::cos(kPi * (i + 0.75) / (n + 0.5));
    double dp = 0.0;
    for (int iter = 0; iter < 100; ++iter) {
      double p0 = 1.0;
      double p1 = z;
      for (int k = 2; k <= n; ++k) {
        const double p2 = ((2.0 * k - 1.0) * z * p1 - (k - 1.0) * p0) / k;
        p0 = p1;
        p1 = p2;
      }
      dp = n * (z * p1 - p0) / (z * z - 1.0);
      const double dz = p1 / dp;
      z -= dz;
      if (std::fabs(dz) < 1e-15) {
        break;
      }
    }
    x[i] = 0.5 * (1.0 - z);
    w[i] = 1.0 / ((1.0 - z * z) * dp * dp);
  }
}

/// Ein(c) = int_0^c (1 - e^-t)/t dt = E1(c) + ln c + gamma.
double ein(double c)
{
  if (c < 8.0) {
    double term = c;
    double sum = c;
    for (int k = 2; k < 200; ++k) {
      term *= -c / k;
      const double add = term / k;
      sum += add;
      if (std::fabs(add) < 1e-16 * std::fabs(sum)) {
        break;
      }
    }
    return sum;
  }
  // E1(c) by its continued fraction (modified Lentz), accurate for c >= 1.
  double b = c + 1.0;
  double cc = 1e300;
  double d = 1.0 / b;
  double h = d;
  for (int i = 1; i < 300; ++i) {
    const double a = -static_cast<double>(i) * i;
    b += 2.0;
    d = 1.0 / (a * d + b);
    cc = b + a / cc;
    const double del = cc * d;
    h *= del;
    if (std::fabs(del - 1.0) < 1e-15) {
      break;
    }
  }
  return h * std::exp(-c) + std::log(c) + kEulerGamma;
}

} // namespace

CoherentVectorMeson::CoherentVectorMeson(const VectorMesonWavefunction& wf, const DipoleParams& dipole,
                                         std::optional<WoodsSaxon> nucleus)
  : mOverlap(wf), mDipole(dipole)
{
  if (nucleus) {
    mNucleus.emplace(*nucleus);
  }

  // r: log grid, 1e-3 .. 25 GeV^-1 (the J/psi overlap dies like exp(-m_c r)).
  {
    const int n = 161;
    const double u0 = std::log(1e-3);
    const double u1 = std::log(25.0);
    const double h = (u1 - u0) / (n - 1);
    const auto w = simpsonWeights(n, h);
    for (int i = 0; i < n; ++i) {
      const double r = std::exp(u0 + i * h);
      mR.push_back(r);
      mRWeight.push_back(w[i] * r * 2.0 * kPi * r); // dr = r du, d^2r = 2 pi r dr
    }
  }
  gaussLegendre(48, mZ, mZWeight);

  // b grid and target thickness; Delta grid for the t integral.
  int nb = 0;
  double bMax = 0.0;
  double deltaMax = 0.0;
  int nDelta = 0;
  if (mNucleus) {
    bMax = mNucleus->maxImpactParameter() / constants::hbarc;
    nb = 1601;
    deltaMax = 0.5; // |t| < 0.25 GeV^2: covers the coherent peak and several dips
    nDelta = 161;
  } else {
    bMax = 12.0 * std::sqrt(mDipole.protonB);
    nb = 601;
    deltaMax = 2.5;
    nDelta = 101;
  }
  {
    const double h = bMax / (nb - 1);
    const auto w = simpsonWeights(nb, h);
    for (int j = 0; j < nb; ++j) {
      const double b = j * h;
      mB.push_back(b);
      mBWeight.push_back(w[j] * 2.0 * kPi * b);
      if (mNucleus) {
        mThickness.push_back(mNucleus->thickness(b * constants::hbarc) * constants::hbarc * constants::hbarc);
      } else {
        mThickness.push_back(std::exp(-b * b / (2.0 * mDipole.protonB)) / (2.0 * kPi * mDipole.protonB));
      }
    }
  }
  {
    const double h = deltaMax / (nDelta - 1);
    const auto w = simpsonWeights(nDelta, h);
    for (int k = 0; k < nDelta; ++k) {
      const double d = k * h;
      mDelta.push_back(d);
      mAbsT.push_back(d * d);
      mTWeight.push_back(w[k] * 2.0 * d); // d|t| = 2 Delta dDelta
    }
  }
  mJ0.resize(mDelta.size() * mB.size());
  for (std::size_t k = 0; k < mDelta.size(); ++k) {
    for (std::size_t j = 0; j < mB.size(); ++j) {
      mJ0[k * mB.size() + j] = besselJ0(mB[j] * mDelta[k]);
    }
  }
}

double CoherentVectorMeson::xPomeron(double W, double Q2, double mass)
{
  return (mass * mass + Q2) / (W * W + Q2);
}

double CoherentVectorMeson::protonDipoleCrossSection(double r, double x) const
{
  // 2 int d^2b [1 - exp(-r^2 Omega T_p(b))] = 4 pi B_p Ein(r^2 Omega / (2 pi B_p)) for a Gaussian T_p.
  const double omega = mDipole.omega0 * std::pow(x, -mDipole.lambda);
  const double B = mDipole.protonB;
  return 4.0 * kPi * B * ein(r * r * omega / (2.0 * kPi * B));
}

void CoherentVectorMeson::dsdtOnGrid(double W, double Q2, std::vector<double>& dsdtT, std::vector<double>& dsdtL) const
{
  const double mass = vectorMesonMass();
  const double x = xPomeron(W, Q2, mass);
  const double omega = mDipole.omega0 * std::pow(x, -mDipole.lambda);
  const std::size_t nr = mR.size();
  const std::size_t nb = mB.size();
  const std::size_t nd = mDelta.size();

  std::vector<double> aT(nd, 0.0);
  std::vector<double> aL(nd, 0.0);
  std::vector<double> twoN(nb);
  std::vector<double> psiT(mZ.size());
  std::vector<double> psiL(mZ.size());
  for (std::size_t i = 0; i < nr; ++i) {
    const double r = mR[i];
    // Dipole amplitude on the b grid: proton directly, nucleus through sigma_p(r).
    const double strength = mNucleus ? 0.5 * protonDipoleCrossSection(r, x) : r * r * omega;
    for (std::size_t j = 0; j < nb; ++j) {
      twoN[j] = 2.0 * (-std::expm1(-strength * mThickness[j]));
    }
    for (std::size_t iz = 0; iz < mZ.size(); ++iz) {
      psiT[iz] = mOverlap.transverse(r, mZ[iz], Q2) * mZWeight[iz] / (4.0 * kPi);
      psiL[iz] = mOverlap.longitudinal(r, mZ[iz], Q2) * mZWeight[iz] / (4.0 * kPi);
    }
    for (std::size_t k = 0; k < nd; ++k) {
      const double* j0 = &mJ0[k * nb];
      double g = 0.0;
      for (std::size_t j = 0; j < nb; ++j) {
        g += mBWeight[j] * j0[j] * twoN[j];
      }
      double iT = 0.0;
      double iL = 0.0;
      for (std::size_t iz = 0; iz < mZ.size(); ++iz) {
        const double phase = besselJ0((1.0 - mZ[iz]) * r * mDelta[k]);
        iT += psiT[iz] * phase;
        iL += psiL[iz] * phase;
      }
      aT[k] += mRWeight[i] * iT * g;
      aL[k] += mRWeight[i] * iL * g;
    }
  }
  dsdtT.resize(nd);
  dsdtL.resize(nd);
  for (std::size_t k = 0; k < nd; ++k) {
    dsdtT[k] = aT[k] * aT[k] / (16.0 * kPi) * constants::gev2ToNb;
    dsdtL[k] = aL[k] * aL[k] / (16.0 * kPi) * constants::gev2ToNb;
  }
}

CoherentVectorMeson::Sigma CoherentVectorMeson::sigma(double W, double Q2) const
{
  std::vector<double> t;
  std::vector<double> l;
  dsdtOnGrid(W, Q2, t, l);
  Sigma s;
  for (std::size_t k = 0; k < t.size(); ++k) {
    s.T += mTWeight[k] * t[k];
    s.L += mTWeight[k] * l[k];
  }
  return s;
}

DipoleParams fitDipoleToPhotoproduction(const VectorMesonWavefunction& wf, double protonB,
                                        const std::function<double(double)>& target, const std::vector<double>& wPoints)
{
  if (wPoints.size() < 3) {
    throw std::invalid_argument("fitDipoleToPhotoproduction: need at least 3 W points");
  }
  const double wLo = wPoints.front();
  const double wMid = wPoints.at(wPoints.size() / 2);
  const double wHi = wPoints.back();
  DipoleParams p;
  p.protonB = protonB;
  CoherentVectorMeson model(wf, p);

  // For a given lambda, Omega_0 such that sigma(wMid) = target(wMid); sigma rises with Omega_0.
  auto normalise = [&](double lambda) {
    double lo = 1e-6;
    double hi = 1e3;
    for (int it = 0; it < 80; ++it) {
      const double mid = std::sqrt(lo * hi);
      model.setDipole({mid, lambda, protonB});
      (model.sigma(wMid, 0.0).T < target(wMid) ? lo : hi) = mid;
    }
    return std::sqrt(lo * hi);
  };
  // lambda such that the ratio sigma(wHi)/sigma(wLo) matches the target's.
  const double wantRatio = target(wHi) / target(wLo);
  double lo = 0.0;
  double hi = 0.6;
  for (int it = 0; it < 40; ++it) {
    const double lambda = 0.5 * (lo + hi);
    const double omega0 = normalise(lambda);
    model.setDipole({omega0, lambda, protonB});
    const double ratio = model.sigma(wHi, 0.0).T / model.sigma(wLo, 0.0).T;
    (ratio < wantRatio ? lo : hi) = lambda;
  }
  p.lambda = 0.5 * (lo + hi);
  p.omega0 = normalise(p.lambda);
  return p;
}

double electronPhotonFlux(double E, double Q2, double electronEnergy)
{
  const double y = E / electronEnergy;
  if (!(y > 0.0 && y < 1.0) || !(Q2 > 0.0)) {
    return 0.0;
  }
  const double q2Min = std::pow(constants::electronMass * E, 2) / (electronEnergy * (electronEnergy - E));
  if (Q2 < q2Min) {
    return 0.0;
  }
  return constants::alphaEm / kPi * (1.0 - y + 0.5 * y * y - (1.0 - y) * q2Min / Q2) / (E * Q2);
}

double photonPolarisation(double E, double electronEnergy)
{
  const double y = E / electronEnergy;
  return (1.0 - y) / (1.0 - y + 0.5 * y * y);
}

ElectronCrossSection electronTargetCrossSection(const CoherentVectorMeson& model, double electronEnergy, double wMin,
                                                double wMax, double q2Min, double q2Max, int threads)
{
  const double mN = constants::protonMass;
  // Largest W reachable: photon energy up to (almost) the electron energy.
  const double wKin = std::sqrt(2.0 * mN * electronEnergy * (1.0 - 1e-6) + mN * mN);
  wMax = std::min(wMax, wKin);
  ElectronCrossSection out;
  out.absT = model.tGrid();
  out.dsdt.assign(out.absT.size(), 0.0);
  if (!(wMax > wMin) || !(q2Max > q2Min)) {
    return out;
  }

  // Integrate in ln W (outer) and ln Q^2 (inner); E = (W^2 + Q^2 - m^2)/(2m), dE = W dW / m.
  const int nW = 33;
  const int nQ = 41;
  const double lw0 = std::log(wMin);
  const double hW = (std::log(wMax) - lw0) / (nW - 1);
  const auto wW = simpsonWeights(nW, hW);
  std::vector<std::vector<double>> rows(nW, std::vector<double>(out.absT.size(), 0.0));

  parallelBlocks(nW, threads, [&](int iw) {
    const double W = std::exp(lw0 + iw * hW);
    // Q^2 range at this W (lower edge: the kinematic minimum for photoproduction windows).
    const double e0 = (W * W - mN * mN) / (2.0 * mN);
    const double y0 = e0 / electronEnergy;
    const double q2Kin = std::pow(constants::electronMass * y0, 2) / (1.0 - y0);
    const double qLo = std::max(q2Min, q2Kin * 1.0000001);
    const double qHi = std::min(q2Max, 4.0 * electronEnergy * (electronEnergy - e0));
    if (!(qHi > qLo)) {
      return;
    }
    const double lq0 = std::log(qLo);
    const double hQ = (std::log(qHi) - lq0) / (nQ - 1);
    const auto wQ = simpsonWeights(nQ, hQ);
    std::vector<double> dT;
    std::vector<double> dL;
    for (int iq = 0; iq < nQ; ++iq) {
      const double Q2 = std::exp(lq0 + iq * hQ);
      const double E = (W * W + Q2 - mN * mN) / (2.0 * mN);
      const double flux = electronPhotonFlux(E, Q2, electronEnergy);
      if (flux <= 0.0) {
        continue;
      }
      const double eps = photonPolarisation(E, electronEnergy);
      model.dsdtOnGrid(W, Q2, dT, dL);
      // dE dQ^2 = (W/m) dW dQ^2 = (W/m) W Q^2 dlnW dlnQ^2
      const double jac = wW[iw] * W * (W / mN) * wQ[iq] * Q2 * flux;
      for (std::size_t k = 0; k < dT.size(); ++k) {
        rows[iw][k] += jac * (dT[k] + eps * dL[k]);
      }
    }
  });

  for (const auto& row : rows) {
    for (std::size_t k = 0; k < row.size(); ++k) {
      out.dsdt[k] += row[k];
    }
  }
  // Total: integrate the flux-weighted spectrum over |t| with the model's own weights.
  const auto& tw = model.tWeights();
  for (std::size_t k = 0; k < out.dsdt.size(); ++k) {
    out.sigma += tw[k] * out.dsdt[k];
  }
  return out;
}

} // namespace triveni
