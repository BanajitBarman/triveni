// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Banajit Barman

#include "triveni/VectorMeson.h"

#include "triveni/Bessel.h"

#include <cmath>

namespace triveni
{

namespace
{

constexpr double kPi = 3.14159265358979323846;
constexpr double kNc = 3.0;
constexpr double kAlphaEm = 1.0 / 137.035999;

} // namespace

double PhotonVectorMesonOverlap::phi(double norm, double r, double z) const
{
  const double zz = z * (1.0 - z);
  const double m2R2 = mWf.quarkMass * mWf.quarkMass * mWf.radius2;
  return norm * zz * std::exp(-m2R2 / (8.0 * zz) - 2.0 * zz * r * r / mWf.radius2 + m2R2 / 2.0);
}

double PhotonVectorMesonOverlap::transverse(double r, double z, double Q2) const
{
  const double zz = z * (1.0 - z);
  const double mf2 = mWf.quarkMass * mWf.quarkMass;
  const double eps = std::sqrt(zz * Q2 + mf2);
  const double p = phi(mWf.normT, r, z);
  const double dphi = -4.0 * zz * r / mWf.radius2 * p; // d phi / d r
  const double e = std::sqrt(4.0 * kPi * kAlphaEm);
  return mWf.quarkCharge * e * kNc / (kPi * zz) *
         (mf2 * besselK0(eps * r) * p - (z * z + (1.0 - z) * (1.0 - z)) * eps * besselK1(eps * r) * dphi);
}

double PhotonVectorMesonOverlap::longitudinal(double r, double z, double Q2) const
{
  if (Q2 <= 0.0) {
    return 0.0;
  }
  const double zz = z * (1.0 - z);
  const double mf2 = mWf.quarkMass * mWf.quarkMass;
  const double eps = std::sqrt(zz * Q2 + mf2);
  const double p = phi(mWf.normL, r, z);
  const double a = 2.0 * zz / mWf.radius2;                // phi ~ exp(-a r^2)
  const double laplacian = (4.0 * a * a * r * r - 4.0 * a) * p;
  const double e = std::sqrt(4.0 * kPi * kAlphaEm);
  return mWf.quarkCharge * e * kNc / kPi * 2.0 * std::sqrt(Q2) * zz * besselK0(eps * r) *
         (mWf.mass * p + (mf2 * p - laplacian) / (mWf.mass * zz));
}

double PhotonVectorMesonOverlap::normalisationT() const
{
  // 1 = N_c/(2 pi) int dz/(z^2 (1-z)^2) int d^2r { m_f^2 phi_T^2 + [z^2 + (1-z)^2] (d_r phi_T)^2 }
  const double mf2 = mWf.quarkMass * mWf.quarkMass;
  const int nz = 400;
  const int nr = 800;
  const double rMax = 20.0;
  double sum = 0.0;
  for (int iz = 0; iz < nz; ++iz) {
    const double z = (iz + 0.5) / nz;
    const double zz = z * (1.0 - z);
    double inner = 0.0;
    for (int ir = 0; ir < nr; ++ir) {
      const double r = (ir + 0.5) * rMax / nr;
      const double p = phi(mWf.normT, r, z);
      const double dphi = -4.0 * zz * r / mWf.radius2 * p;
      inner += 2.0 * kPi * r * (mf2 * p * p + (z * z + (1.0 - z) * (1.0 - z)) * dphi * dphi);
    }
    sum += inner * (rMax / nr) / (zz * zz);
  }
  return kNc / (2.0 * kPi) * sum / nz;
}

} // namespace triveni
