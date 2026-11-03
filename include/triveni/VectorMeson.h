// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Banajit Barman

/// \file VectorMeson.h
/// \brief Overlap of the virtual-photon and vector-meson light-cone wavefunctions.
///
/// Boosted-Gaussian vector meson, Kowalski, Motyka and Watt, PRD 74 (2006) 074016:
///   phi_{T,L}(r,z) = N_{T,L} z(1-z) exp[-m_f^2 R^2/(8z(1-z)) - 2z(1-z) r^2/R^2 + m_f^2 R^2/2]
///   (Psi_V^* Psi)_T = e_f e N_c/(pi z(1-z)) [ m_f^2 K0(eps r) phi_T - (z^2+(1-z)^2) eps K1(eps r) d_r phi_T ]
///   (Psi_V^* Psi)_L = e_f e N_c/pi 2Q z(1-z) K0(eps r) [ M_V phi_L + (m_f^2 - nabla_r^2) phi_L / (M_V z(1-z)) ]
/// with eps^2 = z(1-z)Q^2 + m_f^2. Units: GeV, r in GeV^-1.

#pragma once

namespace triveni
{

struct VectorMesonWavefunction {
  double mass = 3.097;      ///< M_V (GeV)
  double quarkMass = 1.4;   ///< m_f (GeV)
  double quarkCharge = 2.0 / 3.0;
  double normT = 0.578;
  double normL = 0.575;
  double radius2 = 2.3;     ///< R^2 (GeV^-2)

  /// J/psi, boosted Gaussian, m_c = 1.4 GeV (KMW 2006, Table 2).
  static VectorMesonWavefunction jpsi() { return {}; }
};

class PhotonVectorMesonOverlap
{
 public:
  explicit PhotonVectorMesonOverlap(const VectorMesonWavefunction& wf) : mWf(wf) {}

  const VectorMesonWavefunction& wavefunction() const { return mWf; }

  /// (Psi_V^* Psi)_T, summed over helicities and flavours (GeV^2).
  double transverse(double r, double z, double Q2) const;
  /// (Psi_V^* Psi)_L (GeV^2).
  double longitudinal(double r, double z, double Q2) const;

  /// Normalisation integral of the transverse meson wavefunction; equals 1 for consistent parameters.
  double normalisationT() const;

 private:
  double phi(double norm, double r, double z) const;

  VectorMesonWavefunction mWf;
};

} // namespace triveni
