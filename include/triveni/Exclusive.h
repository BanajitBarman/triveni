// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Banajit Barman

/// \file Exclusive.h
/// \brief Coherent exclusive vector-meson production, gamma* + target -> V + target,
///        in the colour-dipole picture (Good–Walker: coherent = |<A>|^2).
///
///   A_{T,L}(Delta) = int d^2r int dz/(4 pi) int d^2b (Psi_V^* Psi)_{T,L} e^{-i[b-(1-z)r].Delta} 2 N(r,b,x)
///   dsigma/dt      = |A|^2 / (16 pi),    t = -Delta^2
///
/// Proton:  N_p(r,b,x) = 1 - exp[-r^2 Omega(x) T_p(b)],  T_p Gaussian of width B_p,
///          Omega(x) = Omega_0 x^-lambda. An eikonalised, impact-parameter dependent
///          saturation model without DGLAP evolution; Omega_0 and lambda are fitted
///          to the gamma p -> V p cross section.
/// Nucleus: <N_A(r,b,x)> = 1 - exp[-sigma_p(r,x) T_A(b) / 2], sigma_p = 2 int d^2b N_p,
///          with a Woods–Saxon T_A: the coherent average with no sub-nucleon structure.
/// x = (M_V^2 + Q^2)/(W^2 + Q^2). Real-part and skewness corrections are not included.

#pragma once

#include "triveni/Nucleus.h"
#include "triveni/VectorMeson.h"

#include <functional>
#include <optional>
#include <vector>

namespace triveni
{

namespace constants
{
constexpr double hbarc = 0.1973269804;       ///< GeV fm
constexpr double gev2ToNb = 0.3893793721e6;  ///< 1 GeV^-2 in nb
constexpr double protonMass = 0.938272046;   ///< GeV (as eSTARlight)
constexpr double electronMass = 0.000510998928;
constexpr double jpsiMass = 3.096916;
constexpr double alphaEm = 1.0 / 137.035999074;
} // namespace constants

struct DipoleParams {
  double omega0 = 0.1;  ///< Omega_0 (dimensionless)
  double lambda = 0.2;  ///< x^-lambda growth
  double protonB = 4.0; ///< B_p (GeV^-2): T_p(b) = exp(-b^2/2B_p)/(2 pi B_p)
};

class CoherentVectorMeson
{
 public:
  /// Proton target when `nucleus` is empty.
  CoherentVectorMeson(const VectorMesonWavefunction& wf, const DipoleParams& dipole,
                      std::optional<WoodsSaxon> nucleus = std::nullopt);

  void setDipole(const DipoleParams& dipole) { mDipole = dipole; }
  const DipoleParams& dipole() const { return mDipole; }
  bool isNucleus() const { return mNucleus.has_value(); }
  double vectorMesonMass() const { return mOverlap.wavefunction().mass; }

  static double xPomeron(double W, double Q2, double mass);

  /// |t| values (GeV^2) of the internal grid used by dsdtOnGrid() and sigma(),
  /// and the weights that integrate a function of |t| over it.
  const std::vector<double>& tGrid() const { return mAbsT; }
  const std::vector<double>& tWeights() const { return mTWeight; }

  /// Transverse and longitudinal dsigma/dt (nb/GeV^2) on tGrid().
  void dsdtOnGrid(double W, double Q2, std::vector<double>& dsdtT, std::vector<double>& dsdtL) const;

  struct Sigma {
    double T = 0.0; ///< nb
    double L = 0.0; ///< nb
  };
  /// t-integrated cross section.
  Sigma sigma(double W, double Q2) const;

  /// Proton dipole cross section sigma_p(r, x) = 2 int d^2b N_p (GeV^-2).
  double protonDipoleCrossSection(double r, double x) const;

 private:
  PhotonVectorMesonOverlap mOverlap;
  DipoleParams mDipole;
  std::optional<Nucleus> mNucleus;

  std::vector<double> mR, mRWeight;     // r grid (GeV^-1); weights include d^2r = 2 pi r dr
  std::vector<double> mZ, mZWeight;     // Gauss–Legendre on [0, 1]
  std::vector<double> mB, mBWeight;     // b grid (GeV^-1); weights include d^2b = 2 pi b db
  std::vector<double> mThickness;       // T_p(b) or T_A(b), GeV^2
  std::vector<double> mDelta, mAbsT, mTWeight; // Delta grid; weights for int d|t|
  std::vector<double> mJ0;              // J0(b_j Delta_k), row k
};

/// Fits Omega_0 and lambda so that the proton cross section at Q^2 = 0 follows
/// `target(W)` (nb): exact at the middle W point, and the same ratio between the
/// first and last W points.
DipoleParams fitDipoleToPhotoproduction(const VectorMesonWavefunction& wf, double protonB,
                                        const std::function<double(double)>& target, const std::vector<double>& wPoints);

/// eSTARlight's electron photon flux in the target rest frame:
/// dN/(dE dQ^2) = alpha/pi [1 - y + y^2/2 - (1-y) Q2min/Q^2] / (E Q^2), y = E/E_e.
double electronPhotonFlux(double E, double Q2, double electronEnergy);

/// Longitudinal-to-transverse flux ratio epsilon = (1-y)/(1 - y + y^2/2).
double photonPolarisation(double E, double electronEnergy);

struct ElectronCrossSection {
  double sigma = 0.0;           ///< nb
  std::vector<double> absT;     ///< |t| grid (GeV^2)
  std::vector<double> dsdt;     ///< flux-weighted dsigma/d|t| (nb/GeV^2)
};

/// e + target -> e + V + target over W in [wMin, wMax] and Q^2 in [q2Min, q2Max]:
/// int dE dQ^2 flux(E,Q^2) [sigma_T + eps sigma_L](W, Q^2), W^2 = 2 m_N E + m_N^2 - Q^2,
/// with E the photon energy and electronEnergy the electron energy in the target rest frame.
ElectronCrossSection electronTargetCrossSection(const CoherentVectorMeson& model, double electronEnergy, double wMin,
                                                double wMax, double q2Min, double q2Max, int threads);

} // namespace triveni
