// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Banajit Barman

/// \file Bessel.h
/// \brief Bessel functions J0, K0, K1 (Abramowitz & Stegun 9.4 and 9.8).
///
/// libc++ does not ship the C++17 special functions, so they are implemented
/// here. Relative accuracy is about 1e-7, far below the other uncertainties
/// of the dipole calculations that use them.

#pragma once

namespace triveni
{

double besselJ0(double x);
double besselK0(double x); ///< x > 0
double besselK1(double x); ///< x > 0

} // namespace triveni
