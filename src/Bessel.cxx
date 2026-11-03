// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Banajit Barman

#include "triveni/Bessel.h"

#include <cmath>
#include <stdexcept>

namespace triveni
{

namespace
{

// Modified Bessel functions I0 and I1 for |x| <= 3.75 (A&S 9.8.1, 9.8.3).
double besselI0Small(double x)
{
  const double t = x / 3.75;
  const double t2 = t * t;
  return 1.0 + t2 * (3.5156229 + t2 * (3.0899424 + t2 * (1.2067492 + t2 * (0.2659732 + t2 * (0.0360768 + t2 * 0.0045813)))));
}

double besselI1Small(double x)
{
  const double t = x / 3.75;
  const double t2 = t * t;
  return x * (0.5 + t2 * (0.87890594 + t2 * (0.51498869 + t2 * (0.15084934 + t2 * (0.02658733 + t2 * (0.00301532 + t2 * 0.00032411))))));
}

} // namespace

double besselJ0(double x)
{
  x = std::fabs(x);
  if (x <= 3.0) { // A&S 9.4.1
    const double y = (x / 3.0) * (x / 3.0);
    return 1.0 + y * (-2.2499997 + y * (1.2656208 + y * (-0.3163866 + y * (0.0444479 + y * (-0.0039444 + y * 0.0002100)))));
  }
  // A&S 9.4.3
  const double y = 3.0 / x;
  const double f0 = 0.79788456 + y * (-0.00000077 + y * (-0.00552740 + y * (-0.00009512 + y * (0.00137237 + y * (-0.00072805 + y * 0.00014476)))));
  const double theta0 = x - 0.78539816 +
                        y * (-0.04166397 + y * (-0.00003954 + y * (0.00262573 + y * (-0.00054125 + y * (-0.00029333 + y * 0.00013558)))));
  return f0 * std::cos(theta0) / std::sqrt(x);
}

double besselK0(double x)
{
  if (!(x > 0.0)) {
    throw std::domain_error("besselK0: x must be positive");
  }
  if (x <= 2.0) { // A&S 9.8.5
    const double y = (x / 2.0) * (x / 2.0);
    return -std::log(x / 2.0) * besselI0Small(x) +
           (-0.57721566 + y * (0.42278420 + y * (0.23069756 + y * (0.03488590 + y * (0.00262698 + y * (0.00010750 + y * 0.00000740))))));
  }
  // A&S 9.8.6
  const double y = 2.0 / x;
  return std::exp(-x) / std::sqrt(x) *
         (1.25331414 + y * (-0.07832358 + y * (0.02189568 + y * (-0.01062446 + y * (0.00587872 + y * (-0.00251540 + y * 0.00053208))))));
}

double besselK1(double x)
{
  if (!(x > 0.0)) {
    throw std::domain_error("besselK1: x must be positive");
  }
  if (x <= 2.0) { // A&S 9.8.7
    const double y = (x / 2.0) * (x / 2.0);
    return (x * std::log(x / 2.0) * besselI1Small(x) +
            (1.0 + y * (0.15443144 + y * (-0.67278579 + y * (-0.18156897 + y * (-0.01919402 + y * (-0.00110404 + y * -0.00004686))))))) /
           x;
  }
  // A&S 9.8.8
  const double y = 2.0 / x;
  return std::exp(-x) / std::sqrt(x) *
         (1.25331414 + y * (0.23498619 + y * (-0.03655620 + y * (0.01504268 + y * (-0.00780353 + y * (0.00325614 + y * -0.00068245))))));
}

} // namespace triveni
