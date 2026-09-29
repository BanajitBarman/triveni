// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Banajit Barman

#include "Check.h"

#include "triveni/Random.h"

#include <cmath>

int main()
{
  constexpr int n = 1000000;

  {
    triveni::Random a(42);
    triveni::Random b(42);
    bool same = true;
    for (int i = 0; i < 1000; ++i) {
      same = same && a.uniform() == b.uniform();
    }
    check::expect(same, "same seed gives the same sequence");
    triveni::Random c(triveni::Random::streamSeed(42, 0));
    triveni::Random d(triveni::Random::streamSeed(42, 1));
    check::expect(c.uniform() != d.uniform(), "different streams differ");
  }

  {
    triveni::Random rng(1);
    double sum = 0.0;
    double sum2 = 0.0;
    double lo = 1.0;
    double hi = 0.0;
    for (int i = 0; i < n; ++i) {
      const double u = rng.uniform();
      sum += u;
      sum2 += u * u;
      lo = std::min(lo, u);
      hi = std::max(hi, u);
    }
    check::near(sum / n, 0.5, 5.0 * std::sqrt(1.0 / 12.0 / n), "uniform mean");
    check::near(sum2 / n - (sum / n) * (sum / n), 1.0 / 12.0, 2e-3, "uniform variance");
    check::expect(lo > 0.0 && hi < 1.0, "uniform stays inside (0, 1)");
  }

  {
    triveni::Random rng(2);
    double sum = 0.0;
    double sum2 = 0.0;
    double sum4 = 0.0;
    for (int i = 0; i < n; ++i) {
      const double x = rng.normal();
      sum += x;
      sum2 += x * x;
      sum4 += x * x * x * x;
    }
    check::near(sum / n, 0.0, 5.0 / std::sqrt(n), "normal mean");
    check::near(sum2 / n, 1.0, 5.0 * std::sqrt(2.0 / n), "normal variance");
    check::near(sum4 / n, 3.0, 5.0 * std::sqrt(96.0 / n), "normal fourth moment");
  }

  for (double k : {0.4, 1.0, 2.5, 10.0}) {
    triveni::Random rng(3);
    double sum = 0.0;
    double sum2 = 0.0;
    bool positive = true;
    for (int i = 0; i < n; ++i) {
      const double g = rng.gamma(k);
      positive = positive && g > 0.0;
      sum += g;
      sum2 += g * g;
    }
    const double mean = sum / n;
    const double variance = sum2 / n - mean * mean;
    const std::string label = "gamma(k=" + std::to_string(k) + ")";
    check::expect(positive, label + " is positive");
    check::near(mean, 1.0, 5.0 * std::sqrt(1.0 / k / n), label + " mean");
    check::near(variance, 1.0 / k, 0.02 / k + 0.02, label + " variance");
  }

  return check::result("test_random");
}
