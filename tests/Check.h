// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Banajit Barman

/// \file Check.h
/// \brief Tiny dependency-free test helpers.

#pragma once

#include <cmath>
#include <iostream>
#include <sstream>
#include <string>

namespace check
{

inline int& failures()
{
  static int count = 0;
  return count;
}

inline void expect(bool ok, const std::string& what)
{
  if (!ok) {
    ++failures();
    std::cerr << "FAIL: " << what << '\n';
  }
}

inline void near(double value, double expected, double tolerance, const std::string& what)
{
  std::ostringstream os;
  os.precision(10);
  os << what << ": got " << value << ", expected " << expected << " +- " << tolerance;
  expect(std::abs(value - expected) <= tolerance, os.str());
}

inline int result(const char* name)
{
  if (failures() > 0) {
    std::cerr << name << ": " << failures() << " failure(s)\n";
    return 1;
  }
  std::cout << name << ": all checks passed\n";
  return 0;
}

} // namespace check
