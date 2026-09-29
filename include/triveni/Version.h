// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Banajit Barman

/// \file Version.h
/// \brief TRIVENI version, set by CMake from project(VERSION ...). Every log records it.

#pragma once

#ifndef TRIVENI_VERSION
#define TRIVENI_VERSION "unknown"
#endif

namespace triveni
{

inline constexpr const char* kVersion = TRIVENI_VERSION;

} // namespace triveni
