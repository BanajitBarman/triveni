// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Banajit Barman

/// \file Settings.h
/// \brief PYTHIA-style "Key = value" settings with registered defaults.
///
/// Every key must be registered (with a default and a description) so that a
/// typo in a card file is an error instead of a silently ignored line.
/// Keys are case-insensitive.

#pragma once

#include <iosfwd>
#include <map>
#include <string>

namespace triveni
{

class Settings
{
 public:
  /// Registers all known keys with their defaults.
  Settings();

  /// Reads a card file: one "Key = value" per line; '#' or '!' start a comment.
  void readFile(const std::string& path);

  /// Parses a single "Key = value" (or "Key=value") line.
  void readString(const std::string& line);

  /// Command line: arguments containing '=' are settings, others are card files,
  /// applied in order. Returns false (after printing usage) for -h / --help.
  bool readArguments(int argc, const char* const* argv, std::ostream& usage);

  std::string word(const std::string& key) const;
  double parm(const std::string& key) const;
  int mode(const std::string& key) const;
  bool flag(const std::string& key) const;

  /// Prints all settings, marking the ones changed from their default.
  void list(std::ostream& os) const;

 private:
  struct Entry {
    std::string name; ///< Key with its original capitalisation.
    std::string value;
    std::string defaultValue;
    std::string help;
  };

  void add(const std::string& key, const std::string& defaultValue, const std::string& help);
  const Entry& find(const std::string& key) const;

  std::map<std::string, Entry> mEntries;
};

} // namespace triveni
