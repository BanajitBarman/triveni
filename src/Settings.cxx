// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Banajit Barman

#include "triveni/Settings.h"

#include "triveni/Version.h"

#include <algorithm>
#include <cctype>
#include <fstream>
#include <iomanip>
#include <ostream>
#include <stdexcept>

namespace triveni
{

namespace
{

std::string lower(std::string s)
{
  std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c) { return std::tolower(c); });
  return s;
}

std::string trim(const std::string& s)
{
  const auto begin = s.find_first_not_of(" \t\r\n");
  if (begin == std::string::npos) {
    return "";
  }
  const auto end = s.find_last_not_of(" \t\r\n");
  return s.substr(begin, end - begin + 1);
}

} // namespace

Settings::Settings()
{
  add("Main:numberOfEvents", "10000", "Number of accepted events (triveni-pp) or sampled protons (triveni-proton)");
  add("Main:threads", "1", "Worker threads; results do not depend on this number");
  add("Main:output", "triveni", "Output path stem: <stem>.txt, <stem>.log, <stem>_maps.txt, ...");
  add("Random:seed", "1", "Master random seed");

  add("Proton:shape", "Y", "Proton model: hotspot | Y | Delta | YDelta");
  add("Proton:quarkSpread", "0.36", "Gaussian sigma (fm) of each Cartesian coordinate of the valence-quark positions");
  add("Proton:quarkWidth", "0.17", "Gaussian width (fm) of the gluon cloud around each valence quark");
  add("Proton:stringWidth", "0.17", "Gaussian width (fm) of the flux tubes");
  add("Proton:stringFraction", "0.5", "Fraction of the proton thickness carried by the flux tubes (Y, Delta)");
  add("Proton:dYDelta", "0.5", "YDelta: mean quark separation (fm) below which the Delta topology is used");
  add("Proton:fluctuationK", "0", "Gamma(k, 1/k) weight fluctuation per hot spot and flux tube; <= 0 disables");

  add("Grid:halfWidth", "3.5", "Half width (fm) of the transverse grid");
  add("Grid:step", "0.05", "Cell size (fm) of the transverse grid");

  add("PP:sigmaInel", "42", "Inelastic pp cross section (mb) used to fix sigma_gg (42 mb: sqrt(s) = 200 GeV)");
  add("PP:bMax", "3.5", "Maximum impact parameter (fm)");
  add("PP:reducedP", "0", "TRENTo p: s = ((T1^p + T2^p)/2)^(1/p); p = 0 gives sqrt(T1 T2)");
  add("PP:calibrationPairs", "100000", "Proton pairs used to calibrate sigma_gg (statistical error ~0.4%)");
  add("PP:nMaps", "4", "Entropy-density maps written for the first accepted events");

  add("ProtonApp:nMaps", "3", "Quark configurations drawn as thickness maps");
  add("ProtonApp:mapHalfWidth", "1.5", "Half width (fm) of the single-proton grid");
  add("ProtonApp:mapStep", "0.02", "Cell size (fm) of the single-proton grid");
}

void Settings::add(const std::string& key, const std::string& defaultValue, const std::string& help)
{
  mEntries[lower(key)] = Entry{key, defaultValue, defaultValue, help};
}

const Settings::Entry& Settings::find(const std::string& key) const
{
  const auto it = mEntries.find(lower(key));
  if (it == mEntries.end()) {
    throw std::invalid_argument("Settings: unknown key '" + key + "'");
  }
  return it->second;
}

void Settings::readFile(const std::string& path)
{
  std::ifstream in(path);
  if (!in) {
    throw std::runtime_error("Settings: cannot open '" + path + "'");
  }
  std::string line;
  while (std::getline(in, line)) {
    readString(line);
  }
}

void Settings::readString(const std::string& rawLine)
{
  std::string line = rawLine;
  const auto comment = line.find_first_of("#!");
  if (comment != std::string::npos) {
    line.erase(comment);
  }
  line = trim(line);
  if (line.empty()) {
    return;
  }
  const auto eq = line.find('=');
  if (eq == std::string::npos) {
    throw std::invalid_argument("Settings: expected 'Key = value', got '" + rawLine + "'");
  }
  const std::string key = lower(trim(line.substr(0, eq)));
  const auto it = mEntries.find(key);
  if (it == mEntries.end()) {
    throw std::invalid_argument("Settings: unknown key '" + trim(line.substr(0, eq)) + "'");
  }
  it->second.value = trim(line.substr(eq + 1));
}

bool Settings::readArguments(int argc, const char* const* argv, std::ostream& usage)
{
  for (int i = 1; i < argc; ++i) {
    const std::string arg = argv[i];
    if (arg == "-h" || arg == "--help") {
      usage << "TRIVENI " << kVersion << "\n\nUsage: " << argv[0] << " [card.cmnd ...] [Key=value ...]\n\nSettings and defaults:\n";
      list(usage);
      return false;
    }
    if (arg.find('=') != std::string::npos) {
      readString(arg);
    } else {
      readFile(arg);
    }
  }
  return true;
}

std::string Settings::word(const std::string& key) const { return find(key).value; }

double Settings::parm(const std::string& key) const
{
  const std::string& value = find(key).value;
  std::size_t used = 0;
  double result = 0.0;
  try {
    result = std::stod(value, &used);
  } catch (const std::exception&) {
    used = 0;
  }
  if (used != value.size()) {
    throw std::invalid_argument("Settings: '" + key + "' is not a number: '" + value + "'");
  }
  return result;
}

int Settings::mode(const std::string& key) const
{
  const std::string& value = find(key).value;
  std::size_t used = 0;
  int result = 0;
  try {
    result = std::stoi(value, &used);
  } catch (const std::exception&) {
    used = 0;
  }
  if (used != value.size()) {
    throw std::invalid_argument("Settings: '" + key + "' is not an integer: '" + value + "'");
  }
  return result;
}

bool Settings::flag(const std::string& key) const
{
  const std::string value = lower(find(key).value);
  if (value == "on" || value == "true" || value == "yes" || value == "1") {
    return true;
  }
  if (value == "off" || value == "false" || value == "no" || value == "0") {
    return false;
  }
  throw std::invalid_argument("Settings: '" + key + "' is not a flag: '" + value + "'");
}

void Settings::list(std::ostream& os) const
{
  for (const auto& [unused, entry] : mEntries) {
    os << (entry.value == entry.defaultValue ? "  " : "* ") << std::left << std::setw(26) << entry.name
       << " = " << std::setw(10) << entry.value << "  ! " << entry.help << '\n';
  }
}

} // namespace triveni
