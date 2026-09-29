#!/usr/bin/env bash
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Banajit Barman
# Reproduces the Phase 1 outputs and figures: build, test, single-proton study,
# pp at 200 GeV for the three proton shapes, ROOT figures in output/figures/.
#
#   scripts/run_phase1.sh [threads]
set -euo pipefail

cd "$(dirname "$0")/.."
threads="${1:-8}"

cmake -S . -B build -DCMAKE_BUILD_TYPE=Release > /dev/null
cmake --build build -j "${threads}"
(cd build && ctest --output-on-failure)

mkdir -p output
./build/triveni-proton Main:numberOfEvents=20000 Main:threads="${threads}" Main:output=output/proton
for shape in hotspot Y Delta; do
  ./build/triveni-pp cards/pp200.cmnd Proton:shape="${shape}" Main:numberOfEvents=50000 \
    Main:threads="${threads}" Main:output="output/pp200_${shape}"
done

root -l -b -q 'macros/plotProton.C("output/proton", "output/figures")'
root -l -b -q 'macros/plotPP.C("output/pp200", "output/figures", "#sqrt{s} = 200 GeV")'
echo "Figures in output/figures/"
