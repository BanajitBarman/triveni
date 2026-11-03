#!/bin/bash
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Banajit Barman
# Reference eSTARlight samples for the TRIVENI benchmark: coherent J/psi (PROD_MODE 12),
# 18 GeV electrons on 110 GeV/u Au and on 275 GeV protons, four Q^2 windows, 50k events each.
#
#   ESTARLIGHT=~/estarlight/build/e_starlight benchmarks/estarlight/run_estarlight.sh <runs dir>
set -euo pipefail
here="$(cd "$(dirname "$0")" && pwd)"
runs="${1:?usage: run_estarlight.sh <runs dir>}"
exe="${ESTARLIGHT:?set ESTARLIGHT to the e_starlight executable}"
mkdir -p "$runs" && cd "$runs"
seed=100
for sys in "eAu 79 197 117.237" "ep 1 1 293.0939"; do
  for win in "0 1" "1 3" "3 10" "10 30"; do
    set -- $sys $win
    dir="$1_q2_$5_$6"; seed=$((seed + 1))
    "$here/make_card.sh" "$dir" "$2" "$3" "$4" "$5" "$6" 50000 "$seed"
    (cd "$dir" && "$exe" > run.log 2>&1)
    echo "$dir: $(grep 'Total cross section' "$dir/run.log")"
  done
done
