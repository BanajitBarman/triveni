#!/bin/bash
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Banajit Barman
# make_card.sh <dir> <Z> <A> <target_gamma> <q2min> <q2max> <nevents> <seed>
mkdir -p "$1"
cat > "$1/slight.in" <<CARD
baseFileName = slight
TARGET_BEAM_Z = $2
TARGET_BEAM_A = $3
ELECTRON_BEAM_GAMMA = 35225.05    # 18 GeV electrons
TARGET_BEAM_GAMMA = $4
W_MAX = -1
W_MIN = -1
W_N_BINS = 50
W_GP_MAX = -1
W_GP_MIN = -1
EGA_N_BINS = 400
CUT_PT = 0
PT_MIN = 1.0
PT_MAX = 3.0
PROD_MODE = 12
N_EVENTS = $7
PROD_PID = 443013
RND_SEED = $8
BREAKUP_MODE = 5
OUTPUT_FORMAT = 0
MIN_GAMMA_Q2 = $5
MAX_GAMMA_Q2 = $6
QUANTUM_GLAUBER = 1
SELECT_IMPULSE_VM = 0
CARD
