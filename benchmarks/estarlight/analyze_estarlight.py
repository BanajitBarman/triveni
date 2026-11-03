# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Banajit Barman
"""Turns eSTARlight coherent J/psi samples into the same tables triveni-exclusive writes.

    python3 analyze_estarlight.py <runs dir> <system: eAu | ep> <out stem> [--wmin 40] [--wmax 1e9]

Expects run directories <system>_q2_<lo>_<hi>/ with slight.out and run.log (as made by
run_estarlight.sh). Every eSTARlight event is accepted, so each one carries the weight
sigma_window / N_events; a bin's cross section is the sum of the weights of its events.

Writes <stem>_q2.txt (q2lo q2hi sigma err, nb) and <stem>_t.txt (bin absT dsdt err, nb/GeV^2)
with the Q^2 bins of triveni-exclusive. W^2 = 2 m_p E + m_p^2 - Q^2 with E the photon energy in
the target rest frame (the GAMMA line), and |t| from the 't:' line.
"""

import argparse
import math
import re
from pathlib import Path

PROTON_MASS = 0.938272046
Q2_EDGES = [0.0, 1e-3, 1e-2, 1e-1, 1.0, 3.0, 10.0, 30.0]


def total_cross_section_nb(log_path):
    """'Total cross section: 1.995 microb.' -> nb (before the decay branching ratio)."""
    text = Path(log_path).read_text()
    m = re.search(r"Total cross section:\s*([0-9.eE+-]+)\s*(\w+)\.", text)
    if not m:
        raise RuntimeError(f"no total cross section in {log_path}")
    value, unit = float(m.group(1)), m.group(2)
    scale = {"millib": 1e6, "microb": 1e3, "nanob": 1.0, "picob": 1e-3, "femtob": 1e-6}
    return value * scale[unit]


def read_events(path):
    """Yields (E_gamma, Q2, |t|) per event."""
    egamma = q2 = None
    with open(path) as f:
        for line in f:
            if line.startswith("GAMMA:"):
                parts = line.split()
                egamma, q2 = float(parts[1]), float(parts[2])
            elif line.startswith("t:"):
                yield egamma, q2, abs(float(line.split()[1]))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("runs")
    ap.add_argument("system")
    ap.add_argument("stem")
    ap.add_argument("--wmin", type=float, default=40.0)
    ap.add_argument("--wmax", type=float, default=1e9)
    ap.add_argument("--tmax", type=float, default=None, help="upper |t| of the histogram (GeV^2)")
    ap.add_argument("--tbins", type=int, default=60)
    args = ap.parse_args()

    tmax = args.tmax if args.tmax else (0.12 if args.system == "eAu" else 1.5)
    dt = tmax / args.tbins
    nbin = len(Q2_EDGES) - 1
    sigma = [0.0] * nbin
    sumw2 = [0.0] * nbin
    hist = [[0.0] * args.tbins for _ in range(nbin)]
    hist_w2 = [[0.0] * args.tbins for _ in range(nbin)]

    runs = sorted(Path(args.runs).glob(f"{args.system}_q2_*"))
    if not runs:
        raise SystemExit(f"no runs {args.system}_q2_* in {args.runs}")
    for run in runs:
        window_sigma = total_cross_section_nb(run / "run.log")
        events = list(read_events(run / "slight.out"))
        weight = window_sigma / len(events)
        for egamma, q2, abst in events:
            w2 = 2.0 * PROTON_MASS * egamma + PROTON_MASS ** 2 - q2
            W = math.sqrt(max(w2, 0.0))
            if not (args.wmin <= W <= args.wmax):
                continue
            ibin = next((i for i in range(nbin) if Q2_EDGES[i] <= q2 < Q2_EDGES[i + 1]), None)
            if ibin is None:
                continue
            sigma[ibin] += weight
            sumw2[ibin] += weight * weight
            it = int(abst / dt)
            if it < args.tbins:
                hist[ibin][it] += weight / dt
                hist_w2[ibin][it] += (weight / dt) ** 2
        print(f"{run.name}: sigma_window = {window_sigma:.4g} nb, {len(events)} events")

    with open(f"{args.stem}_q2.txt", "w") as out:
        out.write("q2lo/D:q2hi/D:sigma/D:err/D\n")
        for i in range(nbin):
            out.write(f"{Q2_EDGES[i]} {Q2_EDGES[i + 1]} {sigma[i]:.7g} {math.sqrt(sumw2[i]):.4g}\n")
    with open(f"{args.stem}_t.txt", "w") as out:
        out.write("bin/I:absT/D:dsdt/D:err/D\n")
        for i in range(nbin):
            for it in range(args.tbins):
                out.write(f"{i} {(it + 0.5) * dt:.6g} {hist[i][it]:.7g} {math.sqrt(hist_w2[i][it]):.4g}\n")
    print("Q2 bin   sigma(W in [%g, %g]) [nb]" % (args.wmin, args.wmax))
    for i in range(nbin):
        print(f"  [{Q2_EDGES[i]:g}, {Q2_EDGES[i + 1]:g}]: {sigma[i]:.4g} +- {math.sqrt(sumw2[i]):.2g}")


if __name__ == "__main__":
    main()
