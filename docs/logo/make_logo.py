# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Banajit Barman
"""Draws the TRIVENI logo as plain SVG geometry (no fonts), plus PNG/PDF exports.

The mark is a proton: three valence quarks in the QCD colour charges red, green
and blue, joined by Y-shaped flux tubes that meet at a colour-neutral (white)
junction. The wordmark is drawn letter by letter, so the files look the same
everywhere.

    python3 docs/logo/make_logo.py      # rsvg-convert (librsvg) is used for PNG/PDF if present
"""

import math
import shutil
import subprocess
from pathlib import Path

OUT = Path(__file__).parent
NAVY, WHITE = "#003366", "#FFFFFF"
RED, GREEN, BLUE = "#E0453A", "#2BA35A", "#2F6FDB"
# Y directions in SVG coordinates (y points down): upper left, upper right, down.
DIRS = ((-math.sqrt(3) / 2, -0.5), (math.sqrt(3) / 2, -0.5), (0.0, 1.0))


def f(v):
    return f"{v:.2f}".rstrip("0").rstrip(".")


def mark(ink):
    """The mark, centred on (0, 0); it reaches radius 87 (ring radius 84 + half its width)."""
    leg = 50.0
    parts = [f'<circle cx="0" cy="0" r="84" fill="none" stroke="{ink}" stroke-width="6"/>']
    for dx, dy in DIRS:
        parts.append(f'<line x1="0" y1="0" x2="{f(leg * dx)}" y2="{f(leg * dy)}" stroke="{ink}" '
                     f'stroke-width="8" stroke-linecap="round"/>')
    for (dx, dy), colour in zip(DIRS, (RED, GREEN, BLUE)):
        parts.append(f'<circle cx="{f(leg * dx)}" cy="{f(leg * dy)}" r="14" fill="{colour}"/>')
    parts.append(f'<circle cx="0" cy="0" r="7" fill="{WHITE}" stroke="{ink}" stroke-width="3"/>')
    return "".join(parts)


# Geometric capitals: cap height 60, stroke 9. (letter, paths, advance width)
LETTERS = [
    ("T", ["M0 0H44V9H26.5V60H17.5V9H0Z"], 44),
    ("R", ["M0 0H9V60H0Z", "M9 0H22A18 18 0 0 1 22 36H9V27H22A9 9 0 0 0 22 9H9Z", "M16 33H26.5L41 60H30.5Z"], 41),
    ("I", ["M0 0H9V60H0Z"], 9),
    ("V", ["M0 0H10L23 44L36 0H46L28.5 60H17.5Z"], 46),
    ("E", ["M0 0H34V9H9V25.5H30V34.5H9V51H34V60H0Z"], 34),
    ("N", ["M0 0H9V60H0Z", "M33 0H42V60H33Z", "M0 0H11L42 60H31Z"], 42),
    ("I", ["M0 0H9V60H0Z"], 9),
]
SPACING = 13


def wordmark(x0, y0, ink):
    parts, x = [], 0.0
    for _, paths, width in LETTERS:
        parts += [f'<path transform="translate({f(x)} 0)" d="{d}"/>' for d in paths]
        x += width + SPACING
    return f'<g transform="translate({f(x0)} {f(y0)})" fill="{ink}">{"".join(parts)}</g>', x - SPACING


def svg(width, height, body, title):
    return (f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {f(width)} {f(height)}" '
            f'width="{f(width)}" height="{f(height)}" role="img" aria-label="{title}">'
            f"<title>{title}</title>{body}</svg>\n")


def lockup(ink):
    """Mark and wordmark side by side; 23 units of margin all round."""
    word, word_width = wordmark(232, 80, ink)
    body = f'<g transform="translate(110 110)">{mark(ink)}</g>{word}'
    return svg(232 + word_width + 23, 220, body, "TRIVENI")


def mark_only(ink):
    return svg(200, 200, f'<g transform="translate(100 100)">{mark(ink)}</g>', "TRIVENI")


def social_preview():
    """1280 x 640 card for GitHub's social preview (Settings > General > Social preview)."""
    word, word_width = wordmark(0, 0, WHITE)
    scale = 1.6
    # Visible extent of the lockup: ring edge at x = 23 to the end of the wordmark.
    visible = (232 + word_width - 23) * scale
    x0 = (1280 - visible) / 2 - 23 * scale
    body = (f'<rect width="1280" height="640" fill="{NAVY}"/>'
            f'<g transform="translate({f(x0)} 100) scale({scale})">'
            f'<g transform="translate(110 110)">{mark(WHITE)}</g>'
            f'<g transform="translate(232 80)">{word}</g></g>'
            f'<text x="640" y="512" text-anchor="middle" fill="#C9D6E8" font-size="34" '
            f'font-family="Avenir Next, Helvetica Neue, Helvetica, Arial, sans-serif">'
            f"A Y-junction Monte Carlo for ep, eA, pp and AA collisions</text>")
    return svg(1280, 640, body, "TRIVENI")


def main():
    files = {
        "triveni-logo.svg": lockup(NAVY),
        "triveni-logo-dark.svg": lockup(WHITE),
        "triveni-mark.svg": mark_only(NAVY),
        "triveni-mark-dark.svg": mark_only(WHITE),
        "triveni-social.svg": social_preview(),
    }
    for name, text in files.items():
        (OUT / name).write_text(text)

    rsvg = shutil.which("rsvg-convert")
    if not rsvg:
        print("rsvg-convert not found: wrote SVG files only")
        return
    exports = [
        ("triveni-logo.svg", ["-w", "1200"], "triveni-logo.png"),
        ("triveni-logo-dark.svg", ["-w", "1200"], "triveni-logo-dark.png"),
        ("triveni-logo.svg", ["-f", "pdf"], "triveni-logo.pdf"),
        ("triveni-mark.svg", ["-w", "512"], "triveni-mark-512.png"),
        ("triveni-mark.svg", ["-w", "64"], "favicon-64.png"),
        ("triveni-social.svg", ["-w", "1280"], "triveni-social.png"),
    ]
    for source, options, target in exports:
        subprocess.run([rsvg, *options, str(OUT / source), "-o", str(OUT / target)], check=True)
    (OUT / "triveni-social.svg").unlink()  # its tagline uses a system font; keep only the PNG
    print("wrote", ", ".join(sorted(p.name for p in OUT.iterdir() if p.suffix in {".svg", ".png", ".pdf"})))


if __name__ == "__main__":
    main()
