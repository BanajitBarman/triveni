<p align="center">
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="docs/logo/triveni-logo-dark.svg">
    <img src="docs/logo/triveni-logo.svg" alt="TRIVENI" width="480">
  </picture>
</p>

# TRIVENI

**A Y-junction Monte Carlo for ep, eA, pp and AA collisions.**

TRIVENI models the proton as three valence quarks joined by a **Y-shaped
gluonic string junction**, and tests that picture with one parameter set
across ep, eA, pp and AA collisions. Established physics (hard scattering,
showers, hydrodynamics, hadronic afterburners) comes from existing codes
(PYTHIA 8, MUSIC, iSS, SMASH); TRIVENI adds the junction.

The name comes from the Sanskrit *tri* + *veṇī*, "three streams". The
Triveni Sangam is the confluence of three rivers, much as three flux tubes
meet at the baryon junction.

**Status:** version 0.1.0, Phase 1 (proton geometry and pp initial state).

## Model

Quark positions are sampled in 3D, each coordinate Gaussian with width
`Proton:quarkSpread`. The flux tubes minimise their total length, so the
junction sits at the **Fermat–Steiner point** of the quark triangle, where the
three legs meet at 120°. If one angle of the triangle is ≥ 120°, the junction
sits on that quark and the Y has two legs ([Steiner.cxx](src/Steiner.cxx)).

The transverse thickness, normalised to 1, is

```
T_p(x) = (1 - f_s)/3 · Σ_quarks G(x - q_i; w_q)  +  f_s · Σ_tubes (L_j / L) · S_j(x; w_s)
```

- A flux tube is its segment convolved with a 3D Gaussian of width `w_s`. Its
  projection along the beam is therefore exactly the projected segment
  convolved with a 2D Gaussian, which has a closed form in erf
  ([Proton.cxx](src/Proton.cxx)).
- Tube weights are proportional to the 3D length, i.e. to the string energy σL.
- The proton is recentred on its centre of mass.

| Shape | Flux tubes |
|---|---|
| `hotspot` | none (the three-hot-spot proton of the literature) |
| `Y` | quarks → junction |
| `Delta` | the three sides of the quark triangle |
| `YDelta` | Δ when the mean quark separation is below `Proton:dYDelta`, else Y |

**pp collisions** ([PPCollision.cxx](src/PPCollision.cxx)):
- Collision probability: `P(b) = 1 − exp(−σ_gg T_pp(b))`, where `T_pp` is the
  overlap of the two sampled protons, event by event.
- `σ_gg` is fixed so that `∫d²b ⟨P⟩ = PP:sigmaInel`.
- Entropy: `s = ((T₁ᵖ + T₂ᵖ)/2)^{1/p}` (TRENTo), with optional Gamma(k, 1/k)
  weight fluctuations per hot spot and tube.
- Eccentricities εₙ are computed about the entropy centroid.

## Build and run

Requires CMake ≥ 3.20 and a C++17 compiler. ROOT is only used for the plotting
macros.

```sh
scripts/run_phase1.sh 8        # build, test, generate, plot (~3 min on 8 cores)
```

Or step by step:

```sh
cmake -S . -B build && cmake --build build -j8 && (cd build && ctest)
./build/triveni-proton Main:numberOfEvents=20000 Main:threads=8 Main:output=output/proton
./build/triveni-pp cards/pp200.cmnd Proton:shape=Y Main:numberOfEvents=50000 Main:threads=8 Main:output=output/pp200_Y
root -l -b -q 'macros/plotPP.C("output/pp200", "output/figures", "#sqrt{s} = 200 GeV")'
```

**Settings.**
- Settings are PYTHIA-style `Key = value`, given as card files
  ([cards/](cards/)) or on the command line, applied in order.
- `./build/triveni-pp --help` lists every key with its default.
- Unknown keys are errors.

**Reproducibility.**
- Every block of events has its own random stream, so output is identical for
  any `Main:threads`.
- The distributions are implemented in [Random.cxx](src/Random.cxx) rather
  than taken from `<random>`, so the same seed gives the same events with
  libc++ and libstdc++.

**Outputs.** Whitespace-separated text whose first line is a ROOT
`TTree::ReadFile` branch descriptor. They also load with
`numpy.loadtxt(..., skiprows=1)`.

## Validation (ctest)

| Test | Checks |
|---|---|
| `test_random` | moments of uniform, normal and Gamma(k) distributions; stream independence |
| `test_steiner` | 120° legs; vertex case at ≥ 120°; collinear and coincident quarks; agreement with an independent Weiszfeld minimiser over 20 000 random 3D triangles |
| `test_proton` | tube profile against direct quadrature (10⁻⁶); normalisation and centring of all shapes; hot-spot ⟨r²⟩ = ⅔ r_q² + w_q² |
| `test_pp` | reduced-thickness limits; two-Gaussian ε₂ against the analytic value; σ_gg calibration reproduced by the independent accepted fraction |

## First results (placeholder parameters, no fluctuations)

- **Proton geometry.**
  - ⟨L_Y⟩ = 1.34 fm and ⟨L_Δ⟩ = 2.43 fm.
  - L_Y/L_Δ lies in [1/2, 1/√3], as it must.
  - The junction sits on a quark in 22% of configurations.
  - The flux tubes fill the proton centre: ⟨T_p⟩(0) is 2.1 fm⁻² for Y against
    1.4 fm⁻² for hot spots.
- **pp at 200 GeV** (σ_inel = 42 mb):
  - At the same σ_inel, the Y proton reduces ε₃{2} by about 25–35% relative to
    hot spots.
  - ε₂{2} drops by about 5–10%.
- **Energy dependence.** At 5.02 TeV (67.6 mb) the default proton is too small:
  the σ_gg calibration becomes black-disk-like and ill-determined. `triveni-pp`
  warns about this. The proton size must grow with √s, which is the next
  calibration task.

## Roadmap

1. **Phase 1 (current):**
   - Energy-dependent proton size.
   - Calibration of `quarkSpread`, `quarkWidth`, `stringWidth` and
     `stringFraction` against HERA incoherent J/ψ and σ_inel(√s).
   - Weight fluctuations and the high-multiplicity tail.
   - Nuclei (Glauber with Y-shaped nucleons) and eccentricities in pA and AA.
2. **Phase 2:** junction baryon transport (stopping slope α_J and baryon
   fraction p_J) in a PYTHIA 8 fork, tested on ep and pp baryon observables.
3. **Phase 3:** net-baryon initial conditions n_B(x) → MUSIC → iSS → SMASH;
   net-proton dN/dy versus √s; the isobar B/Q ratio.
4. **Phase 4:** eA and γA baryon stopping; predictions for the EIC.

## References

- G. C. Rossi, G. Veneziano, Nucl. Phys. B 123 (1977) 507: the baryon junction.
- D. Kharzeev, Phys. Lett. B 378 (1996) 238: baryon number carried by the junction.
- T. Sjöstrand, P. Z. Skands, Nucl. Phys. B 659 (2003) 243: junction strings in PYTHIA.
- H. Mäntysaari, B. Schenke, Phys. Rev. Lett. 117 (2016) 052301: hot-spot proton from HERA J/ψ.
- J. S. Moreland, J. E. Bernhard, S. A. Bass, Phys. Rev. C 92 (2015) 011901: TRENTo.

## Citing

If you use TRIVENI, please cite it using [CITATION.cff](CITATION.cff); GitHub
shows it as "Cite this repository".

## License

Copyright (C) 2026 Banajit Barman.

TRIVENI is free software: you can redistribute it and/or modify it under the
terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version. See [LICENSE](LICENSE).

GPL-3.0-or-later keeps TRIVENI compatible with the codes it builds on:
PYTHIA 8 (GPL-2.0-or-later, forked in Phase 2), and SMASH and JETSCAPE
(GPL-3.0).

## Logo

The logo is a proton: three valence quarks carrying the colour charges red,
green and blue, joined by Y-shaped flux tubes that meet at a colour-neutral
junction. The files are in [docs/logo/](docs/logo/), with light and dark versions:
SVG, PNG and a PDF for LaTeX. `docs/logo/make_logo.py` regenerates them.
