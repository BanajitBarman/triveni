# Benchmark: coherent J/ψ at the EIC, TRIVENI vs eSTARlight

TRIVENI's coherent vector-meson module (`include/triveni/Exclusive.h`) is compared with
[eSTARlight](https://github.com/eic/estarlight) (Lomnitz & Klein, PRC 99 (2019) 015203)
for e + p and e + Au → e + J/ψ + target, at 18 GeV electrons on 275 GeV protons and on
110 GeV/u gold. The comparison is restricted to where the physics overlaps: a Woods–Saxon
nucleus with no sub-nucleon structure, and W > 40 GeV (x ≲ 0.01, where the dipole picture
applies).

## What is compared

| | eSTARlight | TRIVENI |
|---|---|---|
| γp → J/ψ p | HERA fit: 4.06 nb · W^0.65 · threshold factor | dipole model **fitted to the same curve** (Ω₀ = 1.210, λ = 0.194, B_p = 4 GeV⁻²) |
| Q² dependence | (M²/(M²+Q²))ⁿ, n = 2.45 + 0.00084(Q²+M²) | from the photon and boosted-Gaussian J/ψ wavefunctions (KMW 2006) |
| Nucleus | quantum Glauber with σ_VN, Woods–Saxon R = 6.38 fm, a = 0.53 fm | dipole Glauber ⟨N_A⟩ = 1 − exp(−σ_p(r) T_A(b)/2), same Woods–Saxon |
| t-shape | Klein–Nystrand form factor (hard sphere ⊗ Yukawa) | Fourier transform of ⟨N_A⟩ |
| Photon flux | Weizsäcker–Williams, target rest frame | **the same formula**, σ_T + ε σ_L |

Because the e + p input and the photon flux are the same, differences in e + Au measure the
nuclear physics alone.

## Results (W > 40 GeV; σ in nb)

e + p:

| Q² (GeV²) | eSTARlight | TRIVENI | TRIVENI / eSTARlight |
|---|---|---|---|
| 0.001–0.01 | 0.653 ± 0.011 | 0.638 | 0.98 |
| 0.01–0.1 | 0.645 ± 0.010 | 0.632 | 0.98 |
| 0.1–1 | 0.574 ± 0.010 | 0.579 | 1.01 |
| 1–3 | 0.206 ± 0.001 | 0.198 | 0.96 |
| 3–10 | 0.116 ± 0.001 | 0.106 | 0.92 |
| 10–30 | 0.0273 ± 0.0002 | 0.0218 | 0.80 |

e + Au (W from 40 to 89 GeV, the kinematic limit):

| Q² (GeV²) | eSTARlight | eSTARlight, test patch | TRIVENI | TRIVENI / eSTARlight | TRIVENI / test |
|---|---|---|---|---|---|
| 0.001–0.01 | 157.0 ± 2.5 | 157.2 | 120.3 | 0.77 | 0.77 |
| 0.01–0.1 | 156.0 ± 2.5 | 157.0 | 119.3 | 0.76 | 0.76 |
| 0.1–1 | 132.3 ± 2.3 | 142.1 | 109.8 | 0.83 | 0.77 |
| 1–3 | 35.5 ± 0.2 | 51.1 | 38.3 | 1.08 | 0.75 |
| 3–10 | 11.5 ± 0.1 | 30.1 | 21.2 | 1.84 | 0.70 |
| 10–30 | 0.846 ± 0.005 | 7.43 | 4.62 | 5.46 | 0.62 |

Figures: `bench_sigma` (γp fit, dσ/dQ² for e + p and e + Au, Au/p ratio) and `bench_t`
(e + Au dσ/d|t| in three Q² bins).

## Findings

1. **e + p agrees to 2–8% for Q² < 10 GeV²** and to 20% at 10–30 GeV². TRIVENI's Q²
   dependence comes only from the light-cone wavefunctions, eSTARlight's from a HERA fit.
2. **The Woods–Saxon geometry agrees:** the first coherent diffraction minimum is at
   |t| ≈ 0.017 GeV² in both codes.
3. **In photoproduction, TRIVENI's e + Au is 17–24% below eSTARlight.** This is the expected
   difference between dipole-model saturation and VMD Glauber shadowing.
4. **At high Q², eSTARlight's e + Au falls far below** (by 5.5× at 10–30 GeV²), while its
   e + p does not. The σ(e+Au)/σ(e+p) ratio drops from about 390 to 50 between
   photoproduction and Q² ≈ 15 GeV². TRIVENI's ratio rises slightly (299 → 318), as reduced
   shadowing at higher Q² implies.

   In eSTARlight's `photonNucleusCrossSection::getcsgA`, nuclei use
   `cs = getcsgA_Q2_dep(Q2) * sigma_N(Wgp)`: the σ-level Q² factor multiplies σ_VN, and
   dσ/dt(t=0) ∝ σ_VA². In the dilute limit the factor is therefore applied twice for nuclei
   and once for protons, where `csgA = getcsgA_Q2_dep(Q2) * sigmagp(Wgp)`.

   The test patch (branch `triveni-q2-test` in a local eSTARlight clone) uses
   `std::sqrt(getcsgA_Q2_dep(Q2))` instead. With it, eSTARlight's Au/p ratio becomes flat
   to rising (384 → 437), and TRIVENI / eSTARlight is a steady 0.62–0.77. This looks like
   an inconsistency in eSTARlight and should be confirmed with its authors before relying
   on eSTARlight's e + A cross sections at Q² ≳ 1 GeV².

## Numerical precision

TRIVENI's (W, Q²) integration grid (33 × 41 Simpson points in ln W and ln Q²) agrees with
an independent 65 × 81 grid to 0.5% (e + Au, 1 < Q² < 3 GeV²: 38.31 vs 38.49 nb). The
eSTARlight numbers carry the statistical errors of 50,000 events per Q² window.

## Known limitations of the TRIVENI side

- No real-part or skewness corrections; their effect on e + p is absorbed by the fit.
- No DGLAP evolution in the dipole (Ω(x) = Ω₀ x^−λ), so scaling violations at high Q² are
  missing. This is a likely reason for the 20% e + p difference at 10–30 GeV².
- No t_min (longitudinal momentum transfer), which matters only at the largest Q².
- Coherent production only. Incoherent production with Y-junction protons is the next step.

## Reproducing

```sh
# eSTARlight (outside this repository), at commit 939b11a
git clone https://github.com/eic/estarlight.git ~/estarlight
cmake -S ~/estarlight -B ~/estarlight/build -DCMAKE_BUILD_TYPE=Release && cmake --build ~/estarlight/build -j8
ESTARLIGHT=~/estarlight/build/e_starlight benchmarks/estarlight/run_estarlight.sh ~/estarlight/runs

# eSTARlight tables with the same bins and W cut as TRIVENI
python3 benchmarks/estarlight/analyze_estarlight.py ~/estarlight/runs eAu output/bench/estarlight_eAu
python3 benchmarks/estarlight/analyze_estarlight.py ~/estarlight/runs ep  output/bench/estarlight_ep
python3 benchmarks/estarlight/analyze_estarlight.py ~/estarlight/runs ep  output/bench/estarlight_ep_w89 --wmax 88.9936

# TRIVENI
./build/triveni-exclusive Exclusive:target=Au Main:threads=8 Main:output=output/bench/triveni_eAu
./build/triveni-exclusive Exclusive:target=p Exclusive:targetGamma=293.0939 Main:threads=8 Main:output=output/bench/triveni_ep
./build/triveni-exclusive Exclusive:target=p Exclusive:targetGamma=293.0939 Exclusive:wMax=88.9936 Main:threads=8 Main:output=output/bench/triveni_ep_w89

root -l -b -q 'benchmarks/estarlight/compare.C("output/bench", "output/figures")'
```

For the test-patch curve, build the `triveni-q2-test` branch of the eSTARlight clone,
rerun the e + Au windows, and analyse them into `output/bench/estarlight_eAu_q2test`.
