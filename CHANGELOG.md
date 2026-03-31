# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/).

## [1.0.0] — 2026-03-31

### Added

- Core `Star` type with builder pattern, validation, and auto-classification from temperature
- `SpectralClass` enum (OBAFGKM) with Morgan-Keenan temperature boundaries
- `StarBuilder` with optional spectral class/subclass override, luminosity class, and metallicity
- `Star::sun()` convenience constructor with IAU 2015 nominal solar values
- Physical and astronomical constants module (IAU 2015, CODATA 2018)
- `TaraError` enum with `InvalidParameter`, `MathError`, `SpectralError`, `ModelError`, `Io` variants
- Hertzsprung-Russell diagram classification: spectral class from temperature, subclass interpolation, HR region identification, luminosity class from log(g), MK classification string formatting
- Stellar evolution: piecewise mass-luminosity relation (Duric 2004), main-sequence lifetime estimates, remnant type prediction, evolutionary phase determination including PostAgb transition, Chandrasekhar-limited white dwarf initial-final mass relation (Kalirai et al. 2008)
- Spectral analysis: Planck radiance (m and nm), non-relativistic and relativistic Doppler shift, thermal broadening, Gaussian/Lorentzian/pseudo-Voigt (Thompson et al. 1987) line profiles
- Nucleosynthesis: pp-chain, CNO cycle, and triple-alpha energy generation rates with calibrated power-law approximations (Hansen, Kawaler & Trimble 2004; Kippenhahn, Weigert & Weiss 2012), dominant process determination, pp-chain branch fraction model
- Stellar atmosphere: effective temperature (inverse Stefan-Boltzmann), log(g) surface gravity, grey atmosphere temperature profile (Eddington approximation), Kramers' opacity, electron scattering opacity, pressure scale height, linear and quadratic limb darkening
- Luminosity and magnitude: Stefan-Boltzmann luminosity, absolute bolometric magnitude, distance modulus, apparent/absolute magnitude conversion, bolometric correction from temperature
- Cross-crate bridges for falak (gravitational parameter, absolute magnitude from luminosity class), prakash (Wien peak wavelength, limb darkening from log(g)), and tanmatra (composition lifetime scaling, CNO fraction sigmoid)
- Soorat integration module (feature-gated `soorat-compat`): `HrDiagramPoint`, `EvolutionTrack`, `SpectralProfile`, `StarField`, `StarViz` visualization data structures
- Feature gates: `logging`, `optics`, `thermal`, `nuclear`, `soorat-compat`
- Logging module with `try_init()` (production, defaults to `warn`) and `try_init_dev()` (development, compact format with file/line, defaults to `debug`)
- Structured `tracing` instrumentation: `warn` on validation failures and non-physical inputs (star construction, evolution, atmosphere, spectral), `debug` on star construction, `trace` on Planck overflow guard
- Validation error messages include the actual offending value
- `#[non_exhaustive]` on all public enums and soorat structs
- `#[must_use]` on all pure functions
- `#![forbid(unsafe_code)]`
- Serde `Serialize` + `Deserialize` on all public types
- 89 unit tests with serde roundtrip coverage for every public type
- 7 cross-module integration tests (sun end-to-end, massive star, bridge round-trips, spectral pipeline, atmosphere chain, magnitude-distance consistency)
- Criterion benchmarks for star construction, bridge conversions, luminosity, atmosphere, evolution, spectral, and nucleosynthesis functions
