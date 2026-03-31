# Tara Architecture

## Module Map

```
tara (lib)
├── star            Core Star type, SpectralClass enum (OBAFGKM + WLTY), StarBuilder
├── constants       Physical and astronomical constants (IAU 2015, CODATA 2018)
├── classification  HR diagram positioning, spectral classification, luminosity classes
├── evolution       Mass-luminosity relation, MS lifetime, evolutionary phases, remnants
├── sse             Analytic stellar evolution (Tout et al. 1996, Hurley et al. 2000)
│   ├── ZAMS        Zero-age main sequence L, R, T from mass and metallicity
│   ├── TMS         Terminal main sequence anchors
│   ├── MS interp   L(τ), R(τ) with hook corrections (α, β, η, γ, δ)
│   ├── HG/RGB      Hertzsprung gap and red giant branch evolution
│   └── evolve()    Unified evolution: mass + Z + age → StellarState
├── spectral        Planck radiance, Doppler shifts, line profiles (Gaussian/Lorentzian/Voigt)
├── nucleosynthesis pp-chain, CNO cycle, triple-alpha rates and branch fractions
├── atmosphere      T_eff, log(g), grey atmosphere, opacity, limb darkening
├── luminosity      Stefan-Boltzmann, magnitudes, distance modulus, bolometric correction
├── bridge          Cross-crate conversions
│   ├── falak       Mass → μ, luminosity class → abs magnitude
│   ├── prakash     Temperature → peak wavelength, gravity → limb darkening
│   └── tanmatra    Composition → lifetime, core temp → CNO fraction
├── integration/
│   └── soorat      Visualization structs and SSE evolution tracks (feature: soorat-compat)
├── logging         Tracing init with try_init() and try_init_dev() (feature: logging)
└── error           TaraError, Result alias
```

## Data Flow

1. Consumers construct `Star` with physical properties or use `sse::evolve()` for age-dependent evolution
2. SSE module provides metallicity-dependent stellar evolution from ZAMS through remnant formation
3. Bridge functions convert between tara values and other AGNOS crate primitives
4. Soorat integration structs package stellar data for visualization rendering

## Consumers

- **AGNOS science stack** — primary consumer, uses tara as the stellar astrophysics vocabulary
- **kiran** — stellar simulation
- **joshua** — science pipeline
- **Simulation apps** — downstream applications consuming tara types

## Design Decisions

- Flat library crate — no nested sub-crates
- Feature-gated optional dependencies — consumers pull only what they need
- All public types are `Serialize + Deserialize` for data interchange
- `#[non_exhaustive]` on all public enums and structs for forward compatibility
- Zero `unwrap`/`panic` in library code — all fallible operations return `Result`
- SSE coefficients sourced from Hurley SSE Fortran source (COSMIC-PopSynth/COSMIC)
- `#[must_use]` on all pure functions
- `#![forbid(unsafe_code)]`
