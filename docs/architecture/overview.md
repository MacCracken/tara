# Tara Architecture

## Module Map

```
tara (lib)
├── star            Core Star type, SpectralClass enum
├── classification  HR diagram positioning, stellar taxonomy
├── evolution       Stellar evolution tracks (stub)
├── spectral        Spectral line analysis (stub)
├── nucleosynthesis Fusion chains, element production (stub)
├── atmosphere      Atmosphere models (stub)
├── luminosity      Magnitude/distance calculations (stub)
├── bridge          Cross-crate conversions
│   ├── falak       Mass → μ, luminosity class → abs magnitude
│   ├── prakash     Temperature → peak wavelength, gravity → limb darkening
│   └── tanmatra    Composition → lifetime, core temp → CNO fraction
├── integration/
│   └── soorat      Visualization structs (feature: soorat-compat)
├── logging         Tracing init (feature: logging)
└── error           TaraError, Result alias
```

## Data Flow

1. Consumers construct `Star` with physical properties
2. Bridge functions convert between tara values and other AGNOS crate primitives
3. Soorat integration structs package stellar data for visualization rendering

## Consumers

- **AGNOS science stack** — primary consumer, uses tara as the stellar astrophysics vocabulary
- **kiran** — stellar simulation
- **joshua** — science pipeline
- **Simulation apps** — downstream applications consuming tara types

## Design Decisions

- Flat library crate — no nested sub-crates
- Feature-gated optional dependencies — consumers pull only what they need
- All public types are `Serialize + Deserialize` for data interchange
- `#[non_exhaustive]` on all public enums for forward compatibility
- Zero `unwrap`/`panic` in library code — all fallible operations return `Result`
