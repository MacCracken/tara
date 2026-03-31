# Tara

**Tara** (Sanskrit: तारा — star) — stellar astrophysics engine for the AGNOS science stack.

## What it does

Core library for star classification, stellar evolution, nucleosynthesis, spectral analysis, and atmosphere modeling. Provides cross-crate bridges to other AGNOS crates (falak, prakash, tanmatra) and visualization data structures for soorat rendering.

## Modules

| Module | Purpose |
|--------|---------|
| `star` | Core `Star` type and Morgan-Keenan spectral classes (OBAFGKM + WLTY) |
| `constants` | Physical and astronomical constants (IAU 2015, CODATA 2018) |
| `classification` | Hertzsprung-Russell diagram positioning and taxonomy |
| `evolution` | Main sequence, red giant, white dwarf, supernova tracks |
| `sse` | Analytic stellar evolution fitting formulae (Tout et al. 1996, Hurley et al. 2000) |
| `spectral` | Absorption, emission, Doppler shift, line broadening |
| `nucleosynthesis` | Fusion chains, element production, energy release |
| `atmosphere` | Opacity, radiative transfer, limb darkening |
| `luminosity` | Magnitude and distance calculations |
| `bridge` | Cross-crate primitive-value conversions (falak, prakash, tanmatra) |
| `integration::soorat` | Visualization structs and SSE evolution tracks for soorat rendering (feature-gated) |

## Features

| Feature | What it enables |
|---------|----------------|
| `soorat-compat` | Soorat visualization integration |
| `logging` | Tracing-based logging via `tracing-subscriber` |
| `optics` | Prakash optics bridge |
| `thermal` | Ushma thermal bridge |
| `nuclear` | Kimiya nuclear physics bridge |

```toml
[dependencies]
tara = { version = "1.0.0", features = ["soorat-compat"] }
```

## Consumers

AGNOS science stack, kiran, joshua, simulation apps.

## Development

See [CLAUDE.md](CLAUDE.md) for development process, conventions, and work loop.
See [CONTRIBUTING.md](CONTRIBUTING.md) for contribution guidelines.

```sh
make check    # fmt + clippy + test + audit
make bench    # run benchmarks
make doc      # build documentation
```

## License

GPL-3.0-only
