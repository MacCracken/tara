# Tara Roadmap

## Completed

- Initial scaffold
- Star classification (HR diagram)
- Stellar evolution tracks
- Spectral analysis
- Nucleosynthesis chains
- Atmosphere models
- Luminosity/magnitude calculations

## Cross-Crate Bridges

- [x] `bridge.rs` module — primitive-value conversions for cross-crate stellar astrophysics
- [x] **falak bridge**: stellar mass (kg) → gravitational parameter; luminosity class → absolute magnitude
- [x] **prakash bridge**: effective temperature (K) → blackbody spectral energy distribution; surface gravity → limb darkening coefficients
- [x] **tanmatra bridge**: stellar composition (H/He fraction) → nucleosynthesis yield; core temperature → fusion rate

## Soorat Integration

- [x] `integration/soorat.rs` module — feature-gated `soorat-compat`
- [x] **HR diagram data**: temperature vs luminosity points with spectral class for scatter plot rendering
- [x] **Stellar evolution track**: time-series of (T, L, R) for animated line rendering
- [x] **Spectral line profile**: wavelength vs intensity for spectral plot rendering
- [x] **Star field**: positions, magnitudes, color indices for instanced point/billboard rendering

## Backlog

- Soorat integration helper functions (populate structs from Star instances)
- Usage examples in `examples/`
- Hurley/SSE fitting formulae for detailed evolution
- Extended spectral classes (L, T, Y for brown dwarfs; W for Wolf-Rayet)
- Multi-band bolometric corrections (Flower 1996 / Pecaut & Mamajek 2013)
