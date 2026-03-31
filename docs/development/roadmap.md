# Tara Roadmap

## Completed

- Initial scaffold

## Backlog

- Star classification (HR diagram)
- Stellar evolution tracks
- Spectral analysis
- Nucleosynthesis chains
- Atmosphere models
- Luminosity/magnitude calculations

## Cross-Crate Bridges

- [ ] `bridge.rs` module — primitive-value conversions for cross-crate stellar astrophysics
- [ ] **falak bridge**: stellar mass (kg) → gravitational parameter; luminosity class → absolute magnitude
- [ ] **prakash bridge**: effective temperature (K) → blackbody spectral energy distribution; surface gravity → limb darkening coefficients
- [ ] **tanmatra bridge**: stellar composition (H/He fraction) → nucleosynthesis yield; core temperature → fusion rate

## Soorat Integration

- [ ] `integration/soorat.rs` module — feature-gated `soorat-compat`
- [ ] **HR diagram data**: temperature vs luminosity points with spectral class for scatter plot rendering
- [ ] **Stellar evolution track**: time-series of (T, L, R) for animated line rendering
- [ ] **Spectral line profile**: wavelength vs intensity for spectral plot rendering
- [ ] **Star field**: positions, magnitudes, color indices for instanced point/billboard rendering
