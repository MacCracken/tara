//! # Tara
//!
//! **Tara** (Sanskrit: तारा — star) — stellar astrophysics engine for star
//! classification, evolution, nucleosynthesis, and spectral analysis.

#![warn(missing_docs)]
#![forbid(unsafe_code)]

/// Physical and astronomical constants (IAU 2015, CODATA 2018).
pub mod constants;

/// Cross-crate bridges — primitive-value conversions from other AGNOS science crates.
pub mod bridge;
/// Error types for the Tara engine.
pub mod error;
/// Integration APIs for downstream consumers (soorat rendering).
pub mod integration;

/// Core star representation and spectral classes.
pub mod star;

/// Hertzsprung-Russell diagram classification and stellar taxonomy.
pub mod classification;

/// Stellar evolution — main sequence, red giant, white dwarf, supernova tracks.
pub mod evolution;

/// Spectral line analysis — absorption, emission, Doppler shift, broadening.
pub mod spectral;

/// Stellar nucleosynthesis — fusion chains, element production, energy release.
pub mod nucleosynthesis;

/// Stellar atmosphere models — opacity, radiative transfer, limb darkening.
pub mod atmosphere;

/// Luminosity, magnitude, and distance calculations.
pub mod luminosity;

/// Analytic stellar evolution fitting formulae (Tout et al. 1996, Hurley et al. 2000).
pub mod sse;

/// Optional logging initialization (requires `logging` feature).
#[cfg(feature = "logging")]
pub mod logging;
