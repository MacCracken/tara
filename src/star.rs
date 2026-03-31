//! Core star representation and spectral classes.

use std::fmt;

use serde::{Deserialize, Serialize};

use crate::error::{Result, TaraError};

/// A star with fundamental physical properties.
#[derive(Debug, Clone, Serialize, Deserialize)]
#[non_exhaustive]
pub struct Star {
    /// Mass in solar masses (M_sun).
    pub mass_solar: f64,
    /// Radius in solar radii (R_sun).
    pub radius_solar: f64,
    /// Effective surface temperature in Kelvin.
    pub temperature_k: f64,
    /// Luminosity in solar luminosities (L_sun).
    pub luminosity_solar: f64,
    /// Age in years.
    pub age_years: f64,
    /// Spectral classification.
    pub spectral_class: SpectralClass,
}

impl Star {
    /// Create a new `Star` with the given properties.
    ///
    /// # Errors
    ///
    /// Returns [`TaraError::InvalidParameter`] if any value is non-positive.
    pub fn new(
        mass_solar: f64,
        radius_solar: f64,
        temperature_k: f64,
        luminosity_solar: f64,
        age_years: f64,
        spectral_class: SpectralClass,
    ) -> Result<Self> {
        if mass_solar <= 0.0 {
            return Err(TaraError::InvalidParameter(
                "mass_solar must be positive".into(),
            ));
        }
        if radius_solar <= 0.0 {
            return Err(TaraError::InvalidParameter(
                "radius_solar must be positive".into(),
            ));
        }
        if temperature_k <= 0.0 {
            return Err(TaraError::InvalidParameter(
                "temperature_k must be positive".into(),
            ));
        }
        if luminosity_solar <= 0.0 {
            return Err(TaraError::InvalidParameter(
                "luminosity_solar must be positive".into(),
            ));
        }
        if age_years < 0.0 {
            return Err(TaraError::InvalidParameter(
                "age_years must be non-negative".into(),
            ));
        }

        Ok(Self {
            mass_solar,
            radius_solar,
            temperature_k,
            luminosity_solar,
            age_years,
            spectral_class,
        })
    }
}

/// Morgan-Keenan spectral classification of stars.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
#[non_exhaustive]
pub enum SpectralClass {
    /// O-type: hot, blue, >30,000 K.
    O,
    /// B-type: blue-white, 10,000-30,000 K.
    B,
    /// A-type: white, 7,500-10,000 K.
    A,
    /// F-type: yellow-white, 6,000-7,500 K.
    F,
    /// G-type: yellow (Sun-like), 5,200-6,000 K.
    G,
    /// K-type: orange, 3,700-5,200 K.
    K,
    /// M-type: red, <3,700 K.
    M,
}

impl fmt::Display for SpectralClass {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::O => write!(f, "O"),
            Self::B => write!(f, "B"),
            Self::A => write!(f, "A"),
            Self::F => write!(f, "F"),
            Self::G => write!(f, "G"),
            Self::K => write!(f, "K"),
            Self::M => write!(f, "M"),
        }
    }
}
