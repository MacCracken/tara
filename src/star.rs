//! Core star representation and spectral classes.

use std::fmt;

use serde::{Deserialize, Serialize};

use crate::classification::LuminosityClass;
use crate::error::{Result, TaraError};

/// A star with fundamental physical properties.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
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
    /// Spectral classification (OBAFGKM).
    pub spectral_class: SpectralClass,
    /// Spectral subclass (0–9). 0 = hottest in class.
    pub spectral_subclass: u8,
    /// Yerkes luminosity class (I–VII).
    pub luminosity_class: LuminosityClass,
    /// Metallicity [Fe/H] in dex (solar = 0.0).
    pub metallicity: f64,
}

impl Star {
    /// Create a new `Star` with the given properties.
    ///
    /// # Errors
    ///
    /// Returns [`TaraError::InvalidParameter`] if any physical value is non-positive.
    #[must_use = "returns a Result containing the new Star"]
    pub fn new(
        mass_solar: f64,
        radius_solar: f64,
        temperature_k: f64,
        luminosity_solar: f64,
        age_years: f64,
        spectral_class: SpectralClass,
    ) -> Result<Self> {
        Self::builder(
            mass_solar,
            radius_solar,
            temperature_k,
            luminosity_solar,
            age_years,
        )
        .spectral_class(spectral_class)
        .build()
    }

    /// Create a builder for more detailed star construction.
    ///
    /// Spectral class and subclass are auto-derived from temperature
    /// unless overridden.
    #[must_use]
    pub fn builder(
        mass_solar: f64,
        radius_solar: f64,
        temperature_k: f64,
        luminosity_solar: f64,
        age_years: f64,
    ) -> StarBuilder {
        StarBuilder {
            mass_solar,
            radius_solar,
            temperature_k,
            luminosity_solar,
            age_years,
            spectral_class: None,
            spectral_subclass: None,
            luminosity_class: LuminosityClass::V,
            metallicity: 0.0,
        }
    }

    /// Create a Sun-like star with IAU 2015 nominal values.
    ///
    /// # Errors
    ///
    /// Cannot fail with standard solar parameters.
    pub fn sun() -> Result<Self> {
        Self::builder(1.0, 1.0, crate::constants::T_SUN, 1.0, 4.6e9)
            .luminosity_class(LuminosityClass::V)
            .build()
    }
}

/// Builder for constructing [`Star`] instances with optional parameters.
pub struct StarBuilder {
    mass_solar: f64,
    radius_solar: f64,
    temperature_k: f64,
    luminosity_solar: f64,
    age_years: f64,
    spectral_class: Option<SpectralClass>,
    spectral_subclass: Option<u8>,
    luminosity_class: LuminosityClass,
    metallicity: f64,
}

impl StarBuilder {
    /// Override the auto-derived spectral class.
    #[must_use]
    pub fn spectral_class(mut self, class: SpectralClass) -> Self {
        self.spectral_class = Some(class);
        self
    }

    /// Override the auto-derived spectral subclass (0–9).
    #[must_use]
    pub fn spectral_subclass(mut self, sub: u8) -> Self {
        self.spectral_subclass = Some(sub.min(9));
        self
    }

    /// Set the luminosity class (default: V, main sequence).
    #[must_use]
    pub fn luminosity_class(mut self, lc: LuminosityClass) -> Self {
        self.luminosity_class = lc;
        self
    }

    /// Set metallicity [Fe/H] in dex (default: 0.0, solar).
    #[must_use]
    pub fn metallicity(mut self, feh: f64) -> Self {
        self.metallicity = feh;
        self
    }

    /// Build the [`Star`], validating all parameters.
    ///
    /// # Errors
    ///
    /// Returns [`TaraError::InvalidParameter`] if any physical value is non-positive.
    pub fn build(self) -> Result<Star> {
        if self.mass_solar <= 0.0 {
            let err = TaraError::InvalidParameter(format!(
                "mass_solar must be positive, got {}",
                self.mass_solar
            ));
            tracing::warn!(mass_solar = self.mass_solar, "{err}");
            return Err(err);
        }
        if self.radius_solar <= 0.0 {
            let err = TaraError::InvalidParameter(format!(
                "radius_solar must be positive, got {}",
                self.radius_solar
            ));
            tracing::warn!(radius_solar = self.radius_solar, "{err}");
            return Err(err);
        }
        if self.temperature_k <= 0.0 {
            let err = TaraError::InvalidParameter(format!(
                "temperature_k must be positive, got {}",
                self.temperature_k
            ));
            tracing::warn!(temperature_k = self.temperature_k, "{err}");
            return Err(err);
        }
        if self.luminosity_solar <= 0.0 {
            let err = TaraError::InvalidParameter(format!(
                "luminosity_solar must be positive, got {}",
                self.luminosity_solar
            ));
            tracing::warn!(luminosity_solar = self.luminosity_solar, "{err}");
            return Err(err);
        }
        if self.age_years < 0.0 {
            let err = TaraError::InvalidParameter(format!(
                "age_years must be non-negative, got {}",
                self.age_years
            ));
            tracing::warn!(age_years = self.age_years, "{err}");
            return Err(err);
        }

        let spectral_class = self.spectral_class.unwrap_or_else(|| {
            crate::classification::spectral_class_from_temperature(self.temperature_k)
        });
        let spectral_subclass = self
            .spectral_subclass
            .unwrap_or_else(|| crate::classification::spectral_subclass(self.temperature_k));

        tracing::debug!(
            mass = self.mass_solar,
            radius = self.radius_solar,
            temp = self.temperature_k,
            class = %spectral_class,
            sub = spectral_subclass,
            "star constructed"
        );

        Ok(Star {
            mass_solar: self.mass_solar,
            radius_solar: self.radius_solar,
            temperature_k: self.temperature_k,
            luminosity_solar: self.luminosity_solar,
            age_years: self.age_years,
            spectral_class,
            spectral_subclass,
            luminosity_class: self.luminosity_class,
            metallicity: self.metallicity,
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn star_serde_roundtrip() {
        let star = Star::new(1.0, 1.0, 5772.0, 1.0, 4.6e9, SpectralClass::G).unwrap();
        let json = serde_json::to_string(&star).unwrap();
        let back: Star = serde_json::from_str(&json).unwrap();
        assert!((back.mass_solar - 1.0).abs() < f64::EPSILON);
        assert!((back.temperature_k - 5772.0).abs() < f64::EPSILON);
        assert_eq!(back.spectral_class, SpectralClass::G);
    }

    #[test]
    fn spectral_class_serde_roundtrip() {
        for class in [
            SpectralClass::O,
            SpectralClass::B,
            SpectralClass::A,
            SpectralClass::F,
            SpectralClass::G,
            SpectralClass::K,
            SpectralClass::M,
        ] {
            let json = serde_json::to_string(&class).unwrap();
            let back: SpectralClass = serde_json::from_str(&json).unwrap();
            assert_eq!(back, class);
        }
    }

    #[test]
    fn star_rejects_negative_mass() {
        assert!(Star::new(-1.0, 1.0, 5772.0, 1.0, 0.0, SpectralClass::G).is_err());
    }

    #[test]
    fn star_rejects_negative_age() {
        assert!(Star::new(1.0, 1.0, 5772.0, 1.0, -1.0, SpectralClass::G).is_err());
    }

    #[test]
    fn sun_convenience() {
        let sun = Star::sun().unwrap();
        assert_eq!(sun.spectral_class, SpectralClass::G);
        assert_eq!(sun.spectral_subclass, 2);
        assert_eq!(sun.luminosity_class, LuminosityClass::V);
        assert!((sun.metallicity).abs() < f64::EPSILON);
    }

    #[test]
    fn builder_auto_classifies() {
        let star = Star::builder(1.0, 1.0, 5772.0, 1.0, 4.6e9).build().unwrap();
        assert_eq!(star.spectral_class, SpectralClass::G);
        assert_eq!(star.spectral_subclass, 2);
    }

    #[test]
    fn builder_override_class() {
        let star = Star::builder(1.0, 1.0, 5772.0, 1.0, 4.6e9)
            .spectral_class(SpectralClass::F)
            .spectral_subclass(5)
            .luminosity_class(LuminosityClass::IV)
            .metallicity(-0.5)
            .build()
            .unwrap();
        assert_eq!(star.spectral_class, SpectralClass::F);
        assert_eq!(star.spectral_subclass, 5);
        assert_eq!(star.luminosity_class, LuminosityClass::IV);
        assert!((star.metallicity - (-0.5)).abs() < f64::EPSILON);
    }
}
