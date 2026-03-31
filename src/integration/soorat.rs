//! Soorat integration — visualization data structures for stellar astrophysics.

use serde::{Deserialize, Serialize};

/// HR diagram data point for scatter plot rendering.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[non_exhaustive]
pub struct HrDiagramPoint {
    /// Effective temperature (K).
    pub temperature: f64,
    /// Luminosity (solar luminosities).
    pub luminosity: f64,
    /// Spectral class (e.g. "G2V").
    pub spectral_class: String,
    /// RGB color from temperature.
    pub color: [f32; 3],
}

/// Stellar evolution track for animated line rendering.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[non_exhaustive]
pub struct EvolutionTrack {
    /// Time-series of (temperature_K, luminosity_solar, radius_solar).
    pub points: Vec<[f64; 3]>,
    /// Age at each point (years).
    pub ages: Vec<f64>,
    /// Star name/label.
    pub label: String,
}

/// Spectral line profile for plot rendering.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[non_exhaustive]
pub struct SpectralProfile {
    /// Wavelength samples (nm).
    pub wavelengths: Vec<f64>,
    /// Intensity at each wavelength (normalized 0–1).
    pub intensities: Vec<f64>,
    /// Absorption line positions (nm).
    pub absorption_lines: Vec<f64>,
}

/// Star field for instanced point/billboard rendering.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[non_exhaustive]
pub struct StarField {
    /// Stars with position and properties.
    pub stars: Vec<StarViz>,
}

/// A single star for rendering.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[non_exhaustive]
pub struct StarViz {
    /// Position `[x, y, z]` (arbitrary units).
    pub position: [f32; 3],
    /// Apparent magnitude (for brightness).
    pub magnitude: f32,
    /// B-V color index (for tint).
    pub color_index: f32,
    /// Display radius.
    pub radius: f32,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn hr_point_serde_roundtrip() {
        let p = HrDiagramPoint {
            temperature: 5778.0,
            luminosity: 1.0,
            spectral_class: "G2V".into(),
            color: [1.0, 0.95, 0.85],
        };
        let json = serde_json::to_string(&p).unwrap();
        let back: HrDiagramPoint = serde_json::from_str(&json).unwrap();
        assert_eq!(back, p);
    }

    #[test]
    fn evolution_track_serde_roundtrip() {
        let track = EvolutionTrack {
            points: vec![[5778.0, 1.0, 1.0], [4500.0, 100.0, 50.0]],
            ages: vec![0.0, 1e10],
            label: "Solar".into(),
        };
        let json = serde_json::to_string(&track).unwrap();
        let back: EvolutionTrack = serde_json::from_str(&json).unwrap();
        assert_eq!(back, track);
    }

    #[test]
    fn star_field_serde_roundtrip() {
        let field = StarField {
            stars: vec![StarViz {
                position: [0.0, 0.0, 10.0],
                magnitude: -1.46,
                color_index: 0.0,
                radius: 0.5,
            }],
        };
        let json = serde_json::to_string(&field).unwrap();
        let back: StarField = serde_json::from_str(&json).unwrap();
        assert_eq!(back, field);
    }

    #[test]
    fn spectral_profile_serde_roundtrip() {
        let sp = SpectralProfile {
            wavelengths: vec![400.0, 500.0, 600.0, 700.0],
            intensities: vec![0.8, 0.95, 0.9, 0.7],
            absorption_lines: vec![486.1, 656.3],
        };
        let json = serde_json::to_string(&sp).unwrap();
        let back: SpectralProfile = serde_json::from_str(&json).unwrap();
        assert_eq!(back, sp);
    }
}
