//! Soorat integration — visualization data structures and helpers for stellar astrophysics.
//!
//! Provides conversion from [`Star`] instances into rendering-ready
//! data structures for soorat visualization.

use serde::{Deserialize, Serialize};

use crate::classification;
use crate::luminosity;
use crate::spectral;
use crate::star::Star;

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

impl HrDiagramPoint {
    /// Create an HR diagram point from a [`Star`].
    #[must_use]
    pub fn from_star(star: &Star) -> Self {
        Self {
            temperature: star.temperature_k,
            luminosity: star.luminosity_solar,
            spectral_class: classification::format_classification(
                star.temperature_k,
                star.luminosity_class,
            ),
            color: temperature_to_rgb(star.temperature_k),
        }
    }
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

impl EvolutionTrack {
    /// Create a new empty evolution track with the given label.
    #[must_use]
    pub fn new(label: impl Into<String>) -> Self {
        Self {
            points: Vec::new(),
            ages: Vec::new(),
            label: label.into(),
        }
    }

    /// Append a snapshot from a [`Star`] at a given age.
    pub fn push_snapshot(&mut self, star: &Star, age_years: f64) {
        self.points
            .push([star.temperature_k, star.luminosity_solar, star.radius_solar]);
        self.ages.push(age_years);
    }

    /// Number of snapshots in the track.
    #[must_use]
    pub fn len(&self) -> usize {
        self.ages.len()
    }

    /// Whether the track is empty.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.ages.is_empty()
    }

    /// Generate a main-sequence evolution track using SSE fitting formulae.
    ///
    /// Samples `n_points` evenly-spaced steps from ZAMS (τ=0) to TMS (τ=1)
    /// for a star of given mass and metallicity.
    #[must_use]
    pub fn from_sse_ms(mass: f64, z: f64, label: impl Into<String>, n_points: usize) -> Self {
        use crate::sse;

        let mut track = Self::new(label);
        if n_points == 0 || mass <= 0.0 {
            return track;
        }

        let t_ms = sse::ms_lifetime(mass, z);

        for i in 0..n_points {
            let tau = if n_points > 1 {
                i as f64 / (n_points - 1) as f64
            } else {
                0.0
            };
            let p = sse::ms_properties(mass, z, tau * t_ms);
            track
                .points
                .push([p.temperature_k, p.luminosity_solar, p.radius_solar]);
            track.ages.push(tau * t_ms);
        }

        track
    }
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

impl SpectralProfile {
    /// Generate a blackbody spectral profile for a [`Star`] over a wavelength range.
    ///
    /// Samples `n_points` wavelengths from `lambda_min_nm` to `lambda_max_nm`,
    /// normalized so the peak intensity is 1.0.
    ///
    /// `absorption_lines` is a list of absorption line positions (nm) to include.
    #[must_use]
    pub fn from_star_blackbody(
        star: &Star,
        lambda_min_nm: f64,
        lambda_max_nm: f64,
        n_points: usize,
        absorption_lines: Vec<f64>,
    ) -> Self {
        if n_points < 2 || lambda_max_nm <= lambda_min_nm {
            return Self {
                wavelengths: Vec::new(),
                intensities: Vec::new(),
                absorption_lines,
            };
        }

        let step = (lambda_max_nm - lambda_min_nm) / (n_points - 1) as f64;

        let wavelengths: Vec<f64> = (0..n_points)
            .map(|i| lambda_min_nm + i as f64 * step)
            .collect();

        let raw: Vec<f64> = wavelengths
            .iter()
            .map(|&lam| spectral::planck_radiance_nm(lam, star.temperature_k))
            .collect();

        let max = raw.iter().cloned().fold(0.0_f64, f64::max);

        let intensities = if max > 0.0 {
            raw.iter().map(|&v| v / max).collect()
        } else {
            vec![0.0; n_points]
        };

        Self {
            wavelengths,
            intensities,
            absorption_lines,
        }
    }
}

/// Star field for instanced point/billboard rendering.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[non_exhaustive]
pub struct StarField {
    /// Stars with position and properties.
    pub stars: Vec<StarViz>,
}

impl StarField {
    /// Create an empty star field.
    #[must_use]
    pub fn new() -> Self {
        Self { stars: Vec::new() }
    }

    /// Add a star at the given position and distance.
    ///
    /// Computes apparent magnitude and B-V color index from the star's properties.
    pub fn push(&mut self, star: &Star, position: [f32; 3], distance_pc: f64) {
        self.stars
            .push(StarViz::from_star(star, position, distance_pc));
    }

    /// Number of stars in the field.
    #[must_use]
    pub fn len(&self) -> usize {
        self.stars.len()
    }

    /// Whether the field is empty.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.stars.is_empty()
    }
}

impl Default for StarField {
    fn default() -> Self {
        Self::new()
    }
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

impl StarViz {
    /// Create a `StarViz` from a [`Star`], position, and distance.
    #[must_use]
    pub fn from_star(star: &Star, position: [f32; 3], distance_pc: f64) -> Self {
        let abs_bol = luminosity::absolute_bolometric_magnitude(star.luminosity_solar);
        let bc = luminosity::bolometric_correction(star.temperature_k);
        let abs_v = abs_bol - bc;
        let app_mag = luminosity::apparent_magnitude(abs_v, distance_pc);
        let bv = temperature_to_bv(star.temperature_k);

        Self {
            position,
            magnitude: app_mag as f32,
            color_index: bv as f32,
            radius: star.radius_solar as f32,
        }
    }
}

// ── Color conversions ─────────────────────────────────────────────────────

/// Approximate RGB color from stellar effective temperature.
///
/// Based on Ballesteros (2012) CIE chromaticity approximation for blackbody
/// radiators, simplified for rendering use. Returns sRGB in `[0.0, 1.0]`.
#[must_use]
pub fn temperature_to_rgb(temperature_k: f64) -> [f32; 3] {
    // Normalize to thousands of Kelvin, clamp to useful range
    let t = (temperature_k / 100.0).clamp(20.0, 400.0);

    let r = if t <= 66.0 {
        1.0
    } else {
        let x = t - 60.0;
        (329.698_727_446 * x.powf(-0.133_204_759_2) / 255.0).clamp(0.0, 1.0)
    };

    let g = if t <= 66.0 {
        let x = t;
        (99.470_802_586_1 * x.ln() - 161.119_568_166_1) / 255.0
    } else {
        let x = t - 60.0;
        288.122_169_528_3 * x.powf(-0.075_514_849_2) / 255.0
    }
    .clamp(0.0, 1.0);

    let b = if t >= 66.0 {
        1.0
    } else if t <= 19.0 {
        0.0
    } else {
        let x = t - 10.0;
        (138.517_731_223_1 * x.ln() - 305.044_792_730_7) / 255.0
    }
    .clamp(0.0, 1.0);

    [r as f32, g as f32, b as f32]
}

/// Approximate B-V color index from effective temperature.
///
/// Inverts the Ballesteros (2012) relation:
/// T = 4600 × (1/(0.92(B−V)+1.7) + 1/(0.92(B−V)+0.62))
///
/// Uses Newton's method to solve for B-V given T. Valid ~3000–40000 K.
#[must_use]
pub fn temperature_to_bv(temperature_k: f64) -> f64 {
    if temperature_k <= 0.0 {
        return 2.0;
    }

    // Ballesteros T(B-V) function
    let t_from_bv =
        |bv: f64| -> f64 { 4600.0 * (1.0 / (0.92 * bv + 1.7) + 1.0 / (0.92 * bv + 0.62)) };

    // Derivative dT/d(B-V)
    let dt_dbv = |bv: f64| -> f64 {
        let a = 0.92 * bv + 1.7;
        let b = 0.92 * bv + 0.62;
        -4600.0 * 0.92 * (1.0 / (a * a) + 1.0 / (b * b))
    };

    // Newton iteration: solve T(bv) - target = 0
    let mut bv = 0.65; // start near solar
    for _ in 0..20 {
        let f = t_from_bv(bv) - temperature_k;
        if f.abs() < 0.1 {
            break; // converged to <0.1 K
        }
        let df = dt_dbv(bv);
        if df.abs() < 1e-30 {
            break;
        }
        bv -= f / df;
    }

    bv.clamp(-0.4, 2.0)
}

#[cfg(test)]
mod tests {
    use super::*;
    fn sun() -> Star {
        Star::sun().unwrap()
    }

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
    fn hr_point_from_star() {
        let s = sun();
        let pt = HrDiagramPoint::from_star(&s);
        assert_eq!(pt.spectral_class, "G2V");
        assert!((pt.temperature - 5772.0).abs() < f64::EPSILON);
        assert!((pt.luminosity - 1.0).abs() < 0.02);
        // RGB should be warm yellowish
        assert!(pt.color[0] > 0.8, "R={}", pt.color[0]);
        assert!(pt.color[1] > 0.7, "G={}", pt.color[1]);
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
    fn evolution_track_push() {
        let s = sun();
        let mut track = EvolutionTrack::new("Sun");
        assert!(track.is_empty());
        track.push_snapshot(&s, 4.6e9);
        assert_eq!(track.len(), 1);
        assert!((track.ages[0] - 4.6e9).abs() < f64::EPSILON);
        assert!((track.points[0][0] - 5772.0).abs() < f64::EPSILON);
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
    fn star_field_push() {
        let s = sun();
        let mut field = StarField::new();
        assert!(field.is_empty());
        field.push(&s, [0.0, 0.0, 0.0], 10.0);
        assert_eq!(field.len(), 1);
        // At 10 pc, apparent ~= absolute
        assert!(field.stars[0].magnitude.is_finite());
        assert!(field.stars[0].radius > 0.0);
    }

    #[test]
    fn evolution_track_from_sse() {
        let track = EvolutionTrack::from_sse_ms(1.0, 0.02, "Solar", 20);
        assert_eq!(track.len(), 20);
        assert_eq!(track.label, "Solar");
        // First point should be ZAMS, last TMS
        let t_first = track.points[0][0]; // temperature
        let t_last = track.points[19][0];
        let l_first = track.points[0][1]; // luminosity
        let l_last = track.points[19][1];
        // Star should brighten over MS
        assert!(l_last > l_first, "L_TMS > L_ZAMS: {l_last} vs {l_first}");
        // Temperatures should be reasonable
        assert!(t_first > 4000.0 && t_first < 8000.0);
        assert!(t_last > 4000.0 && t_last < 8000.0);
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

    #[test]
    fn spectral_profile_from_star() {
        let s = sun();
        let sp = SpectralProfile::from_star_blackbody(&s, 300.0, 900.0, 100, vec![486.1, 656.3]);
        assert_eq!(sp.wavelengths.len(), 100);
        assert_eq!(sp.intensities.len(), 100);
        // Peak should be normalized to 1.0
        let max = sp.intensities.iter().cloned().fold(0.0_f64, f64::max);
        assert!((max - 1.0).abs() < f64::EPSILON, "max intensity: {max}");
        // Absorption lines preserved
        assert_eq!(sp.absorption_lines, vec![486.1, 656.3]);
    }

    #[test]
    fn spectral_profile_empty_on_bad_input() {
        let s = sun();
        let sp = SpectralProfile::from_star_blackbody(&s, 500.0, 400.0, 100, vec![]);
        assert!(sp.wavelengths.is_empty());
        let sp2 = SpectralProfile::from_star_blackbody(&s, 400.0, 500.0, 1, vec![]);
        assert!(sp2.wavelengths.is_empty());
    }

    #[test]
    fn star_viz_from_star() {
        let s = sun();
        let viz = StarViz::from_star(&s, [1.0, 2.0, 3.0], 10.0);
        assert_eq!(viz.position, [1.0, 2.0, 3.0]);
        assert!(viz.magnitude.is_finite());
        assert!(viz.color_index.is_finite());
        assert!((viz.radius - 1.0).abs() < f32::EPSILON);
    }

    #[test]
    fn temperature_to_rgb_hot_is_blue() {
        let [r, _, b] = temperature_to_rgb(30_000.0);
        assert!(b > r, "Hot star should be bluer: r={r}, b={b}");
    }

    #[test]
    fn temperature_to_rgb_cool_is_red() {
        let [r, _, b] = temperature_to_rgb(3000.0);
        assert!(r > b, "Cool star should be redder: r={r}, b={b}");
    }

    #[test]
    fn temperature_to_bv_sun() {
        let bv = temperature_to_bv(5772.0);
        // Solar B-V ≈ 0.65
        assert!((bv - 0.65).abs() < 0.15, "Solar B-V: {bv}, expected ~0.65");
    }

    #[test]
    fn temperature_to_bv_hot() {
        let bv = temperature_to_bv(30_000.0);
        assert!(bv < 0.0, "Hot star B-V should be negative: {bv}");
    }

    #[test]
    fn temperature_to_bv_cool() {
        let bv = temperature_to_bv(3500.0);
        assert!(bv > 1.0, "Cool star B-V should be > 1.0: {bv}");
    }

    #[test]
    fn star_field_default() {
        let field = StarField::default();
        assert!(field.is_empty());
    }
}
