//! Hertzsprung-Russell diagram classification and stellar taxonomy.
//!
//! Provides spectral classification from temperature, luminosity class
//! determination, and HR diagram region identification.

use serde::{Deserialize, Serialize};

use crate::constants;
use crate::star::SpectralClass;

/// HR diagram position.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[non_exhaustive]
pub struct HrPosition {
    /// Effective temperature in Kelvin.
    pub temperature_k: f64,
    /// Absolute magnitude.
    pub absolute_magnitude: f64,
}

/// Yerkes luminosity class (Morgan-Keenan system).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
#[non_exhaustive]
pub enum LuminosityClass {
    /// Ia — luminous supergiant.
    Ia,
    /// Ib — less luminous supergiant.
    Ib,
    /// II — bright giant.
    II,
    /// III — normal giant.
    III,
    /// IV — subgiant.
    IV,
    /// V — main sequence (dwarf).
    V,
    /// VI — subdwarf.
    VI,
    /// VII — white dwarf.
    VII,
}

/// Region on the Hertzsprung-Russell diagram.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
#[non_exhaustive]
pub enum HrRegion {
    /// Upper left — hot, luminous (O/B supergiants).
    MainSequenceHot,
    /// Middle — solar-type main sequence.
    MainSequenceMid,
    /// Lower right — cool, dim (M dwarfs).
    MainSequenceCool,
    /// Upper right — red giants.
    RedGiantBranch,
    /// Horizontal branch / red clump.
    HorizontalBranch,
    /// Asymptotic giant branch.
    AsymptoticGiantBranch,
    /// Lower left — white dwarfs.
    WhiteDwarfRegion,
    /// Supergiant region.
    Supergiant,
}

/// Determine spectral class from effective temperature (K).
///
/// Uses standard Morgan-Keenan temperature boundaries.
#[must_use]
pub fn spectral_class_from_temperature(temperature_k: f64) -> SpectralClass {
    if temperature_k >= constants::T_W_MIN {
        SpectralClass::W
    } else if temperature_k >= constants::T_O_MIN {
        SpectralClass::O
    } else if temperature_k >= constants::T_B_MIN {
        SpectralClass::B
    } else if temperature_k >= constants::T_A_MIN {
        SpectralClass::A
    } else if temperature_k >= constants::T_F_MIN {
        SpectralClass::F
    } else if temperature_k >= constants::T_G_MIN {
        SpectralClass::G
    } else if temperature_k >= constants::T_K_MIN {
        SpectralClass::K
    } else if temperature_k >= constants::T_M_MIN {
        SpectralClass::M
    } else if temperature_k >= constants::T_L_MIN {
        SpectralClass::L
    } else if temperature_k >= constants::T_T_MIN {
        SpectralClass::T
    } else {
        SpectralClass::Y
    }
}

/// Determine spectral subclass (0–9) from temperature within a spectral type.
///
/// 0 = hottest end of the class, 9 = coolest. Linear interpolation
/// within each class boundary.
#[must_use]
pub fn spectral_subclass(temperature_k: f64) -> u8 {
    let (t_max, t_min) = class_temp_range(spectral_class_from_temperature(temperature_k));

    if t_max <= t_min {
        return 0;
    }

    let frac = (t_max - temperature_k) / (t_max - t_min);
    let sub = (frac * 10.0).clamp(0.0, 9.0) as u8;
    sub.min(9)
}

/// Temperature range for a spectral class, returned as (hottest, coolest).
///
/// The first element is the upper (hotter) boundary, the second is the lower
/// (cooler) boundary. This ordering matches the subclass convention where
/// subclass 0 is the hottest end.
#[must_use]
fn class_temp_range(class: SpectralClass) -> (f64, f64) {
    match class {
        SpectralClass::W => (100_000.0, constants::T_W_MIN),
        SpectralClass::O => (constants::T_W_MIN, constants::T_O_MIN),
        SpectralClass::B => (constants::T_O_MIN, constants::T_B_MIN),
        SpectralClass::A => (constants::T_B_MIN, constants::T_A_MIN),
        SpectralClass::F => (constants::T_A_MIN, constants::T_F_MIN),
        SpectralClass::G => (constants::T_F_MIN, constants::T_G_MIN),
        SpectralClass::K => (constants::T_G_MIN, constants::T_K_MIN),
        SpectralClass::M => (constants::T_K_MIN, constants::T_M_MIN),
        SpectralClass::L => (constants::T_M_MIN, constants::T_L_MIN),
        SpectralClass::T => (constants::T_L_MIN, constants::T_T_MIN),
        SpectralClass::Y => (constants::T_T_MIN, 250.0),
    }
}

/// Format a full MK spectral classification string (e.g., "G2V").
///
/// Combines spectral class, subclass, and luminosity class.
#[must_use]
pub fn format_classification(temperature_k: f64, luminosity_class: LuminosityClass) -> String {
    let class = spectral_class_from_temperature(temperature_k);
    let sub = spectral_subclass(temperature_k);
    let lc = match luminosity_class {
        LuminosityClass::Ia => "Ia",
        LuminosityClass::Ib => "Ib",
        LuminosityClass::II => "II",
        LuminosityClass::III => "III",
        LuminosityClass::IV => "IV",
        LuminosityClass::V => "V",
        LuminosityClass::VI => "VI",
        LuminosityClass::VII => "VII",
    };
    format!("{class}{sub}{lc}")
}

/// Estimate luminosity class from surface gravity log(g) (CGS).
///
/// Simplified mapping based on typical log(g) ranges:
/// - log(g) < 1.0 → Ia supergiant
/// - 1.0–2.0 → II bright giant
/// - 2.0–3.5 → III giant
/// - 3.5–4.0 → IV subgiant
/// - 4.0–5.5 → V main sequence
/// - 5.5–7.0 → VI subdwarf
/// - > 7.0 → VII white dwarf
#[must_use]
pub fn luminosity_class_from_log_g(log_g: f64) -> LuminosityClass {
    if log_g < 1.0 {
        LuminosityClass::Ia
    } else if log_g < 2.0 {
        LuminosityClass::II
    } else if log_g < 3.5 {
        LuminosityClass::III
    } else if log_g < 4.0 {
        LuminosityClass::IV
    } else if log_g < 5.5 {
        LuminosityClass::V
    } else if log_g < 7.0 {
        LuminosityClass::VI
    } else {
        LuminosityClass::VII
    }
}

/// Determine HR diagram region from temperature and luminosity (solar units).
#[must_use]
pub fn hr_region(temperature_k: f64, luminosity_solar: f64) -> HrRegion {
    let abs_mag = crate::luminosity::absolute_bolometric_magnitude(luminosity_solar);

    // White dwarf region: hot but very dim
    if abs_mag > 10.0 {
        return HrRegion::WhiteDwarfRegion;
    }

    // Supergiants: very luminous
    if abs_mag < -4.0 {
        return HrRegion::Supergiant;
    }

    // Giant region: luminous + cool
    if temperature_k < 5500.0 && abs_mag < 2.0 {
        if abs_mag < 0.0 {
            return HrRegion::AsymptoticGiantBranch;
        }
        return HrRegion::RedGiantBranch;
    }

    // Horizontal branch: intermediate luminosity, moderate temperature
    if abs_mag > -1.0 && abs_mag < 2.0 && temperature_k > 5500.0 && temperature_k < 8000.0 {
        return HrRegion::HorizontalBranch;
    }

    // Main sequence
    if temperature_k >= 10_000.0 {
        HrRegion::MainSequenceHot
    } else if temperature_k >= 3700.0 {
        HrRegion::MainSequenceMid
    } else {
        HrRegion::MainSequenceCool
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn hr_position_serde_roundtrip() {
        let pos = HrPosition {
            temperature_k: 5772.0,
            absolute_magnitude: 4.83,
        };
        let json = serde_json::to_string(&pos).unwrap();
        let back: HrPosition = serde_json::from_str(&json).unwrap();
        assert!((back.temperature_k - 5772.0).abs() < f64::EPSILON);
        assert!((back.absolute_magnitude - 4.83).abs() < f64::EPSILON);
    }

    #[test]
    fn sun_spectral_class() {
        let class = spectral_class_from_temperature(constants::T_SUN);
        assert_eq!(class, SpectralClass::G);
    }

    #[test]
    fn hot_star_class() {
        assert_eq!(spectral_class_from_temperature(35_000.0), SpectralClass::O);
        assert_eq!(spectral_class_from_temperature(15_000.0), SpectralClass::B);
    }

    #[test]
    fn cool_star_class() {
        assert_eq!(spectral_class_from_temperature(3000.0), SpectralClass::M);
        assert_eq!(spectral_class_from_temperature(4000.0), SpectralClass::K);
    }

    #[test]
    fn wolf_rayet_class() {
        assert_eq!(spectral_class_from_temperature(60_000.0), SpectralClass::W);
        assert_eq!(spectral_class_from_temperature(50_000.0), SpectralClass::W);
    }

    #[test]
    fn brown_dwarf_classes() {
        assert_eq!(spectral_class_from_temperature(1_800.0), SpectralClass::L);
        assert_eq!(spectral_class_from_temperature(1_000.0), SpectralClass::T);
        assert_eq!(spectral_class_from_temperature(400.0), SpectralClass::Y);
    }

    #[test]
    fn solar_subclass() {
        // Sun is G2 (T = 5772 K)
        let sub = spectral_subclass(constants::T_SUN);
        // G range: 6000–5200 K, Sun at 5772 → frac = (6000-5772)/(6000-5200) = 228/800 = 0.285 → sub=2
        assert_eq!(sub, 2, "Solar subclass: G{sub}, expected G2");
    }

    #[test]
    fn solar_classification_string() {
        let s = format_classification(constants::T_SUN, LuminosityClass::V);
        assert_eq!(s, "G2V");
    }

    #[test]
    fn solar_luminosity_class() {
        let lc = luminosity_class_from_log_g(4.44);
        assert_eq!(lc, LuminosityClass::V);
    }

    #[test]
    fn giant_luminosity_class() {
        let lc = luminosity_class_from_log_g(2.5);
        assert_eq!(lc, LuminosityClass::III);
    }

    #[test]
    fn white_dwarf_luminosity_class() {
        let lc = luminosity_class_from_log_g(8.0);
        assert_eq!(lc, LuminosityClass::VII);
    }

    #[test]
    fn sun_hr_region() {
        let region = hr_region(constants::T_SUN, 1.0);
        assert_eq!(region, HrRegion::MainSequenceMid);
    }

    #[test]
    fn luminosity_class_serde_roundtrip() {
        for lc in [
            LuminosityClass::Ia,
            LuminosityClass::Ib,
            LuminosityClass::II,
            LuminosityClass::III,
            LuminosityClass::IV,
            LuminosityClass::V,
            LuminosityClass::VI,
            LuminosityClass::VII,
        ] {
            let json = serde_json::to_string(&lc).unwrap();
            let back: LuminosityClass = serde_json::from_str(&json).unwrap();
            assert_eq!(back, lc);
        }
    }

    #[test]
    fn hr_region_serde_roundtrip() {
        for r in [
            HrRegion::MainSequenceHot,
            HrRegion::MainSequenceMid,
            HrRegion::MainSequenceCool,
            HrRegion::RedGiantBranch,
            HrRegion::WhiteDwarfRegion,
            HrRegion::Supergiant,
        ] {
            let json = serde_json::to_string(&r).unwrap();
            let back: HrRegion = serde_json::from_str(&json).unwrap();
            assert_eq!(back, r);
        }
    }
}
