//! Analytic stellar evolution fitting formulae.
//!
//! Implements the ZAMS (zero-age main sequence) luminosity and radius fits
//! from Tout et al. (1996, MNRAS 281, 257) and the main-sequence lifetime
//! from Hurley, Pols & Tout (2000, MNRAS 315, 543).
//!
//! All formulae are parameterised by metallicity Z through ζ = log₁₀(Z/0.02).
//! Valid for M ∈ \[0.1, 100\] M☉ and Z ∈ \[0.0001, 0.03\].

use serde::{Deserialize, Serialize};

/// Solar metallicity (mass fraction of metals).
pub const Z_SUN: f64 = 0.02;

/// Metallicity-dependent coefficient: α + βζ + γζ² + δζ³ + εζ⁴.
#[must_use]
fn coeff(c: &[f64; 5], zeta: f64) -> f64 {
    c[0] + zeta * (c[1] + zeta * (c[2] + zeta * (c[3] + zeta * c[4])))
}

/// Compute ζ = log₁₀(Z / 0.02).
#[must_use]
#[inline]
fn zeta(z: f64) -> f64 {
    (z / Z_SUN).log10()
}

// ── Tout et al. (1996) Table 1 — ZAMS luminosity coefficients ────────────

const A1: [f64; 5] = [
    0.397_041_70,
    -0.329_135_74,
    0.347_766_88,
    0.374_708_51,
    0.090_119_15,
];
const A2: [f64; 5] = [
    8.527_626_00,
    -24.412_259_73,
    56.435_971_07,
    37.061_525_75,
    5.456_240_60,
];
const A3: [f64; 5] = [
    0.000_255_46,
    -0.001_234_61,
    -0.000_232_46,
    0.000_455_19,
    0.000_161_76,
];
const A4: [f64; 5] = [
    5.432_889_00,
    -8.621_578_06,
    13.442_020_49,
    14.515_841_35,
    3.397_930_84,
];
const A5: [f64; 5] = [
    5.563_579_00,
    -10.323_452_24,
    19.443_229_80,
    18.973_613_47,
    4.169_030_97,
];
const A6: [f64; 5] = [
    0.788_660_60,
    -2.908_709_42,
    6.547_135_31,
    4.056_066_57,
    0.532_873_22,
];
const A7: [f64; 5] = [
    0.005_866_85,
    -0.017_042_37,
    0.038_723_48,
    0.025_700_41,
    0.003_833_76,
];

// ── Tout et al. (1996) Table 2 — ZAMS radius coefficients ────────────────

const A8: [f64; 5] = [
    1.715_359_00,
    0.622_462_12,
    -0.925_577_61,
    -1.169_969_66,
    -0.306_314_91,
];
const A9: [f64; 5] = [
    6.597_788_00,
    -0.424_500_44,
    -12.133_394_27,
    -10.735_094_84,
    -2.514_870_77,
];
const A10: [f64; 5] = [
    10.088_550_00,
    -7.117_270_86,
    31.671_194_79,
    24.248_483_22,
    3.781_420_98,
];
const A11: [f64; 5] = [
    1.012_495_00,
    0.326_996_90,
    -0.009_234_18,
    -0.038_768_58,
    -0.004_127_50,
];
const A12: [f64; 5] = [
    0.074_901_66,
    0.024_104_13,
    0.072_336_64,
    0.030_404_67,
    0.001_977_41,
];
const A13: [f64; 5] = [0.010_774_22, 0.0, 0.0, 0.0, 0.0];
const A14: [f64; 5] = [
    3.082_234_00,
    0.944_720_50,
    -2.152_008_82,
    -2.492_194_96,
    -0.638_487_38,
];
const A15: [f64; 5] = [
    17.847_780_00,
    -7.453_456_90,
    48.250_603_16,
    37.269_119_24,
    5.995_741_20,
];
const A16: [f64; 5] = [
    0.000_225_82,
    -0.001_868_99,
    0.003_887_83,
    0.001_424_02,
    -0.000_076_71,
];

/// ZAMS luminosity (L☉) from mass (M☉) and metallicity Z.
///
/// Tout et al. (1996) Eq. 1. Valid for M ∈ \[0.1, 100\] M☉, Z ∈ \[0.0001, 0.03\].
#[must_use]
pub fn zams_luminosity(mass: f64, z: f64) -> f64 {
    if mass <= 0.0 {
        return 0.0;
    }
    let ze = zeta(z.clamp(0.0001, 0.03));
    let a1 = coeff(&A1, ze);
    let a2 = coeff(&A2, ze);
    let a3 = coeff(&A3, ze);
    let a4 = coeff(&A4, ze);
    let a5 = coeff(&A5, ze);
    let a6 = coeff(&A6, ze);
    let a7 = coeff(&A7, ze);

    let m55 = mass.powf(5.5);
    let m11 = mass.powi(11);
    let m3 = mass.powi(3);
    let m5 = mass.powi(5);
    let m7 = mass.powi(7);
    let m8 = mass.powi(8);
    let m95 = mass.powf(9.5);

    let num = a1 * m55 + a2 * m11;
    let den = a3 + m3 + a4 * m5 + a5 * m7 + a6 * m8 + a7 * m95;

    if den <= 0.0 { 0.0 } else { num / den }
}

/// ZAMS radius (R☉) from mass (M☉) and metallicity Z.
///
/// Tout et al. (1996) Eq. 2. Valid for M ∈ \[0.1, 100\] M☉, Z ∈ \[0.0001, 0.03\].
#[must_use]
pub fn zams_radius(mass: f64, z: f64) -> f64 {
    if mass <= 0.0 {
        return 0.0;
    }
    let ze = zeta(z.clamp(0.0001, 0.03));
    let a8 = coeff(&A8, ze);
    let a9 = coeff(&A9, ze);
    let a10 = coeff(&A10, ze);
    let a11 = coeff(&A11, ze);
    let a12 = coeff(&A12, ze);
    let a13 = coeff(&A13, ze);
    let a14 = coeff(&A14, ze);
    let a15 = coeff(&A15, ze);
    let a16 = coeff(&A16, ze);

    let m25 = mass.powf(2.5);
    let m65 = mass.powf(6.5);
    let m11 = mass.powi(11);
    let m19 = mass.powi(19);
    let m195 = mass.powf(19.5);
    let m2 = mass * mass;
    let m85 = mass.powf(8.5);
    let m185 = mass.powf(18.5);

    let num = a8 * m25 + a9 * m65 + a10 * m11 + a11 * m19 + a12 * m195;
    let den = a13 + a14 * m2 + a15 * m85 + m185 + a16 * m195;

    if den <= 0.0 { 0.0 } else { num / den }
}

/// ZAMS effective temperature (K) from mass and metallicity.
///
/// Derived from L_ZAMS and R_ZAMS via Stefan-Boltzmann:
/// T = T☉ × (L/L☉)^0.25 / (R/R☉)^0.5.
#[must_use]
pub fn zams_temperature(mass: f64, z: f64) -> f64 {
    let l = zams_luminosity(mass, z);
    let r = zams_radius(mass, z);
    if l <= 0.0 || r <= 0.0 {
        return 0.0;
    }
    crate::constants::T_SUN * l.powf(0.25) / r.sqrt()
}

/// ZAMS properties: (luminosity L☉, radius R☉, temperature K).
///
/// Convenience function returning all three ZAMS anchors.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct ZamsProperties {
    /// Luminosity in solar luminosities.
    pub luminosity_solar: f64,
    /// Radius in solar radii.
    pub radius_solar: f64,
    /// Effective temperature in Kelvin.
    pub temperature_k: f64,
}

/// Compute all ZAMS properties for a given mass and metallicity.
#[must_use]
pub fn zams_properties(mass: f64, z: f64) -> ZamsProperties {
    let luminosity_solar = zams_luminosity(mass, z);
    let radius_solar = zams_radius(mass, z);
    let temperature_k = if luminosity_solar > 0.0 && radius_solar > 0.0 {
        crate::constants::T_SUN * luminosity_solar.powf(0.25) / radius_solar.sqrt()
    } else {
        0.0
    };
    ZamsProperties {
        luminosity_solar,
        radius_solar,
        temperature_k,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn zams_luminosity_solar() {
        // 1 M_sun at Z=0.02: ZAMS luminosity ~0.7 L_sun (faint young Sun)
        let l = zams_luminosity(1.0, Z_SUN);
        assert!(
            (l - 0.74).abs() < 0.1,
            "Solar ZAMS luminosity: {l}, expected ~0.74"
        );
    }

    #[test]
    fn zams_luminosity_increases_with_mass() {
        let l1 = zams_luminosity(1.0, Z_SUN);
        let l5 = zams_luminosity(5.0, Z_SUN);
        let l20 = zams_luminosity(20.0, Z_SUN);
        assert!(l5 > l1, "5 M_sun should be brighter than 1 M_sun");
        assert!(l20 > l5, "20 M_sun should be brighter than 5 M_sun");
    }

    #[test]
    fn zams_luminosity_low_mass() {
        let l = zams_luminosity(0.1, Z_SUN);
        assert!(l > 0.0 && l < 0.01, "0.1 M_sun ZAMS L: {l}");
    }

    #[test]
    fn zams_radius_solar() {
        // 1 M_sun at Z=0.02 should give ~1 R_sun
        let r = zams_radius(1.0, Z_SUN);
        assert!(
            (r - 1.0).abs() < 0.15,
            "Solar ZAMS radius: {r}, expected ~1.0"
        );
    }

    #[test]
    fn zams_radius_increases_with_mass() {
        let r1 = zams_radius(1.0, Z_SUN);
        let r5 = zams_radius(5.0, Z_SUN);
        let r20 = zams_radius(20.0, Z_SUN);
        assert!(r5 > r1);
        assert!(r20 > r5);
    }

    #[test]
    fn zams_temperature_solar() {
        let t = zams_temperature(1.0, Z_SUN);
        // Solar ZAMS T should be near 5772 K (within ~10%)
        assert!(
            (t - 5772.0).abs() < 600.0,
            "Solar ZAMS T: {t} K, expected ~5772"
        );
    }

    #[test]
    fn zams_temperature_hot_for_massive() {
        let t1 = zams_temperature(1.0, Z_SUN);
        let t20 = zams_temperature(20.0, Z_SUN);
        assert!(t20 > t1, "Massive stars should be hotter on ZAMS");
    }

    #[test]
    fn zams_metallicity_effect() {
        // Lower metallicity → slightly higher luminosity (less opacity)
        let l_low_z = zams_luminosity(1.0, 0.001);
        let l_solar = zams_luminosity(1.0, Z_SUN);
        assert!(
            l_low_z > l_solar,
            "Low-Z should be brighter: Z=0.001 → L={l_low_z}, Z=0.02 → L={l_solar}"
        );
    }

    #[test]
    fn zams_properties_roundtrip() {
        let p = zams_properties(1.0, Z_SUN);
        let json = serde_json::to_string(&p).unwrap();
        let back: ZamsProperties = serde_json::from_str(&json).unwrap();
        assert_eq!(back, p);
    }

    #[test]
    fn zams_zero_mass() {
        assert!(zams_luminosity(0.0, Z_SUN).abs() < f64::EPSILON);
        assert!(zams_radius(0.0, Z_SUN).abs() < f64::EPSILON);
        assert!(zams_temperature(0.0, Z_SUN).abs() < f64::EPSILON);
    }
}
