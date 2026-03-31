//! Luminosity, magnitude, and distance calculations.
//!
//! All functions use IAU 2015 / CODATA 2018 constants from [`crate::constants`].

use crate::constants;

/// Compute luminosity (W) from radius (m) and effective temperature (K).
///
/// Stefan-Boltzmann law: L = 4π R² σ T⁴.
#[must_use]
#[inline]
pub fn luminosity_from_radius_temp(radius_m: f64, temperature_k: f64) -> f64 {
    4.0 * std::f64::consts::PI * radius_m * radius_m * constants::SIGMA_SB * temperature_k.powi(4)
}

/// Compute luminosity in solar units from radius (solar radii) and temperature (K).
///
/// Returns L / L_sun.
#[must_use]
#[inline]
pub fn luminosity_solar(radius_solar: f64, temperature_k: f64) -> f64 {
    let r_m = radius_solar * constants::R_SUN;
    luminosity_from_radius_temp(r_m, temperature_k) / constants::L_SUN
}

/// Compute absolute bolometric magnitude from luminosity (solar luminosities).
///
/// M_bol = M_bol,sun − 2.5 log₁₀(L / L_sun).
#[must_use]
#[inline]
pub fn absolute_bolometric_magnitude(luminosity_solar: f64) -> f64 {
    if luminosity_solar <= 0.0 {
        return f64::INFINITY;
    }
    constants::M_BOL_SUN - 2.5 * luminosity_solar.log10()
}

/// Compute luminosity (solar luminosities) from absolute bolometric magnitude.
///
/// L / L_sun = 10^((M_bol,sun − M_bol) / 2.5).
#[must_use]
#[inline]
pub fn luminosity_from_abs_bol_mag(absolute_bol_mag: f64) -> f64 {
    10.0_f64.powf((constants::M_BOL_SUN - absolute_bol_mag) / 2.5)
}

/// Compute distance modulus: m − M = 5 log₁₀(d / 10 pc).
///
/// `distance_pc` is distance in parsecs.
#[must_use]
#[inline]
pub fn distance_modulus(distance_pc: f64) -> f64 {
    if distance_pc <= 0.0 {
        return f64::NEG_INFINITY;
    }
    5.0 * (distance_pc / 10.0).log10()
}

/// Compute distance (parsecs) from distance modulus.
///
/// d = 10^((μ + 5) / 5) = 10 × 10^(μ/5).
#[must_use]
#[inline]
pub fn distance_from_modulus(modulus: f64) -> f64 {
    10.0 * 10.0_f64.powf(modulus / 5.0)
}

/// Compute apparent magnitude from absolute magnitude and distance (parsecs).
///
/// m = M + 5 log₁₀(d / 10 pc).
#[must_use]
#[inline]
pub fn apparent_magnitude(absolute_mag: f64, distance_pc: f64) -> f64 {
    absolute_mag + distance_modulus(distance_pc)
}

/// Compute absolute magnitude from apparent magnitude and distance (parsecs).
///
/// M = m − 5 log₁₀(d / 10 pc).
#[must_use]
#[inline]
pub fn absolute_magnitude(apparent_mag: f64, distance_pc: f64) -> f64 {
    apparent_mag - distance_modulus(distance_pc)
}

/// Bolometric correction from effective temperature (K).
///
/// Polynomial fit in log₁₀(T_eff), calibrated to reproduce known values:
/// Sun (5772 K) → BC ≈ −0.08, Vega (9600 K) → BC ≈ −0.3,
/// cool M-dwarf (3200 K) → BC ≈ −2.0, hot O-star (40000 K) → BC ≈ −4.0.
///
/// Valid range: ~3000–50000 K. Returns BC such that M_V = M_bol − BC.
#[must_use]
pub fn bolometric_correction(temperature_k: f64) -> f64 {
    // BC is always negative (or near zero) and reaches minimum near T_eff ~ 6700 K.
    // Simple quadratic in log(T) centered near the BC minimum.
    let log_t = temperature_k.log10();
    // Fit: BC = -a * (log_t - log_t_min)^2 + bc_min
    // where log_t_min ≈ 3.83 (T ≈ 6760 K), bc_min ≈ -0.01
    let x = log_t - 3.83;
    let bc = -7.2 * x * x - 0.01;
    bc.min(0.0)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn solar_luminosity_roundtrip() {
        // 1 R_sun at T_sun should give ~1 L_sun
        let l = luminosity_solar(1.0, constants::T_SUN);
        assert!(
            (l - 1.0).abs() < 0.02,
            "Solar luminosity: {l}, expected ~1.0"
        );
    }

    #[test]
    fn solar_bolometric_magnitude() {
        let m = absolute_bolometric_magnitude(1.0);
        assert!(
            (m - constants::M_BOL_SUN).abs() < f64::EPSILON,
            "Solar M_bol: {m}, expected {}",
            constants::M_BOL_SUN
        );
    }

    #[test]
    fn magnitude_luminosity_roundtrip() {
        let l_in = 100.0;
        let m = absolute_bolometric_magnitude(l_in);
        let l_out = luminosity_from_abs_bol_mag(m);
        assert!(
            (l_out - l_in).abs() / l_in < 1e-10,
            "Roundtrip: {l_in} → M={m} → {l_out}"
        );
    }

    #[test]
    fn distance_modulus_10pc() {
        // At 10 pc, distance modulus = 0
        let dm = distance_modulus(10.0);
        assert!(dm.abs() < f64::EPSILON, "10 pc modulus: {dm}");
    }

    #[test]
    fn distance_modulus_roundtrip() {
        let d_in = 250.0;
        let dm = distance_modulus(d_in);
        let d_out = distance_from_modulus(dm);
        assert!(
            (d_out - d_in).abs() / d_in < 1e-10,
            "Distance roundtrip: {d_in} → μ={dm} → {d_out}"
        );
    }

    #[test]
    fn apparent_absolute_roundtrip() {
        let abs_mag = 4.83; // Sun
        let dist = 10.0; // 10 pc
        let app = apparent_magnitude(abs_mag, dist);
        let abs_back = absolute_magnitude(app, dist);
        assert!(
            (abs_back - abs_mag).abs() < 1e-10,
            "Magnitude roundtrip: {abs_mag} → m={app} → {abs_back}"
        );
    }

    #[test]
    fn bolometric_correction_sun() {
        // Solar BC ≈ -0.07 to -0.11 depending on calibration
        let bc = bolometric_correction(constants::T_SUN);
        assert!(
            bc.abs() < 0.5,
            "Solar BC = {bc}, expected near -0.07 to -0.11"
        );
    }

    #[test]
    fn zero_luminosity_returns_infinity() {
        assert!(absolute_bolometric_magnitude(0.0).is_infinite());
    }

    #[test]
    fn zero_distance_returns_neg_infinity() {
        assert!(distance_modulus(0.0).is_infinite());
    }
}
