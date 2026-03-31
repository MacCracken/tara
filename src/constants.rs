//! Physical and astronomical constants for stellar astrophysics.
//!
//! Values follow IAU 2015 Resolution B3 (nominal solar values) and
//! CODATA 2018 (physical constants) unless otherwise noted.

// ── IAU 2015 nominal solar values ─────────────────────────────────────────

/// Nominal solar luminosity (W). IAU 2015 Resolution B3.
pub const L_SUN: f64 = 3.828e26;

/// Nominal solar effective temperature (K). IAU 2015 Resolution B3.
pub const T_SUN: f64 = 5772.0;

/// Nominal solar radius (m). IAU 2015 Resolution B3.
pub const R_SUN: f64 = 6.957e8;

/// Solar absolute bolometric magnitude. IAU 2015 Resolution B2.
pub const M_BOL_SUN: f64 = 4.74;

/// Heliocentric gravitational constant GM_sun (m³/s²).
///
/// Known to much higher precision than G or M_sun individually.
/// TCB-compatible value from IAU 2015 Resolution B3.
pub const GM_SUN: f64 = 1.327_124_4e20;

/// Nominal solar mass (kg), derived from GM_sun / G.
///
/// Prefer [`GM_SUN`] directly for gravitational calculations.
pub const M_SUN: f64 = 1.988_409_87e30;

// ── CODATA 2018 physical constants ────────────────────────────────────────

/// Newtonian gravitational constant (m³ kg⁻¹ s⁻²). CODATA 2018.
pub const G: f64 = 6.674_30e-11;

/// Stefan-Boltzmann constant (W m⁻² K⁻⁴). CODATA 2018.
pub const SIGMA_SB: f64 = 5.670_374_419e-8;

/// Wien displacement constant (m·K). CODATA 2018.
pub const WIEN_B: f64 = 2.897_771_955e-3;

/// Wien displacement constant in nanometer·Kelvin.
pub const WIEN_B_NM: f64 = 2_897_771.955;

/// Boltzmann constant (J/K). CODATA 2018 (exact by definition).
pub const K_B: f64 = 1.380_649e-23;

/// Speed of light in vacuum (m/s). Exact by definition.
pub const C: f64 = 2.997_924_58e8;

/// Planck constant (J·s). CODATA 2018 (exact by definition).
pub const H: f64 = 6.626_070_15e-34;

/// Atomic mass unit (kg). CODATA 2018.
pub const AMU: f64 = 1.660_539_066_6e-27;

// ── Spectral class temperature boundaries (K) ────────────────────────────
// Approximate boundaries for Morgan-Keenan classes (main sequence).
// Based on Pecaut & Mamajek (2013) and Habets & Heintze (1981).

/// Minimum effective temperature for O-type stars (K).
pub const T_O_MIN: f64 = 30_000.0;
/// Minimum effective temperature for B-type stars (K).
pub const T_B_MIN: f64 = 10_000.0;
/// Minimum effective temperature for A-type stars (K).
pub const T_A_MIN: f64 = 7_500.0;
/// Minimum effective temperature for F-type stars (K).
pub const T_F_MIN: f64 = 6_000.0;
/// Minimum effective temperature for G-type stars (K).
pub const T_G_MIN: f64 = 5_200.0;
/// Minimum effective temperature for K-type stars (K).
pub const T_K_MIN: f64 = 3_700.0;
// M-type: below T_K_MIN

/// Minimum effective temperature for M-type stars (K).
/// L-type begins below this boundary.
pub const T_M_MIN: f64 = 2_100.0;
/// Minimum effective temperature for L-type (brown dwarfs, K).
pub const T_L_MIN: f64 = 1_300.0;
/// Minimum effective temperature for T-type (methane dwarfs, K).
pub const T_T_MIN: f64 = 500.0;
// Y-type: below T_T_MIN

/// Minimum effective temperature for Wolf-Rayet (WR) stars (K).
/// WR stars are hotter than O-type.
pub const T_W_MIN: f64 = 50_000.0;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn gm_sun_consistent() {
        // GM_sun should approximately equal G * M_sun
        let computed = G * M_SUN;
        let rel_err = (computed - GM_SUN).abs() / GM_SUN;
        assert!(rel_err < 1e-6, "GM_sun consistency: rel_err = {rel_err}");
    }

    #[test]
    fn stefan_boltzmann_solar_luminosity() {
        // L_sun = 4π R_sun² σ T_sun⁴
        let l = 4.0 * std::f64::consts::PI * R_SUN * R_SUN * SIGMA_SB * T_SUN.powi(4);
        let rel_err = (l - L_SUN).abs() / L_SUN;
        assert!(
            rel_err < 0.01,
            "Stefan-Boltzmann solar check: L = {l:.3e}, expected {L_SUN:.3e}, rel_err = {rel_err}"
        );
    }

    #[test]
    fn wien_solar_peak() {
        // Sun peak wavelength ~ 502 nm (IAU T_sun = 5772 K)
        let peak_nm = WIEN_B_NM / T_SUN;
        assert!(
            (peak_nm - 502.0).abs() < 5.0,
            "Solar peak wavelength: {peak_nm} nm"
        );
    }

    #[test]
    fn spectral_boundaries_ordered() {
        const {
            assert!(T_O_MIN > T_B_MIN);
            assert!(T_B_MIN > T_A_MIN);
            assert!(T_A_MIN > T_F_MIN);
            assert!(T_F_MIN > T_G_MIN);
            assert!(T_G_MIN > T_K_MIN);
        }
    }
}
