//! Stellar atmosphere models — opacity, radiative transfer, limb darkening.
//!
//! Provides fundamental atmosphere calculations using standard models.

use crate::constants;

/// Compute effective temperature (K) from luminosity (W) and radius (m).
///
/// Inverse Stefan-Boltzmann: T_eff = (L / (4π R² σ))^(1/4).
#[must_use]
#[inline]
pub fn effective_temperature(luminosity_w: f64, radius_m: f64) -> f64 {
    if luminosity_w <= 0.0 || radius_m <= 0.0 {
        return 0.0;
    }
    let flux = luminosity_w / (4.0 * std::f64::consts::PI * radius_m * radius_m);
    (flux / constants::SIGMA_SB).powf(0.25)
}

/// Compute surface gravity log₁₀(g) in CGS units (cm/s²).
///
/// log(g) = log₁₀(G M / R²), converted to CGS (multiply by 100).
/// Solar value: log(g) ≈ 4.44.
#[must_use]
#[inline]
pub fn log_surface_gravity(mass_solar: f64, radius_solar: f64) -> f64 {
    if mass_solar <= 0.0 || radius_solar <= 0.0 {
        return 0.0;
    }
    let m_kg = mass_solar * constants::M_SUN;
    let r_m = radius_solar * constants::R_SUN;
    let g_si = constants::G * m_kg / (r_m * r_m);
    // Convert m/s² to cm/s² (×100) then take log₁₀
    (g_si * 100.0).log10()
}

/// Grey atmosphere temperature profile (Eddington approximation).
///
/// T⁴(τ) = (3/4) T_eff⁴ (τ + 2/3), where τ is optical depth.
/// Returns temperature at optical depth τ.
#[must_use]
#[inline]
pub fn grey_atmosphere_temperature(t_eff: f64, optical_depth: f64) -> f64 {
    (0.75 * t_eff.powi(4) * (optical_depth + 2.0 / 3.0)).powf(0.25)
}

/// Kramers' opacity (bound-free + free-free), approximate.
///
/// κ ∝ ρ T^(−3.5). Returns opacity in cm²/g.
/// Normalization: κ_bf+ff ≈ 4.34×10²⁵ (1 + X) Z ρ T^(−3.5),
/// where X = hydrogen fraction, Z = metallicity.
#[must_use]
#[inline]
pub fn kramers_opacity(density_gcc: f64, temperature_k: f64, x: f64, z: f64) -> f64 {
    if temperature_k <= 0.0 {
        return 0.0;
    }
    4.34e25 * (1.0 + x) * z * density_gcc * temperature_k.powf(-3.5)
}

/// Electron scattering opacity (cm²/g).
///
/// κ_es = 0.2 (1 + X), where X is hydrogen mass fraction.
/// Dominant in hot stellar interiors.
#[must_use]
#[inline]
pub fn electron_scattering_opacity(hydrogen_fraction: f64) -> f64 {
    0.2 * (1.0 + hydrogen_fraction)
}

/// Pressure scale height (m).
///
/// H = k_B T / (μ m_H g), where μ is mean molecular weight
/// and g is surface gravity (m/s²).
#[must_use]
#[inline]
pub fn scale_height(temperature_k: f64, mean_molecular_weight: f64, gravity_ms2: f64) -> f64 {
    if mean_molecular_weight <= 0.0 || gravity_ms2 <= 0.0 {
        return 0.0;
    }
    constants::K_B * temperature_k / (mean_molecular_weight * constants::AMU * gravity_ms2)
}

/// Linear limb darkening: I(μ)/I(1) = 1 − u(1 − μ).
///
/// `mu` is cos(θ), the angle from disk center. `u` is the limb darkening coefficient.
#[must_use]
#[inline]
pub fn limb_darkening_linear(mu: f64, u: f64) -> f64 {
    (1.0 - u * (1.0 - mu)).clamp(0.0, 1.0)
}

/// Quadratic limb darkening: I(μ)/I(1) = 1 − a(1−μ) − b(1−μ)².
///
/// `mu` is cos(θ). `a` and `b` are the quadratic coefficients.
#[must_use]
#[inline]
pub fn limb_darkening_quadratic(mu: f64, a: f64, b: f64) -> f64 {
    let omu = 1.0 - mu;
    (1.0 - a * omu - b * omu * omu).clamp(0.0, 1.0)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn solar_effective_temperature() {
        let t = effective_temperature(constants::L_SUN, constants::R_SUN);
        assert!(
            (t - constants::T_SUN).abs() < 10.0,
            "Solar T_eff: {t} K, expected ~{}",
            constants::T_SUN
        );
    }

    #[test]
    fn solar_log_g() {
        let lg = log_surface_gravity(1.0, 1.0);
        assert!(
            (lg - 4.44).abs() < 0.02,
            "Solar log(g): {lg}, expected ~4.44"
        );
    }

    #[test]
    fn grey_atmosphere_surface() {
        // At τ = 2/3, T should equal T_eff
        let t = grey_atmosphere_temperature(5772.0, 2.0 / 3.0);
        let expected = (0.75 * 5772.0_f64.powi(4) * (2.0 / 3.0 + 2.0 / 3.0)).powf(0.25);
        assert!(
            (t - expected).abs() < 1.0,
            "Grey atmosphere at τ=2/3: {t}, expected {expected}"
        );
    }

    #[test]
    fn grey_atmosphere_deep() {
        // Deep in atmosphere, T should increase with τ
        let t_surface = grey_atmosphere_temperature(5772.0, 0.0);
        let t_deep = grey_atmosphere_temperature(5772.0, 10.0);
        assert!(t_deep > t_surface);
    }

    #[test]
    fn electron_scattering_solar() {
        // Solar X ≈ 0.74, κ_es ≈ 0.348 cm²/g
        let k = electron_scattering_opacity(0.74);
        assert!((k - 0.348).abs() < 0.01, "Solar electron scattering: {k}");
    }

    #[test]
    fn limb_darkening_center() {
        // At disk center (μ=1), intensity should be 1.0
        assert!((limb_darkening_linear(1.0, 0.6) - 1.0).abs() < f64::EPSILON);
        assert!((limb_darkening_quadratic(1.0, 0.4, 0.2) - 1.0).abs() < f64::EPSILON);
    }

    #[test]
    fn limb_darkening_edge() {
        // At limb (μ=0), intensity should be 1 - u (linear)
        let i = limb_darkening_linear(0.0, 0.6);
        assert!((i - 0.4).abs() < f64::EPSILON, "Limb at μ=0: {i}");
    }

    #[test]
    fn scale_height_positive() {
        // Solar photosphere: T ≈ 5772 K, μ ≈ 1.3, g ≈ 274 m/s²
        let h = scale_height(5772.0, 1.3, 274.0);
        // Expected ~150 km = 1.5e5 m
        assert!(h > 1e5 && h < 3e5, "Solar scale height: {h} m");
    }

    #[test]
    fn zero_inputs_return_zero() {
        assert!(effective_temperature(0.0, 1.0).abs() < f64::EPSILON);
        assert!(log_surface_gravity(0.0, 1.0).abs() < f64::EPSILON);
        assert!(kramers_opacity(1.0, 0.0, 0.74, 0.02).abs() < f64::EPSILON);
        assert!(scale_height(5772.0, 0.0, 274.0).abs() < f64::EPSILON);
    }
}
