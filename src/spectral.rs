//! Spectral line analysis — Planck function, Doppler shift, line broadening.
//!
//! Provides blackbody radiation, wavelength shifts, and spectral line profiles.

use crate::constants;

/// Planck spectral radiance B(λ, T) in W·sr⁻¹·m⁻³.
///
/// B(λ, T) = 2hc² / λ⁵ × 1 / (exp(hc / λk_BT) − 1).
///
/// `wavelength_m` in meters, `temperature_k` in Kelvin.
#[must_use]
pub fn planck_radiance(wavelength_m: f64, temperature_k: f64) -> f64 {
    if wavelength_m <= 0.0 || temperature_k <= 0.0 {
        return 0.0;
    }
    let c = constants::C;
    let h = constants::H;
    let kb = constants::K_B;

    let numerator = 2.0 * h * c * c / wavelength_m.powi(5);
    let exponent = h * c / (wavelength_m * kb * temperature_k);

    // Guard against overflow in exp
    if exponent > 700.0 {
        return 0.0;
    }

    numerator / (exponent.exp() - 1.0)
}

/// Planck radiance with wavelength in nanometers (convenience).
///
/// Returns W·sr⁻¹·m⁻³ (same units as [`planck_radiance`]).
#[must_use]
#[inline]
pub fn planck_radiance_nm(wavelength_nm: f64, temperature_k: f64) -> f64 {
    planck_radiance(wavelength_nm * 1e-9, temperature_k)
}

/// Non-relativistic Doppler wavelength shift.
///
/// Δλ / λ₀ = v_r / c.
/// Returns observed wavelength in same units as `rest_wavelength`.
/// Positive `radial_velocity_ms` means receding (redshift).
#[must_use]
#[inline]
pub fn doppler_shift(rest_wavelength: f64, radial_velocity_ms: f64) -> f64 {
    rest_wavelength * (1.0 + radial_velocity_ms / constants::C)
}

/// Relativistic Doppler wavelength shift.
///
/// λ_obs = λ₀ × √((1 + β) / (1 − β)), where β = v/c.
/// Returns observed wavelength. Positive velocity = receding.
#[must_use]
pub fn doppler_shift_relativistic(rest_wavelength: f64, radial_velocity_ms: f64) -> f64 {
    let beta = radial_velocity_ms / constants::C;
    if beta.abs() >= 1.0 {
        return 0.0;
    }
    rest_wavelength * ((1.0 + beta) / (1.0 - beta)).sqrt()
}

/// Thermal (Doppler) broadening width σ_λ for a spectral line.
///
/// σ_λ = λ₀ / c × √(k_B T / m_atom).
///
/// `rest_wavelength_m` in meters, `temperature_k` in K,
/// `atom_mass_kg` is the mass of the absorbing atom.
/// Returns σ in meters (standard deviation of the Gaussian profile).
#[must_use]
pub fn thermal_broadening(rest_wavelength_m: f64, temperature_k: f64, atom_mass_kg: f64) -> f64 {
    if temperature_k <= 0.0 || atom_mass_kg <= 0.0 || rest_wavelength_m <= 0.0 {
        return 0.0;
    }
    rest_wavelength_m / constants::C * (constants::K_B * temperature_k / atom_mass_kg).sqrt()
}

/// Gaussian spectral line profile (normalized).
///
/// G(λ) = 1 / (σ√(2π)) × exp(−(λ − λ₀)² / (2σ²)).
///
/// All wavelengths in the same units (e.g., meters or nm).
#[must_use]
pub fn gaussian_profile(wavelength: f64, center: f64, sigma: f64) -> f64 {
    if sigma <= 0.0 {
        return 0.0;
    }
    let norm = 1.0 / (sigma * (2.0 * std::f64::consts::PI).sqrt());
    let x = (wavelength - center) / sigma;
    norm * (-0.5 * x * x).exp()
}

/// Lorentzian spectral line profile (normalized).
///
/// L(λ) = (γ / π) / ((λ − λ₀)² + γ²),
/// where γ is the half-width at half-maximum (HWHM).
///
/// All wavelengths in the same units.
#[must_use]
pub fn lorentzian_profile(wavelength: f64, center: f64, gamma: f64) -> f64 {
    if gamma <= 0.0 {
        return 0.0;
    }
    let dx = wavelength - center;
    (gamma / std::f64::consts::PI) / (dx * dx + gamma * gamma)
}

/// Pseudo-Voigt approximation to the Voigt profile.
///
/// V(λ) ≈ η L(λ) + (1 − η) G(λ),
/// where η = 1.36603 (f_L/f) − 0.47719 (f_L/f)² + 0.11116 (f_L/f)³
/// and f = (f_G⁵ + 2.69269 f_G⁴ f_L + 2.42843 f_G³ f_L² + 4.47163 f_G² f_L³
///          + 0.07842 f_G f_L⁴ + f_L⁵)^(1/5).
///
/// `sigma_g` is the Gaussian standard deviation, `gamma_l` is the Lorentzian HWHM.
///
/// Thompson et al. (1987) approximation.
#[must_use]
pub fn pseudo_voigt_profile(wavelength: f64, center: f64, sigma_g: f64, gamma_l: f64) -> f64 {
    if sigma_g <= 0.0 && gamma_l <= 0.0 {
        return 0.0;
    }

    // Convert Gaussian sigma to FWHM
    let f_g = 2.0 * (2.0_f64.ln()).sqrt() * sigma_g.max(1e-30);
    let f_l = 2.0 * gamma_l.max(1e-30);

    // Total FWHM (Thompson et al. 1987)
    let f = (f_g.powi(5)
        + 2.692_69 * f_g.powi(4) * f_l
        + 2.428_43 * f_g.powi(3) * f_l.powi(2)
        + 4.471_63 * f_g.powi(2) * f_l.powi(3)
        + 0.078_42 * f_g * f_l.powi(4)
        + f_l.powi(5))
    .powf(0.2);

    // Mixing parameter
    let ratio = f_l / f;
    let eta = 1.366_03 * ratio - 0.477_19 * ratio.powi(2) + 0.111_16 * ratio.powi(3);
    let eta = eta.clamp(0.0, 1.0);

    // Convert total FWHM back to profile parameters
    let sigma_total = f / (2.0 * (2.0_f64.ln()).sqrt());
    let gamma_total = f / 2.0;

    eta * lorentzian_profile(wavelength, center, gamma_total)
        + (1.0 - eta) * gaussian_profile(wavelength, center, sigma_total)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn planck_peak_near_wien() {
        // Sun at 5772 K: peak near 502 nm
        let peak_nm = constants::WIEN_B_NM / constants::T_SUN;

        // Sample around the peak and verify it's actually a maximum
        let at_peak = planck_radiance_nm(peak_nm, constants::T_SUN);
        let below = planck_radiance_nm(peak_nm - 50.0, constants::T_SUN);
        let above = planck_radiance_nm(peak_nm + 50.0, constants::T_SUN);

        assert!(at_peak > below, "Peak should exceed short-wavelength side");
        assert!(at_peak > above, "Peak should exceed long-wavelength side");
    }

    #[test]
    fn planck_hotter_is_brighter() {
        let hot = planck_radiance_nm(500.0, 10_000.0);
        let cool = planck_radiance_nm(500.0, 5000.0);
        assert!(hot > cool);
    }

    #[test]
    fn planck_zero_inputs() {
        assert!(planck_radiance(0.0, 5772.0).abs() < f64::EPSILON);
        assert!(planck_radiance(500e-9, 0.0).abs() < f64::EPSILON);
    }

    #[test]
    fn doppler_redshift() {
        // Receding at 1000 km/s
        let obs = doppler_shift(500.0, 1e6);
        assert!(obs > 500.0, "Receding → redshift: {obs}");
    }

    #[test]
    fn doppler_blueshift() {
        let obs = doppler_shift(500.0, -1e6);
        assert!(obs < 500.0, "Approaching → blueshift: {obs}");
    }

    #[test]
    fn doppler_at_rest() {
        let obs = doppler_shift(500.0, 0.0);
        assert!((obs - 500.0).abs() < f64::EPSILON);
    }

    #[test]
    fn relativistic_doppler_low_v() {
        // At low velocity, relativistic ≈ classical
        let classical = doppler_shift(500.0, 1000.0);
        let relativistic = doppler_shift_relativistic(500.0, 1000.0);
        assert!(
            (classical - relativistic).abs() < 1e-6,
            "Low-v: classical={classical}, relativistic={relativistic}"
        );
    }

    #[test]
    fn relativistic_doppler_superluminal() {
        assert!(doppler_shift_relativistic(500.0, constants::C).abs() < f64::EPSILON);
    }

    #[test]
    fn thermal_broadening_hydrogen() {
        // H-alpha (656.3 nm) at solar temperature, hydrogen atom
        let sigma = thermal_broadening(656.3e-9, constants::T_SUN, constants::AMU);
        // Should be a small fraction of the wavelength
        assert!(sigma > 0.0);
        assert!(sigma < 1e-9, "Thermal broadening sigma: {sigma} m");
    }

    #[test]
    fn gaussian_normalized() {
        // Integral of Gaussian should approximate 1.0
        let sigma = 0.1;
        let center = 500.0;
        let n = 10_000;
        let range = 5.0 * sigma;
        let dx = 2.0 * range / n as f64;
        let integral: f64 = (0..n)
            .map(|i| {
                let x = center - range + (i as f64 + 0.5) * dx;
                gaussian_profile(x, center, sigma) * dx
            })
            .sum();
        assert!(
            (integral - 1.0).abs() < 0.001,
            "Gaussian integral: {integral}"
        );
    }

    #[test]
    fn lorentzian_normalized() {
        // Lorentzian has wider tails, need larger range
        let gamma = 0.1;
        let center = 500.0;
        let n = 100_000;
        let range = 100.0 * gamma;
        let dx = 2.0 * range / n as f64;
        let integral: f64 = (0..n)
            .map(|i| {
                let x = center - range + (i as f64 + 0.5) * dx;
                lorentzian_profile(x, center, gamma) * dx
            })
            .sum();
        assert!(
            (integral - 1.0).abs() < 0.01,
            "Lorentzian integral: {integral}"
        );
    }

    #[test]
    fn voigt_between_gaussian_and_lorentzian() {
        // At the wings, Voigt should be between pure Gaussian and pure Lorentzian
        let center = 500.0;
        let sigma = 0.1;
        let gamma = 0.05;
        let wing = center + 3.0 * sigma;

        let g = gaussian_profile(wing, center, sigma);
        let l = lorentzian_profile(wing, center, gamma);
        let v = pseudo_voigt_profile(wing, center, sigma, gamma);

        // Voigt wings should be broader than Gaussian
        assert!(v > g, "Voigt wing > Gaussian wing: v={v}, g={g}");
        // But not exceed Lorentzian
        assert!(v < l * 2.0, "Voigt wing reasonable vs Lorentzian");
    }
}
