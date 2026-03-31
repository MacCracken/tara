//! Cross-crate bridges — convert primitive values from other AGNOS science crates
//! into tara stellar astrophysics parameters and vice versa.

// ── Falak bridges (orbital mechanics) ──────────────────────────────────────

/// Convert stellar mass (solar masses) to gravitational parameter μ (m³/s²).
///
/// μ = G × M, M_sun = 1.989e30 kg.
#[must_use]
#[inline]
pub fn stellar_mass_solar_to_mu(mass_solar: f64) -> f64 {
    const G: f64 = 6.674_30e-11;
    const M_SUN: f64 = 1.989e30;
    G * mass_solar * M_SUN
}

/// Convert luminosity class to approximate absolute magnitude.
///
/// Rough mapping for main-sequence stars.
#[must_use]
pub fn luminosity_class_to_abs_magnitude(luminosity_class: u8) -> f64 {
    match luminosity_class {
        1 => -6.0, // Ia supergiant
        2 => -3.0, // II bright giant
        3 => 0.0,  // III giant
        4 => 2.0,  // IV subgiant
        5 => 5.0,  // V main sequence (solar-like)
        6 => 10.0, // VI subdwarf
        7 => 15.0, // VII white dwarf
        _ => 5.0,  // default solar
    }
}

// ── Prakash bridges (optics) ───────────────────────────────────────────────

/// Convert effective temperature (K) to peak blackbody wavelength (nm).
///
/// Wien's displacement law: λ_max = 2.898e6 / T (nm).
#[must_use]
#[inline]
pub fn temperature_to_peak_wavelength_nm(temperature_k: f64) -> f64 {
    if temperature_k <= 0.0 {
        return 0.0;
    }
    2_898_000.0 / temperature_k
}

/// Convert surface gravity log(g) to limb darkening coefficient (linear).
///
/// Simplified: u ≈ 0.6 - 0.1 × (log_g - 4.0). Solar log_g ≈ 4.44.
#[must_use]
#[inline]
pub fn surface_gravity_to_limb_darkening(log_g: f64) -> f64 {
    (0.6 - 0.1 * (log_g - 4.0)).clamp(0.2, 0.9)
}

// ── Tanmatra bridges (atomic/nuclear physics) ──────────────────────────────

/// Convert stellar H/He mass fraction to nucleosynthesis yield scaling.
///
/// Higher hydrogen fraction → more fuel → longer main sequence.
/// Returns relative main-sequence lifetime scaling (1.0 for solar composition X=0.74).
#[must_use]
#[inline]
pub fn composition_to_lifetime_scale(hydrogen_fraction: f64) -> f64 {
    if hydrogen_fraction <= 0.0 {
        return 0.0;
    }
    (hydrogen_fraction / 0.74).clamp(0.0, 2.0)
}

/// Convert core temperature (K) to pp-chain vs CNO-cycle dominance factor.
///
/// Below ~17 MK: pp-chain dominates. Above: CNO dominates.
/// Returns CNO fraction (0.0 = pure pp, 1.0 = pure CNO).
#[must_use]
#[inline]
pub fn core_temp_to_cno_fraction(core_temperature_k: f64) -> f64 {
    let t_mk = core_temperature_k / 1e6;
    // Sigmoid transition around 17 MK
    let x = (t_mk - 17.0) / 3.0;
    1.0 / (1.0 + (-x).exp())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn sun_mu() {
        let mu = stellar_mass_solar_to_mu(1.0);
        assert!((mu - 1.327e20).abs() < 1e18);
    }

    #[test]
    fn luminosity_main_sequence() {
        let m = luminosity_class_to_abs_magnitude(5);
        assert!((m - 5.0).abs() < 0.1);
    }

    #[test]
    fn wien_sun() {
        // Sun: T ≈ 5778K → λ_max ≈ 501 nm
        let nm = temperature_to_peak_wavelength_nm(5778.0);
        assert!((nm - 501.0).abs() < 5.0);
    }

    #[test]
    fn wien_zero() {
        assert!(temperature_to_peak_wavelength_nm(0.0).abs() < f64::EPSILON);
    }

    #[test]
    fn limb_darkening_solar() {
        let u = surface_gravity_to_limb_darkening(4.44);
        assert!(u > 0.5 && u < 0.7);
    }

    #[test]
    fn composition_solar() {
        let s = composition_to_lifetime_scale(0.74);
        assert!((s - 1.0).abs() < 0.01);
    }

    #[test]
    fn cno_cool_star() {
        let f = core_temp_to_cno_fraction(10e6);
        assert!(f < 0.1, "cool star should be pp-dominated: {f}");
    }

    #[test]
    fn cno_hot_star() {
        let f = core_temp_to_cno_fraction(25e6);
        assert!(f > 0.9, "hot star should be CNO-dominated: {f}");
    }
}
