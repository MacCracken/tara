//! Cross-crate bridges — convert primitive values from other AGNOS science crates
//! into tara stellar astrophysics parameters and vice versa.

use crate::constants;

// ── Falak bridges (orbital mechanics) ──────────────────────────────────────

/// Convert stellar mass (solar masses) to gravitational parameter μ (m³/s²).
///
/// Uses the heliocentric gravitational constant GM_sun directly
/// (IAU 2015), which is known to higher precision than G × M_sun.
#[must_use]
#[inline]
pub fn stellar_mass_solar_to_mu(mass_solar: f64) -> f64 {
    mass_solar * constants::GM_SUN
}

/// Convert luminosity class to approximate absolute magnitude.
///
/// Rough mapping for typical stars of each luminosity class.
#[must_use]
pub fn luminosity_class_to_abs_magnitude(
    luminosity_class: crate::classification::LuminosityClass,
) -> f64 {
    use crate::classification::LuminosityClass;
    match luminosity_class {
        LuminosityClass::Ia => -7.0,  // luminous supergiant
        LuminosityClass::Ib => -5.0,  // less luminous supergiant
        LuminosityClass::II => -3.0,  // bright giant
        LuminosityClass::III => 0.0,  // normal giant
        LuminosityClass::IV => 2.0,   // subgiant
        LuminosityClass::V => 5.0,    // main sequence (solar-like)
        LuminosityClass::VI => 10.0,  // subdwarf
        LuminosityClass::VII => 15.0, // white dwarf
    }
}

// ── Prakash bridges (optics) ───────────────────────────────────────────────

/// Convert effective temperature (K) to peak blackbody wavelength (nm).
///
/// Wien's displacement law: λ_max = b / T, using CODATA 2018 Wien constant.
#[must_use]
#[inline]
pub fn temperature_to_peak_wavelength_nm(temperature_k: f64) -> f64 {
    if temperature_k <= 0.0 {
        return 0.0;
    }
    constants::WIEN_B_NM / temperature_k
}

/// Convert surface gravity log(g) to limb darkening coefficient (linear).
///
/// Simplified approximation: u ≈ 0.6 - 0.1 × (log_g - 4.0).
/// Solar log_g ≈ 4.44. This is a rough bridge value; real limb darkening
/// depends on wavelength, temperature, and metallicity (see Claret 2000).
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
        // Should be exactly GM_SUN for 1 solar mass
        assert!(
            (mu - constants::GM_SUN).abs() < f64::EPSILON,
            "1 M_sun → μ = {mu}, expected {}",
            constants::GM_SUN
        );
    }

    #[test]
    fn luminosity_main_sequence() {
        use crate::classification::LuminosityClass;
        let m = luminosity_class_to_abs_magnitude(LuminosityClass::V);
        assert!((m - 5.0).abs() < 0.1);
    }

    #[test]
    fn luminosity_supergiant_ia_brighter_than_ib() {
        use crate::classification::LuminosityClass;
        let ia = luminosity_class_to_abs_magnitude(LuminosityClass::Ia);
        let ib = luminosity_class_to_abs_magnitude(LuminosityClass::Ib);
        assert!(ia < ib, "Ia should be brighter (lower mag) than Ib");
    }

    #[test]
    fn wien_sun() {
        // Sun (IAU 2015): T = 5772 K → λ_max ≈ 502 nm
        let nm = temperature_to_peak_wavelength_nm(constants::T_SUN);
        assert!(
            (nm - 502.0).abs() < 5.0,
            "Solar peak: {nm} nm, expected ~502"
        );
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
