//! Stellar evolution — main sequence lifetime, mass-luminosity relation,
//! evolutionary phases, and remnant type prediction.
//!
//! Uses piecewise mass-luminosity relation from Duric (2004) and
//! standard remnant mass thresholds.

use serde::{Deserialize, Serialize};

/// Stellar evolutionary phase.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
#[non_exhaustive]
pub enum EvolutionaryPhase {
    /// Pre-main-sequence (contracting toward ZAMS).
    PreMainSequence,
    /// Zero-age main sequence.
    ZeroAgeMainSequence,
    /// Core hydrogen burning.
    MainSequence,
    /// Subgiant branch (hydrogen shell burning).
    Subgiant,
    /// Red giant branch (deep convective envelope).
    RedGiant,
    /// Horizontal branch / red clump (core helium burning).
    HorizontalBranch,
    /// Asymptotic giant branch (double shell burning).
    AsymptoticGiantBranch,
    /// Post-AGB / planetary nebula phase.
    PostAgb,
    /// White dwarf remnant.
    WhiteDwarf,
    /// Neutron star remnant.
    NeutronStar,
    /// Black hole remnant.
    BlackHole,
}

/// Stellar remnant type based on initial mass.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
#[non_exhaustive]
pub enum RemnantType {
    /// White dwarf (initial mass < ~8 M_sun).
    WhiteDwarf,
    /// Neutron star (initial mass ~8–25 M_sun).
    NeutronStar,
    /// Black hole (initial mass > ~25 M_sun).
    BlackHole,
}

/// Mass-luminosity relation for main-sequence stars.
///
/// Piecewise power law (Duric 2004):
/// - M < 0.43 M_sun: L/L_sun = 0.23 (M/M_sun)^2.3
/// - 0.43 ≤ M < 2 M_sun: L/L_sun = (M/M_sun)^4.0
/// - 2 ≤ M < 55 M_sun: L/L_sun = 1.4 (M/M_sun)^3.5
/// - M ≥ 55 M_sun: L/L_sun = 32000 (M/M_sun)
///
/// Returns luminosity in solar units.
#[must_use]
pub fn mass_luminosity(mass_solar: f64) -> f64 {
    if mass_solar <= 0.0 {
        tracing::warn!(
            mass_solar,
            "non-positive mass in mass_luminosity — returning 0"
        );
        return 0.0;
    }
    if mass_solar < 0.43 {
        0.23 * mass_solar.powf(2.3)
    } else if mass_solar < 2.0 {
        mass_solar.powf(4.0)
    } else if mass_solar < 55.0 {
        1.4 * mass_solar.powf(3.5)
    } else {
        32_000.0 * mass_solar
    }
}

/// Main-sequence lifetime estimate (years).
///
/// t_MS ≈ t_sun × (M/M_sun) / (L/L_sun), using the mass-luminosity relation.
/// Solar main-sequence lifetime ≈ 10 Gyr.
///
/// Returns lifetime in years.
#[must_use]
pub fn main_sequence_lifetime(mass_solar: f64) -> f64 {
    if mass_solar <= 0.0 {
        return 0.0;
    }
    let l = mass_luminosity(mass_solar);
    if l <= 0.0 {
        return 0.0;
    }
    1e10 * mass_solar / l
}

/// Predict remnant type from initial stellar mass (solar masses).
///
/// Thresholds based on standard stellar evolution theory:
/// - < 8 M_sun → white dwarf
/// - 8–25 M_sun → neutron star (core-collapse supernova)
/// - > 25 M_sun → black hole (direct collapse or fallback)
#[must_use]
pub fn remnant_type(initial_mass_solar: f64) -> RemnantType {
    if initial_mass_solar < 8.0 {
        RemnantType::WhiteDwarf
    } else if initial_mass_solar < 25.0 {
        RemnantType::NeutronStar
    } else {
        RemnantType::BlackHole
    }
}

/// Chandrasekhar mass limit for white dwarfs (solar masses).
pub const CHANDRASEKHAR_LIMIT: f64 = 1.44;

/// White dwarf final mass from initial mass (solar masses).
///
/// Initial-final mass relation from Kalirai et al. (2008):
/// M_final = 0.109 M_initial + 0.394 M_sun.
///
/// Valid for initial mass ~1–8 M_sun. Returns final mass in solar masses,
/// clamped to the Chandrasekhar limit (~1.44 M_sun).
#[must_use]
#[inline]
pub fn white_dwarf_mass(initial_mass_solar: f64) -> f64 {
    (0.109 * initial_mass_solar + 0.394).clamp(0.0, CHANDRASEKHAR_LIMIT)
}

/// Determine approximate evolutionary phase from mass (solar masses) and age (years).
///
/// Simplified model: compares age against main-sequence lifetime to estimate phase.
/// For detailed evolution, use the Hurley/SSE fitting formulae.
#[must_use]
pub fn evolutionary_phase(mass_solar: f64, age_years: f64) -> EvolutionaryPhase {
    if mass_solar <= 0.0 || age_years < 0.0 {
        tracing::warn!(
            mass_solar,
            age_years,
            "non-physical input to evolutionary_phase — defaulting to PreMainSequence"
        );
        return EvolutionaryPhase::PreMainSequence;
    }

    let t_ms = main_sequence_lifetime(mass_solar);

    // Post-MS phase boundaries are approximate fractions of t_MS.
    // These will be replaced by SSE-based phase determination once
    // post-MS fitting formulae are implemented.
    const F_ZAMS: f64 = 0.001; // ZAMS → MS transition
    const F_SUBGIANT: f64 = 1.1; // MS → subgiant
    const F_RGB: f64 = 1.3; // subgiant → RGB
    const F_HB: f64 = 1.4; // RGB → horizontal branch
    const F_AGB: f64 = 1.5; // HB → AGB
    const F_POST_AGB: f64 = 1.6; // AGB → post-AGB (low-mass only)

    if age_years < 1e6 {
        // Very young — still contracting
        EvolutionaryPhase::PreMainSequence
    } else if age_years < t_ms * F_ZAMS {
        EvolutionaryPhase::ZeroAgeMainSequence
    } else if age_years < t_ms {
        EvolutionaryPhase::MainSequence
    } else if age_years < t_ms * F_SUBGIANT {
        EvolutionaryPhase::Subgiant
    } else if age_years < t_ms * F_RGB {
        EvolutionaryPhase::RedGiant
    } else if age_years < t_ms * F_HB {
        EvolutionaryPhase::HorizontalBranch
    } else if age_years < t_ms * F_AGB {
        EvolutionaryPhase::AsymptoticGiantBranch
    } else if age_years < t_ms * F_POST_AGB && mass_solar < 8.0 {
        EvolutionaryPhase::PostAgb
    } else if mass_solar < 8.0 {
        EvolutionaryPhase::WhiteDwarf
    } else if mass_solar < 25.0 {
        EvolutionaryPhase::NeutronStar
    } else {
        EvolutionaryPhase::BlackHole
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn solar_luminosity_from_mass() {
        let l = mass_luminosity(1.0);
        assert!(
            (l - 1.0).abs() < 0.01,
            "Solar mass-luminosity: {l}, expected 1.0"
        );
    }

    #[test]
    fn massive_star_brighter() {
        assert!(mass_luminosity(10.0) > mass_luminosity(1.0));
        assert!(mass_luminosity(50.0) > mass_luminosity(10.0));
    }

    #[test]
    fn low_mass_dimmer() {
        assert!(mass_luminosity(0.1) < mass_luminosity(1.0));
    }

    #[test]
    fn solar_lifetime() {
        let t = main_sequence_lifetime(1.0);
        // Solar MS lifetime ≈ 10 Gyr
        assert!((t - 1e10).abs() / 1e10 < 0.01, "Solar lifetime: {t:.2e} yr");
    }

    #[test]
    fn massive_stars_shorter_lived() {
        assert!(main_sequence_lifetime(10.0) < main_sequence_lifetime(1.0));
    }

    #[test]
    fn low_mass_longer_lived() {
        assert!(main_sequence_lifetime(0.1) > main_sequence_lifetime(1.0));
    }

    #[test]
    fn remnant_types() {
        assert_eq!(remnant_type(1.0), RemnantType::WhiteDwarf);
        assert_eq!(remnant_type(5.0), RemnantType::WhiteDwarf);
        assert_eq!(remnant_type(10.0), RemnantType::NeutronStar);
        assert_eq!(remnant_type(20.0), RemnantType::NeutronStar);
        assert_eq!(remnant_type(30.0), RemnantType::BlackHole);
    }

    #[test]
    fn solar_evolution_main_sequence() {
        let phase = evolutionary_phase(1.0, 4.6e9);
        assert_eq!(phase, EvolutionaryPhase::MainSequence);
    }

    #[test]
    fn old_solar_mass_becomes_wd() {
        // Well past PostAgb (> 1.6× t_ms ≈ 16 Gyr)
        let phase = evolutionary_phase(1.0, 20e9);
        assert_eq!(phase, EvolutionaryPhase::WhiteDwarf);
    }

    #[test]
    fn evolutionary_phase_serde_roundtrip() {
        let phases = [
            EvolutionaryPhase::PreMainSequence,
            EvolutionaryPhase::MainSequence,
            EvolutionaryPhase::RedGiant,
            EvolutionaryPhase::WhiteDwarf,
            EvolutionaryPhase::NeutronStar,
            EvolutionaryPhase::BlackHole,
        ];
        for phase in &phases {
            let json = serde_json::to_string(phase).unwrap();
            let back: EvolutionaryPhase = serde_json::from_str(&json).unwrap();
            assert_eq!(&back, phase);
        }
    }

    #[test]
    fn remnant_type_serde_roundtrip() {
        for rt in [
            RemnantType::WhiteDwarf,
            RemnantType::NeutronStar,
            RemnantType::BlackHole,
        ] {
            let json = serde_json::to_string(&rt).unwrap();
            let back: RemnantType = serde_json::from_str(&json).unwrap();
            assert_eq!(back, rt);
        }
    }

    #[test]
    fn white_dwarf_mass_solar() {
        // 1 M_sun → WD mass ≈ 0.503 M_sun
        let m = white_dwarf_mass(1.0);
        assert!(
            (m - 0.503).abs() < 0.01,
            "Solar WD mass: {m}, expected ~0.503"
        );
    }

    #[test]
    fn white_dwarf_mass_clamped_to_chandrasekhar() {
        // Very massive progenitor should clamp to Chandrasekhar limit
        let m = white_dwarf_mass(100.0);
        assert!(
            (m - CHANDRASEKHAR_LIMIT).abs() < f64::EPSILON,
            "WD mass should clamp to Chandrasekhar limit: {m}, expected {CHANDRASEKHAR_LIMIT}"
        );
    }

    #[test]
    fn post_agb_phase() {
        // A solar-mass star just past AGB should be PostAgb before WD
        let t_ms = main_sequence_lifetime(1.0);
        let phase = evolutionary_phase(1.0, t_ms * 1.55);
        assert_eq!(phase, EvolutionaryPhase::PostAgb);
    }
}
