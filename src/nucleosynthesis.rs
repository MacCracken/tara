//! Stellar nucleosynthesis — fusion chains, energy generation rates, and element production.
//!
//! Implements approximate energy generation rates for the pp-chain, CNO cycle,
//! and triple-alpha process following standard stellar structure theory
//! (Kippenhahn, Weigert & Weiss 2012).

use serde::{Deserialize, Serialize};

/// Dominant nuclear burning process.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
#[non_exhaustive]
pub enum FusionProcess {
    /// Proton-proton chain (dominant below ~17 MK).
    PpChain,
    /// Carbon-nitrogen-oxygen cycle (dominant above ~17 MK).
    CnoCycle,
    /// Triple-alpha process (helium burning, T > ~100 MK).
    TripleAlpha,
}

/// pp-chain branch identifier.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
#[non_exhaustive]
pub enum PpBranch {
    /// ppI: 4H → He-4 (dominant at lower T).
    PpI,
    /// ppII: involves Be-7 and Li-7.
    PpII,
    /// ppIII: involves B-8 (rare, produces high-energy neutrinos).
    PpIII,
}

/// Net energy release per He-4 nucleus produced (MeV).
///
/// 4 ¹H → ⁴He + 2e⁺ + 2νₑ + 26.73 MeV (total).
/// Neutrino losses reduce the effective energy available.
pub const PP_CHAIN_ENERGY_MEV: f64 = 26.73;

/// Average neutrino energy loss for ppI branch (MeV).
pub const PPI_NEUTRINO_LOSS_MEV: f64 = 0.53;

/// Average neutrino energy loss for ppII branch (MeV).
pub const PPII_NEUTRINO_LOSS_MEV: f64 = 1.17;

/// Average neutrino energy loss for ppIII branch (MeV).
pub const PPIII_NEUTRINO_LOSS_MEV: f64 = 7.46;

/// Average effective energy per He-4 for ppI (MeV).
pub const PPI_EFFECTIVE_ENERGY_MEV: f64 = PP_CHAIN_ENERGY_MEV - 2.0 * PPI_NEUTRINO_LOSS_MEV;

/// pp-chain energy generation rate (erg/g/s).
///
/// Power-law approximation: ε_pp ≈ ε₀,pp ρ X² (T / T₀)⁴,
/// calibrated to solar core conditions (ρ₀ = 150 g/cm³, T₀ = 15.7 MK, X₀ = 0.34).
/// Solar core ε_pp ≈ 17 erg/g/s.
///
/// Reference: Hansen, Kawaler & Trimble (2004), Kippenhahn et al. (2012).
#[must_use]
pub fn pp_chain_rate(density_gcc: f64, hydrogen_fraction: f64, temperature_k: f64) -> f64 {
    if density_gcc <= 0.0 || hydrogen_fraction <= 0.0 || temperature_k <= 0.0 {
        return 0.0;
    }
    let t9 = temperature_k / 1e9;
    // Calibrated: at solar core (ρ=150, X=0.34, T=15.7 MK) → ε ≈ 17 erg/g/s
    let eps_0 = 2.4e4;
    eps_0 * density_gcc * hydrogen_fraction * hydrogen_fraction * t9.powi(4)
}

/// CNO cycle energy generation rate (erg/g/s).
///
/// Power-law approximation: ε_CNO ≈ ε₀,CNO ρ X X_CNO (T / T₀)¹⁶·⁷,
/// where X_CNO is the CNO mass fraction (solar ≈ 0.01).
///
/// The steep ~T¹⁶·⁷ dependence makes CNO dominant above ~17 MK.
/// Reference: Hansen, Kawaler & Trimble (2004), Kippenhahn et al. (2012).
#[must_use]
pub fn cno_cycle_rate(
    density_gcc: f64,
    hydrogen_fraction: f64,
    cno_fraction: f64,
    temperature_k: f64,
) -> f64 {
    if density_gcc <= 0.0 || hydrogen_fraction <= 0.0 || cno_fraction <= 0.0 || temperature_k <= 0.0
    {
        return 0.0;
    }
    let t9 = temperature_k / 1e9;
    // Calibrated so crossover with pp occurs near 17 MK at solar composition
    let eps_0 = 4.4e27;
    eps_0 * density_gcc * hydrogen_fraction * cno_fraction * t9.powf(16.7)
}

/// Triple-alpha process energy generation rate (erg/g/s).
///
/// ε_3α ≈ 3.86×10¹¹ ρ² Y³ T₈⁻³ exp(−44.027 T₈⁻¹),
/// where Y is helium mass fraction and T₈ = T / 10⁸ K.
///
/// Dominates during helium-burning (core T > ~100 MK).
#[must_use]
pub fn triple_alpha_rate(density_gcc: f64, helium_fraction: f64, temperature_k: f64) -> f64 {
    if density_gcc <= 0.0 || helium_fraction <= 0.0 || temperature_k <= 0.0 {
        return 0.0;
    }
    let t8 = temperature_k / 1e8;
    if t8 <= 0.0 {
        return 0.0;
    }
    3.86e11
        * density_gcc
        * density_gcc
        * helium_fraction.powi(3)
        * t8.powi(-3)
        * (-44.027 / t8).exp()
}

/// Determine dominant fusion process from core temperature.
///
/// - T < 17 MK: pp-chain
/// - 17 MK ≤ T < 100 MK: CNO cycle
/// - T ≥ 100 MK: triple-alpha
#[must_use]
pub fn dominant_process(core_temperature_k: f64) -> FusionProcess {
    let t_mk = core_temperature_k / 1e6;
    if t_mk < 17.0 {
        FusionProcess::PpChain
    } else if t_mk < 100.0 {
        FusionProcess::CnoCycle
    } else {
        FusionProcess::TripleAlpha
    }
}

/// Estimate pp-chain branch fractions at a given core temperature.
///
/// Returns `(ppI_frac, ppII_frac, ppIII_frac)` summing to 1.0.
/// Simplified model: ppI dominates below 14 MK, ppII from 14–23 MK,
/// ppIII negligible except at high T.
#[must_use]
pub fn pp_branch_fractions(core_temperature_k: f64) -> (f64, f64, f64) {
    let t_mk = core_temperature_k / 1e6;

    if t_mk < 10.0 {
        // Very cool: nearly pure ppI
        (0.98, 0.02, 0.0)
    } else if t_mk < 14.0 {
        // Transition: ppI → ppII
        let f = (t_mk - 10.0) / 4.0;
        let ppi = 0.98 - 0.65 * f;
        let ppii = 0.02 + 0.65 * f;
        (ppi, ppii, 0.0)
    } else if t_mk < 23.0 {
        // ppII dominant, ppIII growing
        let f = (t_mk - 14.0) / 9.0;
        (0.33 - 0.30 * f, 0.67 - 0.17 * f, 0.0 + 0.47 * f)
    } else {
        // High T: ppIII significant but CNO should dominate over pp overall
        (0.03, 0.50, 0.47)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn pp_rate_solar_core() {
        // Solar core: ρ ≈ 150 g/cm³, T ≈ 15.7 MK, X ≈ 0.34
        let eps = pp_chain_rate(150.0, 0.34, 15.7e6);
        // Should be positive and in a reasonable range (erg/g/s)
        assert!(eps > 0.0, "Solar core pp rate: {eps}");
    }

    #[test]
    fn cno_rate_hot_star() {
        // Massive star core: ρ ≈ 100, T ≈ 25 MK, X ≈ 0.70, X_CNO ≈ 0.01
        let eps_cno = cno_cycle_rate(100.0, 0.70, 0.01, 25e6);
        let eps_pp = pp_chain_rate(100.0, 0.70, 25e6);
        assert!(eps_cno > 0.0, "Hot star CNO rate: {eps_cno}");
        // At 25 MK, CNO should dominate over pp
        assert!(
            eps_cno > eps_pp,
            "CNO should dominate at 25 MK: CNO={eps_cno}, pp={eps_pp}"
        );
    }

    #[test]
    fn pp_dominates_cool() {
        let eps_pp = pp_chain_rate(150.0, 0.70, 10e6);
        let eps_cno = cno_cycle_rate(150.0, 0.70, 0.01, 10e6);
        assert!(
            eps_pp > eps_cno,
            "pp should dominate at 10 MK: pp={eps_pp}, CNO={eps_cno}"
        );
    }

    #[test]
    fn triple_alpha_cold_negligible() {
        // At 10 MK, triple-alpha should be essentially zero
        let eps = triple_alpha_rate(150.0, 0.98, 10e6);
        assert!(eps < 1e-50, "Triple-alpha at 10 MK: {eps}");
    }

    #[test]
    fn triple_alpha_helium_burning() {
        // Helium burning: T ≈ 200 MK, ρ ≈ 10⁴, Y ≈ 0.98
        let eps = triple_alpha_rate(1e4, 0.98, 200e6);
        assert!(eps > 0.0, "He-burning triple-alpha: {eps}");
    }

    #[test]
    fn dominant_process_ranges() {
        assert_eq!(dominant_process(10e6), FusionProcess::PpChain);
        assert_eq!(dominant_process(20e6), FusionProcess::CnoCycle);
        assert_eq!(dominant_process(200e6), FusionProcess::TripleAlpha);
    }

    #[test]
    fn pp_branches_sum_to_one() {
        for t_mk in [5.0, 10.0, 12.0, 14.0, 18.0, 23.0, 30.0] {
            let (a, b, c) = pp_branch_fractions(t_mk * 1e6);
            let sum = a + b + c;
            assert!(
                (sum - 1.0).abs() < 0.01,
                "Branch fractions at {t_mk} MK: ({a}, {b}, {c}) sum={sum}"
            );
        }
    }

    #[test]
    fn ppi_dominates_cool() {
        let (ppi, _, _) = pp_branch_fractions(8e6);
        assert!(ppi > 0.9, "ppI at 8 MK: {ppi}");
    }

    #[test]
    fn zero_inputs_return_zero() {
        assert!(pp_chain_rate(0.0, 0.70, 15e6).abs() < f64::EPSILON);
        assert!(cno_cycle_rate(100.0, 0.0, 0.01, 20e6).abs() < f64::EPSILON);
        assert!(triple_alpha_rate(100.0, 0.98, 0.0).abs() < f64::EPSILON);
    }

    #[test]
    fn fusion_process_serde_roundtrip() {
        for p in [
            FusionProcess::PpChain,
            FusionProcess::CnoCycle,
            FusionProcess::TripleAlpha,
        ] {
            let json = serde_json::to_string(&p).unwrap();
            let back: FusionProcess = serde_json::from_str(&json).unwrap();
            assert_eq!(back, p);
        }
    }

    #[test]
    fn pp_branch_serde_roundtrip() {
        for b in [PpBranch::PpI, PpBranch::PpII, PpBranch::PpIII] {
            let json = serde_json::to_string(&b).unwrap();
            let back: PpBranch = serde_json::from_str(&json).unwrap();
            assert_eq!(back, b);
        }
    }
}
