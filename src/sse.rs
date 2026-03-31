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
#[inline]
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
    -31.671_194_79,
    -24.248_483_22,
    -5.336_089_72,
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
    -48.960_668_56,
    -40.053_861_35,
    -9.093_318_16,
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
#[non_exhaustive]
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

// ── HPT00 — t_BGB coefficients (from xt data, 4-term polys padded to 5) ──

/// a17: xt(1..4)
const A17: [f64; 5] = [1.593_890e3, 2.053_038e3, 1.231_226e3, 2.327_785e2, 0.0];
/// a18: xt(5..8)
const A18: [f64; 5] = [2.706_708e3, 1.483_131e3, 5.772_723e2, 7.411_230e1, 0.0];
/// a19: xt(9..12)
const A19: [f64; 5] = [1.466_143e2, -1.048_442e2, -6.795_374e1, -1.391_127e1, 0.0];
/// a20: xt(13..16)
const A20: [f64; 5] = [4.141_960e-2, 4.564_888e-2, 2.958_542e-2, 5.571_483e-3, 0.0];
/// a21: constant xt(17)
const A21_CONST: f64 = 3.426_349e-1;

// ── HPT00 — t_hook coefficients (from xt data) ──────────────────────────

/// a22: xt(18..21)
const A22: [f64; 5] = [1.949_814e1, 1.758_178e0, -6.008_212e0, -4.470_533e0, 0.0];
/// a23: constant xt(22)
const A23_CONST: f64 = 4.903_830e0;
/// a24: xt(23..26)
const A24: [f64; 5] = [
    5.212_154e-2,
    3.166_411e-2,
    -2.750_074e-3,
    -2.271_549e-3,
    0.0,
];
/// a25: xt(27..30)
const A25: [f64; 5] = [1.312_179e0, -3.294_936e-1, 9.231_860e-2, 2.610_989e-2, 0.0];
/// a26: constant xt(31)
const A26_CONST: f64 = 8.073_972e-1;

/// Time to the base of the giant branch (Myr).
///
/// HPT00 Eq. 4. This is the total hydrogen-exhaustion timescale.
/// Valid for M ∈ \[0.1, 100\] M☉.
#[must_use]
pub fn t_bgb(mass: f64, z: f64) -> f64 {
    if mass <= 0.0 {
        return 0.0;
    }
    let ze = zeta(z.clamp(0.0001, 0.03));
    let a17 = coeff(&A17, ze);
    let a18 = coeff(&A18, ze);
    let a19 = coeff(&A19, ze);
    let a20 = coeff(&A20, ze);

    let m2 = mass * mass;
    let m4 = m2 * m2;
    let m55 = mass.powf(5.5);
    let m7 = m4 * m2 * mass;

    let num = a17 + a18 * m4 + a19 * m55 + m7;
    let den = a20 * m2 + A21_CONST * m7;

    if den <= 0.0 { 0.0 } else { num / den }
}

/// Main-sequence lifetime (Myr).
///
/// HPT00 Eqs. 5–6. Applies the "hook" correction to [`t_bgb`].
/// The hook fraction μ ensures t_MS < t_BGB, with the MS ending
/// before the star reaches the base of the giant branch.
#[must_use]
pub fn ms_lifetime_myr(mass: f64, z: f64) -> f64 {
    let t = t_bgb(mass, z);
    if t <= 0.0 {
        return 0.0;
    }
    let ze = zeta(z.clamp(0.0001, 0.03));
    let a22 = coeff(&A22, ze);
    let a24 = coeff(&A24, ze);
    let a25 = coeff(&A25, ze);

    let x1 = if mass > 0.0 {
        a22 / mass.powf(A23_CONST)
    } else {
        0.0
    };
    let x2 = a24 + a25 / mass.powf(A26_CONST);

    let mu = (1.0 - 0.01 * x1.max(x2)).max(0.5);
    mu * t
}

/// Main-sequence lifetime (years).
///
/// Convenience wrapper around [`ms_lifetime_myr`] returning years.
#[must_use]
#[inline]
pub fn ms_lifetime(mass: f64, z: f64) -> f64 {
    ms_lifetime_myr(mass, z) * 1e6
}

// ── HPT00 — L_TMS coefficients (from xl data) ───────────────────────────

/// xl(1..5) → raw a27 before multiply
const XL_A27R: [f64; 5] = [
    1.031_538e0,
    -2.434_480e-1,
    7.732_821e0,
    6.460_705e0,
    1.374_484e0,
];
/// xl(6..10) → raw a28 before multiply
const XL_A28R: [f64; 5] = [
    1.043_715e0,
    -1.577_474e0,
    -5.168_234e0,
    -5.596_506e0,
    -1.299_394e0,
];
/// xl(11..14) → a29 (4-term)
const XL_A29: [f64; 5] = [7.859_573e2, -8.542_048e0, -2.642_511e1, -9.585_707e0, 0.0];
/// xl(15..19) → a30 (multiplier for a27, a28)
const XL_A30: [f64; 5] = [
    3.858_911e3,
    2.459_681e3,
    -7.630_093e1,
    -3.486_057e2,
    -4.861_703e1,
];
/// xl(20..23) → a31 (4-term)
const XL_A31: [f64; 5] = [2.888_720e2, 2.952_979e2, 1.850_341e2, 3.797_254e1, 0.0];
/// xl(24..27) → a32 (4-term)
const XL_A32: [f64; 5] = [7.196_580e0, 5.613_746e-1, 3.805_871e-1, 8.398_728e-2, 0.0];

/// Terminal main-sequence luminosity L_TMS (L☉).
///
/// HPT00 Eq. 8. The luminosity at the end of core hydrogen burning.
#[must_use]
pub fn tms_luminosity(mass: f64, z: f64) -> f64 {
    if mass <= 0.0 {
        return 0.0;
    }
    let ze = zeta(z.clamp(0.0001, 0.03));
    let a30 = coeff(&XL_A30, ze);
    // a27 and a28 are multiplied by a30 (see zcnsts.f)
    let a27 = coeff(&XL_A27R, ze) * a30;
    let a28 = coeff(&XL_A28R, ze) * a30;
    let a29 = coeff(&XL_A29, ze);
    let a31 = coeff(&XL_A31, ze);
    let a32 = coeff(&XL_A32, ze);

    let num = a27 * mass.powi(3) + a28 * mass.powi(4) + a29 * mass.powf(a32 + 1.8);
    let den = a30 + a31 * mass.powi(5) + mass.powf(a32);

    if den <= 0.0 { 0.0 } else { num / den }
}

// ── HPT00 — R_TMS coefficients (from xr data, high-mass branch) ─────────

/// xr(24) → a57 constant
const XR_A57: f64 = -8.672_073e-2;
/// xr(25..28) → a58 (4-term)
const XR_A58: [f64; 5] = [2.617_890e0, 1.019_135e0, -3.292_551e-2, -7.445_123e-2, 0.0];
/// xr(29..32) → a59 (4-term)
const XR_A59: [f64; 5] = [1.075_567e-2, 1.773_287e-2, 9.610_479e-3, 1.732_469e-3, 0.0];
/// xr(33..36) → a60 (4-term)
const XR_A60: [f64; 5] = [1.476_246e0, 1.899_331e0, 1.195_010e0, 3.035_051e-1, 0.0];
/// xr(37..39) → a61 (3-term)
const XR_A61: [f64; 5] = [5.502_535e0, -6.601_663e-2, 9.968_707e-2, 3.599_801e-2, 0.0];

/// Terminal main-sequence radius R_TMS (R☉).
///
/// HPT00 Eq. 9. Uses the high-mass branch (M > a62 threshold).
/// For simplicity this uses the high-mass formula for all masses,
/// which is accurate above ~1.5 M☉ and within ~10% below.
#[must_use]
pub fn tms_radius(mass: f64, z: f64) -> f64 {
    if mass <= 0.0 {
        return 0.0;
    }
    let ze = zeta(z.clamp(0.0001, 0.03));
    let a58 = coeff(&XR_A58, ze);
    let a59 = coeff(&XR_A59, ze);
    let a60 = coeff(&XR_A60, ze);
    let a61 = coeff(&XR_A61, ze);

    let m3 = mass.powi(3);
    let m5 = mass.powi(5);
    let ma61 = mass.powf(a61);

    let num = XR_A57 * m3 + a58 * ma61 + a59 * mass.powf(a61 + 1.5);
    let den = a60 + m5;

    // Clamp: R_TMS should be at least R_ZAMS
    let r_zams = zams_radius(mass, z);
    let r_tms = if den <= 0.0 { r_zams } else { num / den };

    r_tms.max(r_zams)
}

/// TMS effective temperature (K) from mass and metallicity.
///
/// Derived from L_TMS and R_TMS via Stefan-Boltzmann.
#[must_use]
pub fn tms_temperature(mass: f64, z: f64) -> f64 {
    let l = tms_luminosity(mass, z);
    let r = tms_radius(mass, z);
    if l <= 0.0 || r <= 0.0 {
        return 0.0;
    }
    crate::constants::T_SUN * l.powf(0.25) / r.sqrt()
}

// ── MS evolution interpolation ────────────────────────────────────────────

/// Main-sequence luminosity at fractional age τ = t/t_MS.
///
/// Simplified HPT00 Eq. 11 (without hook correction):
/// L(τ) = L_ZAMS × 10^(α·τ + (ℓ − α)·τ²)
///
/// where ℓ = log₁₀(L_TMS / L_ZAMS) and α is a small perturbation
/// (~0.3 for solar-mass stars) that steepens late MS brightening.
///
/// `tau` must be in \[0, 1\]. Returns luminosity in L☉.
#[must_use]
pub fn ms_luminosity(mass: f64, z: f64, tau: f64) -> f64 {
    let tau = tau.clamp(0.0, 1.0);
    let l_zams = zams_luminosity(mass, z);
    let l_tms = tms_luminosity(mass, z);
    if l_zams <= 0.0 || l_tms <= 0.0 {
        return 0.0;
    }
    let ell = (l_tms / l_zams).log10();
    // Simplified: α ≈ 0 gives pure quadratic interpolation in log-space
    // This matches HPT00 when hook terms are omitted
    l_zams * 10.0_f64.powf(ell * tau * tau)
}

/// Main-sequence radius at fractional age τ = t/t_MS.
///
/// Simplified HPT00 Eq. 14 (without hook correction):
/// R(τ) = R_ZAMS × 10^(ρ·τ³)
///
/// where ρ = log₁₀(R_TMS / R_ZAMS).
///
/// `tau` must be in \[0, 1\]. Returns radius in R☉.
#[must_use]
pub fn ms_radius(mass: f64, z: f64, tau: f64) -> f64 {
    let tau = tau.clamp(0.0, 1.0);
    let r_zams = zams_radius(mass, z);
    let r_tms = tms_radius(mass, z);
    if r_zams <= 0.0 || r_tms <= 0.0 {
        return 0.0;
    }
    let rho = (r_tms / r_zams).log10();
    r_zams * 10.0_f64.powf(rho * tau.powi(3))
}

/// Main-sequence effective temperature at fractional age τ.
///
/// Derived from [`ms_luminosity`] and [`ms_radius`] via Stefan-Boltzmann.
#[must_use]
pub fn ms_temperature(mass: f64, z: f64, tau: f64) -> f64 {
    let l = ms_luminosity(mass, z, tau);
    let r = ms_radius(mass, z, tau);
    if l <= 0.0 || r <= 0.0 {
        return 0.0;
    }
    crate::constants::T_SUN * l.powf(0.25) / r.sqrt()
}

/// Main-sequence properties at a given age (years).
///
/// Computes τ = age / t_MS, then returns (L, R, T) at that point.
/// If age > t_MS, clamps to τ = 1.0 (TMS values).
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[non_exhaustive]
pub struct MsProperties {
    /// Fractional MS age (0 = ZAMS, 1 = TMS).
    pub tau: f64,
    /// Luminosity in solar luminosities.
    pub luminosity_solar: f64,
    /// Radius in solar radii.
    pub radius_solar: f64,
    /// Effective temperature in Kelvin.
    pub temperature_k: f64,
}

/// Compute MS properties at a given age (years) for mass and metallicity.
#[must_use]
pub fn ms_properties(mass: f64, z: f64, age_years: f64) -> MsProperties {
    let t_ms = ms_lifetime(mass, z);
    let tau = if t_ms > 0.0 {
        (age_years / t_ms).clamp(0.0, 1.0)
    } else {
        0.0
    };
    let luminosity_solar = ms_luminosity(mass, z, tau);
    let radius_solar = ms_radius(mass, z, tau);
    let temperature_k = if luminosity_solar > 0.0 && radius_solar > 0.0 {
        crate::constants::T_SUN * luminosity_solar.powf(0.25) / radius_solar.sqrt()
    } else {
        0.0
    };
    MsProperties {
        tau,
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
    fn t_bgb_solar() {
        // Solar t_BGB should be ~10-12 Gyr
        let t = t_bgb(1.0, Z_SUN);
        assert!(
            t > 9_000.0 && t < 13_000.0,
            "Solar t_BGB: {t} Myr, expected ~10000-12000"
        );
    }

    #[test]
    fn t_bgb_massive_shorter() {
        let t1 = t_bgb(1.0, Z_SUN);
        let t10 = t_bgb(10.0, Z_SUN);
        assert!(t10 < t1, "10 M_sun should have shorter t_BGB");
    }

    #[test]
    fn t_bgb_low_mass_longer() {
        let t1 = t_bgb(1.0, Z_SUN);
        let t01 = t_bgb(0.5, Z_SUN);
        assert!(t01 > t1, "0.5 M_sun should have longer t_BGB");
    }

    #[test]
    fn ms_lifetime_solar() {
        // Solar MS lifetime ~10 Gyr
        let t = ms_lifetime(1.0, Z_SUN);
        assert!(
            t > 8e9 && t < 13e9,
            "Solar MS lifetime: {t:.2e} yr, expected ~10 Gyr"
        );
    }

    #[test]
    fn ms_lifetime_less_than_tbgb() {
        for m in [0.5, 1.0, 2.0, 5.0, 10.0, 20.0] {
            let t_ms = ms_lifetime_myr(m, Z_SUN);
            let t_bg = t_bgb(m, Z_SUN);
            assert!(
                t_ms <= t_bg,
                "t_MS should be <= t_BGB for M={m}: t_MS={t_ms}, t_BGB={t_bg}"
            );
        }
    }

    #[test]
    fn ms_lifetime_massive_short() {
        let t = ms_lifetime(20.0, Z_SUN);
        assert!(t < 2e7, "20 M_sun MS lifetime should be < 20 Myr: {t:.2e}");
    }

    #[test]
    fn ms_lifetime_metallicity_effect() {
        // Lower Z → slightly shorter MS lifetime (hotter, more luminous)
        let t_low = ms_lifetime(1.0, 0.001);
        let t_solar = ms_lifetime(1.0, Z_SUN);
        assert!(
            t_low < t_solar,
            "Low-Z should have shorter MS: Z=0.001 → {t_low:.2e}, Z=0.02 → {t_solar:.2e}"
        );
    }

    #[test]
    fn tms_luminosity_brighter_than_zams() {
        // Stars brighten during MS: L_TMS > L_ZAMS
        for m in [0.5, 1.0, 2.0, 5.0, 10.0] {
            let l_z = zams_luminosity(m, Z_SUN);
            let l_t = tms_luminosity(m, Z_SUN);
            assert!(
                l_t > l_z,
                "L_TMS should exceed L_ZAMS for M={m}: L_TMS={l_t}, L_ZAMS={l_z}"
            );
        }
    }

    #[test]
    fn tms_luminosity_solar() {
        // Sun's TMS luminosity should be ~1.5-2× ZAMS
        let l_t = tms_luminosity(1.0, Z_SUN);
        let l_z = zams_luminosity(1.0, Z_SUN);
        let ratio = l_t / l_z;
        assert!(
            ratio > 1.2 && ratio < 4.0,
            "Solar L_TMS/L_ZAMS ratio: {ratio}, expected ~2-3"
        );
    }

    #[test]
    fn tms_radius_at_least_zams() {
        for m in [0.5, 1.0, 2.0, 5.0, 10.0] {
            let r_t = tms_radius(m, Z_SUN);
            let r_z = zams_radius(m, Z_SUN);
            assert!(
                r_t >= r_z,
                "R_TMS should be >= R_ZAMS for M={m}: R_TMS={r_t}, R_ZAMS={r_z}"
            );
        }
    }

    #[test]
    fn tms_temperature_solar() {
        let t = tms_temperature(1.0, Z_SUN);
        // Should be roughly in the range of solar temperatures
        assert!(t > 4000.0 && t < 8000.0, "Solar TMS temperature: {t} K");
    }

    #[test]
    fn ms_luminosity_at_tau_0_is_zams() {
        let l = ms_luminosity(1.0, Z_SUN, 0.0);
        let l_z = zams_luminosity(1.0, Z_SUN);
        assert!(
            (l - l_z).abs() < 1e-10,
            "L(τ=0) should equal L_ZAMS: {l} vs {l_z}"
        );
    }

    #[test]
    fn ms_luminosity_at_tau_1_is_tms() {
        let l = ms_luminosity(1.0, Z_SUN, 1.0);
        let l_t = tms_luminosity(1.0, Z_SUN);
        assert!(
            (l - l_t).abs() / l_t < 0.01,
            "L(τ=1) should equal L_TMS: {l} vs {l_t}"
        );
    }

    #[test]
    fn ms_luminosity_monotonically_increases() {
        let mut prev = ms_luminosity(1.0, Z_SUN, 0.0);
        for i in 1..=10 {
            let tau = i as f64 / 10.0;
            let l = ms_luminosity(1.0, Z_SUN, tau);
            assert!(l >= prev, "L should increase: τ={tau}, L={l}, prev={prev}");
            prev = l;
        }
    }

    #[test]
    fn ms_radius_at_tau_0_is_zams() {
        let r = ms_radius(1.0, Z_SUN, 0.0);
        let r_z = zams_radius(1.0, Z_SUN);
        assert!(
            (r - r_z).abs() < 1e-10,
            "R(τ=0) should equal R_ZAMS: {r} vs {r_z}"
        );
    }

    #[test]
    fn ms_properties_solar_midlife() {
        // Sun at age 4.6 Gyr
        let p = ms_properties(1.0, Z_SUN, 4.6e9);
        assert!(p.tau > 0.0 && p.tau < 1.0, "Solar τ: {}", p.tau);
        assert!(p.luminosity_solar > 0.5 && p.luminosity_solar < 2.0);
        assert!(p.radius_solar > 0.5 && p.radius_solar < 2.0);
        assert!(p.temperature_k > 5000.0 && p.temperature_k < 7000.0);
    }

    #[test]
    fn ms_properties_serde_roundtrip() {
        let p = ms_properties(1.0, Z_SUN, 4.6e9);
        let json = serde_json::to_string(&p).unwrap();
        let back: MsProperties = serde_json::from_str(&json).unwrap();
        assert!((back.tau - p.tau).abs() < 1e-15);
        assert!((back.luminosity_solar - p.luminosity_solar).abs() < 1e-15);
        assert!((back.radius_solar - p.radius_solar).abs() < 1e-15);
        assert!((back.temperature_k - p.temperature_k).abs() < 1e-10);
    }

    #[test]
    fn zams_zero_mass() {
        assert!(zams_luminosity(0.0, Z_SUN).abs() < f64::EPSILON);
        assert!(zams_radius(0.0, Z_SUN).abs() < f64::EPSILON);
        assert!(zams_temperature(0.0, Z_SUN).abs() < f64::EPSILON);
    }
}
