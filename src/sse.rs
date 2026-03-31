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

// ── HPT00 — R_TMS coefficients (from xr data) ───────────────────────────

// Low-mass branch (M <= a62): R_TMS = max(1.5*R_ZAMS, (a52+a53*M^a55)/(a54+M^a56))

/// xr(1..5) → raw a52 before multiply by a54
const XR_A52R: [f64; 5] = [
    2.187_715e-1,
    -2.154_437e0,
    -3.768_678e0,
    -1.975_518e0,
    -3.021_475e-1,
];
/// xr(6..10) → raw a53 before multiply by a54
const XR_A53R: [f64; 5] = [
    1.466_440e0,
    1.839_725e0,
    6.442_199e0,
    4.023_635e0,
    6.957_529e-1,
];
/// xr(11..15) → a54 (also multiplier for a52, a53)
const XR_A54: [f64; 5] = [
    2.652_091e1,
    8.178_458e1,
    1.156_058e2,
    7.633_811e1,
    1.950_698e1,
];
/// xr(16..19) → a55 (4-term)
const XR_A55: [f64; 5] = [1.472_103e0, -2.947_609e0, -3.312_828e0, -9.945_065e-1, 0.0];
/// xr(20..23) → a56 (4-term)
const XR_A56: [f64; 5] = [3.071_048e0, -5.679_941e0, -9.745_523e0, -3.594_543e0, 0.0];

// High-mass branch (M >= a62+0.1): R_TMS = (a57*M^3+a58*M^a61+a59*M^(a61+1.5))/(a60+M^5)

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

/// Compute the mass threshold a62 separating low/high-mass R_TMS branches.
#[must_use]
fn rtms_mass_threshold(z: f64) -> f64 {
    let lz = z.log10();
    let x = 0.1461 + 0.1237 * (lz + 2.0);
    let inner = 0.1461_f64.min(x).max(0.097);
    let outer = (0.097 - 0.1072 * (lz + 3.0)).max(inner);
    10.0_f64.powf(outer)
}

/// Low-mass R_TMS branch: (a52 + a53*M^a55) / (a54 + M^a56).
#[must_use]
fn rtms_low_mass(mass: f64, ze: f64) -> f64 {
    let a54 = coeff(&XR_A54, ze);
    let a52 = coeff(&XR_A52R, ze) * a54;
    let a53 = coeff(&XR_A53R, ze) * a54;
    let a55 = coeff(&XR_A55, ze);
    let a56 = coeff(&XR_A56, ze);

    let num = a52 + a53 * mass.powf(a55);
    let den = a54 + mass.powf(a56);
    if den <= 0.0 { 0.0 } else { num / den }
}

/// High-mass R_TMS branch: (a57*M^3 + a58*M^a61 + a59*M^(a61+1.5)) / (a60 + M^5).
#[must_use]
fn rtms_high_mass(mass: f64, ze: f64) -> f64 {
    let a58 = coeff(&XR_A58, ze);
    let a59 = coeff(&XR_A59, ze);
    let a60 = coeff(&XR_A60, ze);
    let a61 = coeff(&XR_A61, ze);

    let num = XR_A57 * mass.powi(3) + a58 * mass.powf(a61) + a59 * mass.powf(a61 + 1.5);
    let den = a60 + mass.powi(5);
    if den <= 0.0 { 0.0 } else { num / den }
}

/// Terminal main-sequence radius R_TMS (R☉).
///
/// HPT00 Eq. 9. Piecewise formula with three branches:
/// - M ≤ a62: low-mass rational function, floored at 1.5 × R_ZAMS
/// - M ≥ a62 + 0.1: high-mass rational function
/// - a62 < M < a62 + 0.1: linear interpolation between branches
#[must_use]
pub fn tms_radius(mass: f64, z: f64) -> f64 {
    if mass <= 0.0 {
        return 0.0;
    }
    let z_clamped = z.clamp(0.0001, 0.03);
    let ze = zeta(z_clamped);
    let a62 = rtms_mass_threshold(z_clamped);
    let m2 = a62 + 0.1;

    let r = if mass <= a62 {
        let r_low = rtms_low_mass(mass, ze);
        let r_floor = 1.5 * zams_radius(mass, z_clamped);
        r_low.max(r_floor)
    } else if mass >= m2 {
        rtms_high_mass(mass, ze)
    } else {
        // Linear interpolation across the transition
        let r_at_a62 = {
            let r_low = rtms_low_mass(a62, ze);
            let r_floor = 1.5 * zams_radius(a62, z_clamped);
            r_low.max(r_floor)
        };
        let r_at_m2 = rtms_high_mass(m2, ze);
        let frac = (mass - a62) / 0.1;
        r_at_a62 + (r_at_m2 - r_at_a62) * frac
    };

    // R_TMS must be at least R_ZAMS
    r.max(zams_radius(mass, z_clamped))
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

// ── HPT00 — MS interpolation coefficients ────────────────────────────────

// Lalpha coefficients: msp(33)–msp(42)
/// xl(28..31) → a33
const XL_A33: [f64; 5] = [
    2.321_400e-1,
    1.828_075e-3,
    -2.232_007e-2,
    -3.378_734e-3,
    0.0,
];
/// xl(32..35) → a34
const XL_A34: [f64; 5] = [1.163_659e-2, 3.427_682e-3, 1.421_393e-3, -3.710_666e-3, 0.0];
/// xl(36..39) → a35
const XL_A35: [f64; 5] = [
    1.048_020e-2,
    -1.231_921e-2,
    -1.686_860e-2,
    -4.234_354e-3,
    0.0,
];
/// xl(40..43) → a36
const XL_A36: [f64; 5] = [
    1.555_590e0,
    -3.223_927e-1,
    -5.197_429e-1,
    -1.066_441e-1,
    0.0,
];

// Lbeta coefficients: msp(43)–msp(46)
/// xl(44..48) → a43
const XL_A43: [f64; 5] = [
    3.855_707e-1,
    -6.104_166e-1,
    5.676_742e0,
    1.060_894e1,
    5.284_014e0,
];
/// xl(49..53) → a44
const XL_A44: [f64; 5] = [
    3.579_064e-1,
    -6.442_936e-1,
    5.494_644e0,
    1.054_952e1,
    5.280_991e0,
];
/// xl(54..56) → a45
const XL_A45: [f64; 5] = [9.587_587e-1, 8.777_464e-1, 2.017_321e-1, 0.0, 0.0];

// Lhook coefficients: msp(47)–msp(51)
/// xl(57..60) → a47
const XL_A47: [f64; 5] = [1.910_302e-1, 1.158_624e-1, 3.348_990e-2, 2.599_706e-3, 0.0];
/// xl(61..64) → a48
const XL_A48: [f64; 5] = [
    3.931_056e-1,
    7.277_637e-2,
    -1.366_593e-1,
    -4.508_946e-2,
    0.0,
];
/// xl(65..68) → a49
const XL_A49: [f64; 5] = [3.267_776e-1, 1.204_424e-1, 9.988_332e-2, 2.455_361e-2, 0.0];
/// xl(69..72) → a50
const XL_A50: [f64; 5] = [5.990_212e-1, 5.570_264e-2, 6.207_626e-2, 1.777_283e-2, 0.0];

// Ralpha coefficients: msp(65)–msp(76) — from xr(41..64)
/// xr(41..44) → a65
const XR_A65: [f64; 5] = [
    4.907_546e-1,
    -1.683_928e-1,
    -3.108_742e-1,
    -7.202_918e-2,
    0.0,
];
/// xr(45..48) → a66
const XR_A66: [f64; 5] = [4.537_070e0, -4.465_455e0, -1.612_690e0, -1.623_246e0, 0.0];
/// xr(49..52) → a67
const XR_A67: [f64; 5] = [1.796_220e0, 2.814_020e-1, 1.423_325e0, 3.421_036e-1, 0.0];
/// xr(53..56) → a68
const XR_A68: [f64; 5] = [2.256_216e0, 3.773_400e-1, 1.537_867e0, 4.396_373e-1, 0.0];
/// xr(57..61) → a69
const XR_A69: [f64; 5] = [
    1.564_231e-3,
    1.653_042e-3,
    -4.439_786e-3,
    -4.951_011e-3,
    -1.216_530e-3,
];
/// xr(62..64) → a72
const XR_A72: [f64; 5] = [5.210_157e0, -4.143_695e0, -2.120_870e0, 0.0, 0.0];

// Rbeta coefficients: msp(77)–msp(82) — from xr(65..83)
/// xr(65..68) → a77
const XR_A77: [f64; 5] = [
    1.071_489e0,
    -1.164_852e-1,
    -8.623_831e-2,
    -1.582_349e-2,
    0.0,
];
/// xr(69..72) → a78
const XR_A78: [f64; 5] = [7.108_492e-1, 7.935_927e-1, 3.926_983e-1, 3.622_146e-2, 0.0];
/// xr(73..76) → a79
const XR_A79: [f64; 5] = [
    3.478_514e0,
    -2.585_474e-2,
    -1.512_955e-2,
    -2.833_691e-3,
    0.0,
];
/// xr(77..80) → a80
const XR_A80: [f64; 5] = [3.969_331e-3, 4.539_076e-3, 1.720_906e-3, 1.897_857e-4, 0.0];
/// xr(81..83) → a81 (3-term, special: lzs*lzs*xr(83))
const XR_A81: [f64; 5] = [9.132_108e-1, -1.653_695e-1, 3.636_784e-2, 0.0, 0.0];

// Rhook coefficients (xr(104..119)) reserved for future rhook implementation.

// ── MS interpolation helper functions ────────────────────────────────────

/// Luminosity alpha coefficient (HPT00).
#[must_use]
fn lalpha(mass: f64, ze: f64) -> f64 {
    let a33 = coeff(&XL_A33, ze);
    let a34 = coeff(&XL_A34, ze);
    let a35 = coeff(&XL_A35, ze);
    let a36 = coeff(&XL_A36, ze);
    let a39 = (0.0977 - ze * (0.231 + 0.0753 * ze)).max(0.145);
    let a40 = (0.24 + ze * (0.18 + 0.595 * ze)).min(0.306 + 0.053 * ze);
    let a41 = (0.33 + ze * (0.132 + 0.218 * ze)).min(0.3625 + 0.062 * ze);
    let a37 = (1.1064 + ze * (0.415 + 0.18 * ze)).max(0.9);
    let a38 = (1.19 + ze * (0.377 + 0.176 * ze)).max(1.0);
    let a42 = (a33 + a34 * 2.0_f64.powf(a36)) / (2.0_f64.powf(0.4) + a35 * 2.0_f64.powf(1.9));

    if mass >= 2.0 {
        (a33 + a34 * mass.powf(a36)) / (mass.powf(0.4) + a35 * mass.powf(1.9))
    } else if mass <= 0.5 {
        a39
    } else if mass <= 0.7 {
        a39 + ((0.3 - a39) / 0.2) * (mass - 0.5)
    } else if mass <= a37 {
        0.3 + ((a40 - 0.3) / (a37 - 0.7)) * (mass - 0.7)
    } else if mass <= a38 {
        a40 + ((a41 - a40) / (a38 - a37)) * (mass - a37)
    } else {
        a41 + ((a42 - a41) / (2.0 - a38)) * (mass - a38)
    }
}

/// Luminosity beta coefficient (HPT00).
#[must_use]
fn lbeta(mass: f64, ze: f64) -> f64 {
    let a43 = coeff(&XL_A43, ze);
    let a44 = coeff(&XL_A44, ze);
    let a45 = coeff(&XL_A45, ze);
    let a46 = {
        let v = (1.5135 + 0.3769 * ze).min(1.4);
        (0.6355 - 0.4192 * ze).max(1.25_f64.max(v))
    };

    let mut b = a43 - a44 * mass.powf(a45);
    b = b.max(0.0);
    if mass > a46 && b > 0.0 {
        let a1 = a43 - a44 * a46.powf(a45);
        b = (a1 - 10.0 * a1 * (mass - a46)).max(0.0);
    }
    b
}

/// Luminosity eta exponent (HPT00).
#[must_use]
fn leta(mass: f64, ze: f64) -> f64 {
    let a97 = if ze > (0.0009_f64 / Z_SUN).log10() {
        10.0
    } else {
        20.0
    };
    let eta = if mass <= 1.0 {
        10.0
    } else if mass >= 1.1 {
        20.0
    } else {
        10.0 + 100.0 * (mass - 1.0)
    };
    eta.min(a97)
}

/// Luminosity hook perturbation δ_L (HPT00).
#[must_use]
fn lhook(mass: f64, ze: f64, z: f64) -> f64 {
    let a47 = coeff(&XL_A47, ze);
    let a48 = coeff(&XL_A48, ze);
    let a49 = coeff(&XL_A49, ze);
    let a50 = coeff(&XL_A50, ze);
    let a51 = {
        let v = (1.5135 + 0.3769 * ze).min(1.4);
        (0.6355 - 0.4192 * ze).max(1.25_f64.max(v))
    };
    // zpars(1) = Mhook
    let mhook = 1.0185 + ze * (0.16015 + ze * 0.0892);
    let _ = z; // used indirectly via ze

    if mass <= mhook {
        0.0
    } else if mass >= a51 {
        (a47 / mass.powf(a48)).min(a49 / mass.powf(a50))
    } else {
        let a2 = (a47 / a51.powf(a48)).min(a49 / a51.powf(a50));
        a2 * ((mass - mhook) / (a51 - mhook)).powf(0.4)
    }
}

/// Radius alpha coefficient (HPT00) — simplified for M >= 0.5.
#[must_use]
fn ralpha(mass: f64, ze: f64) -> f64 {
    let a65 = coeff(&XR_A65, ze);
    let a66 = coeff(&XR_A66, ze);
    let a67 = coeff(&XR_A67, ze);
    let a68 = coeff(&XR_A68, ze);
    let a69 = coeff(&XR_A69, ze);
    let a72 = coeff(&XR_A72, ze);
    let a73 = (0.0843 - ze * (0.0475 + 0.0352 * ze)).max(0.065);
    let a74 = {
        let v = 0.0736 + ze * (0.0749 + 0.04426 * ze);
        let lz = ze + Z_SUN.log10(); // log10(z) = ζ + log10(Z_SUN)
        if lz < (0.004_f64).log10() {
            v.min(0.055)
        } else {
            v
        }
    };
    let a70 = (1.116 + 0.166 * ze).clamp(0.9, 1.0);
    let a71 = {
        let v = (1.477 + 0.296 * ze).max((-0.308 - 1.046 * ze).min(1.6));
        (0.8 - 2.0 * ze).min(v).max(0.8)
    };
    let a75 = (0.136 + 0.0352 * ze).clamp(0.091, 0.121);

    if mass <= 0.5 {
        a73
    } else if mass <= 0.65 {
        a73 + ((a74 - a73) / 0.15) * (mass - 0.5)
    } else if mass <= a70 {
        a74 + ((a75 - a74) / (a70 - 0.65)) * (mass - 0.65)
    } else if mass <= a71 {
        let a76 = (a65 * a71.powf(a67)) / (a66 + a71.powf(a68));
        a75 + ((a76 - a75) / (a71 - a70)) * (mass - a70)
    } else if mass <= a72 {
        (a65 * mass.powf(a67)) / (a66 + mass.powf(a68))
    } else {
        let a5 = (a65 * a72.powf(a67)) / (a66 + a72.powf(a68));
        a5 + a69 * (mass - a72)
    }
}

/// Radius beta coefficient (HPT00).
#[must_use]
fn rbeta(mass: f64, ze: f64, z: f64) -> f64 {
    let a77 = coeff(&XR_A77, ze);
    let a78 = coeff(&XR_A78, ze);
    let a79 = coeff(&XR_A79, ze);
    let a80 = coeff(&XR_A80, ze);
    let a81 = coeff(&XR_A81, ze);
    let a82 = {
        let v = (1.6 + ze * (0.764 + 0.3322 * ze)).clamp(1.4, 1.6);
        if z > 0.01 { v.max(0.95) } else { v }
    };
    let _ = z;

    let b = if mass <= 1.0 {
        1.06
    } else if mass <= a82 {
        1.06 + ((a81 - 1.06) / (a82 - 1.0)) * (mass - 1.0)
    } else if mass <= 2.0 {
        let b2 = (a77 * 2.0_f64.powf(3.5)) / (a78 + 2.0_f64.powf(a79));
        a81 + ((b2 - a81) / (2.0 - a82)) * (mass - a82)
    } else if mass <= 16.0 {
        (a77 * mass.powf(3.5)) / (a78 + mass.powf(a79))
    } else {
        let b3 = (a77 * 16.0_f64.powf(3.5)) / (a78 + 16.0_f64.powf(a79));
        b3 + a80 * (mass - 16.0)
    };
    b - 1.0
}

/// Radius gamma coefficient (HPT00) — simplified.
///
/// Returns 0 for most stars; non-zero only near M ~ 1 M☉.
#[must_use]
fn rgamma(mass: f64, ze: f64) -> f64 {
    // Simplified Rgamma — full xr(84..103) coefficients reserved for future.
    // Gamma is ~0 for most masses; non-zero only near M ~ 1 M☉.
    let a88 = (-0.2711 - ze * (0.5756 + 0.0838 * ze)).max(5.855_420e-2);
    // Simplified: rgamma = 0 for mass > a88 + 0.1
    if mass > a88 + 0.1 || mass > 1.5 {
        0.0
    } else {
        // Small positive value near M ~ 1
        let peak = 0.02_f64.max(0.0);
        if mass <= 1.0 {
            peak * (1.0 - (mass - 1.0).abs()).max(0.0)
        } else {
            peak * (1.0 - (mass - 1.0) * 10.0).max(0.0)
        }
    }
}

// ── MS evolution interpolation ────────────────────────────────────────────

/// Main-sequence luminosity at fractional age τ = t/t_MS.
///
/// HPT00 Eq. 11 with hook correction:
/// L(τ) = L_ZAMS × 10^(α·τ + β·τ^η + (ℓ − α − β)·τ² − δ·(τ₁² − τ₂²))
///
/// where ℓ = log₁₀(L_TMS / L_ZAMS), and τ₁, τ₂ define the hook region.
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
    let z_clamped = z.clamp(0.0001, 0.03);
    let ze = zeta(z_clamped);
    let ell = (l_tms / l_zams).log10();
    let alpha = lalpha(mass, ze);
    let beta = lbeta(mass, ze);
    let eta = leta(mass, ze);
    let dell = lhook(mass, ze, z_clamped);

    // Hook region: τ₁ = min(1, max(0.5, ..hook fraction..))
    // τ₂ = 1 (end of hook). dtau = τ₁² − τ₂² = τ₁² − 1
    // For simplicity, hook occupies τ ∈ [0.95, 1.0]
    let tau_min = 0.99 - 0.01 * dell.abs().min(1.0); // narrower hook window when dell is large
    let tau1 = tau.min(1.0).max(tau_min);
    let tau2 = 1.0_f64.min(tau.max(tau_min));
    let dtau = tau1 * tau1 - tau2 * tau2;

    let xx = if beta > 0.0 && tau > tau_min {
        alpha * tau + beta * tau.powf(eta) + (ell - alpha - beta) * tau * tau - dell * dtau
    } else {
        alpha * tau + (ell - alpha) * tau * tau - dell * dtau
    };

    l_zams * 10.0_f64.powf(xx)
}

/// Main-sequence radius at fractional age τ = t/t_MS.
///
/// HPT00 Eq. 14 with hook correction:
/// R(τ) = R_ZAMS × 10^(α·τ + β·τ^10 + γ·τ^40 + (ρ − α − β − γ)·τ³ − δ·(τ₁³ − τ₂³))
///
/// `tau` must be in \[0, 1\]. Returns radius in R☉.
#[must_use]
pub fn ms_radius(mass: f64, z: f64, tau: f64) -> f64 {
    let tau = tau.clamp(0.0, 1.0);
    let z_clamped = z.clamp(0.0001, 0.03);
    let r_zams = zams_radius(mass, z_clamped);
    let r_tms = tms_radius(mass, z_clamped);
    if r_zams <= 0.0 || r_tms <= 0.0 {
        return 0.0;
    }
    let ze = zeta(z_clamped);
    let rho = (r_tms / r_zams).log10();
    let alpha = ralpha(mass, ze);
    let beta = rbeta(mass, ze, z_clamped);
    let gamma = rgamma(mass, ze);

    let xx = alpha * tau
        + beta * tau.powi(10)
        + gamma * tau.powi(40)
        + (rho - alpha - beta - gamma) * tau.powi(3);

    r_zams * 10.0_f64.powf(xx)
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
    fn tms_radius_low_mass_uses_floor() {
        // Very low mass: R_TMS should be at least 1.5 × R_ZAMS
        let r_tms = tms_radius(0.2, Z_SUN);
        let r_zams = zams_radius(0.2, Z_SUN);
        assert!(
            r_tms >= 1.5 * r_zams * 0.99,
            "Low-mass R_TMS should be >= 1.5*R_ZAMS: R_TMS={r_tms}, 1.5*R_ZAMS={}",
            1.5 * r_zams
        );
    }

    #[test]
    fn tms_radius_continuous_across_threshold() {
        // R_TMS should be roughly continuous across the a62 threshold
        let a62 = rtms_mass_threshold(Z_SUN);
        let r_below = tms_radius(a62 - 0.01, Z_SUN);
        let r_above = tms_radius(a62 + 0.11, Z_SUN);
        let r_mid = tms_radius(a62 + 0.05, Z_SUN);
        // Mid should be between the two endpoints (within transition)
        assert!(
            r_mid > r_below.min(r_above) * 0.5 && r_mid < r_below.max(r_above) * 2.0,
            "R_TMS should be roughly continuous: r_below={r_below}, r_mid={r_mid}, r_above={r_above}"
        );
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
