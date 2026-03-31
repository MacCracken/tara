use tara::atmosphere;
use tara::bridge;
use tara::classification::{self, LuminosityClass};
use tara::evolution;
use tara::luminosity;
use tara::nucleosynthesis;
use tara::spectral;
use tara::star::{SpectralClass, Star};

#[test]
fn crate_loads() {
    let _ = std::any::type_name::<tara::error::TaraError>();
}

/// End-to-end: build a Sun, verify classification, evolution, and luminosity agree.
#[test]
fn sun_end_to_end() {
    let sun = Star::sun().unwrap();

    // Classification matches
    assert_eq!(sun.spectral_class, SpectralClass::G);
    assert_eq!(sun.spectral_subclass, 2);
    assert_eq!(sun.luminosity_class, LuminosityClass::V);

    // MK string
    let mk = classification::format_classification(sun.temperature_k, sun.luminosity_class);
    assert_eq!(mk, "G2V");

    // Luminosity from radius + temp should be ~1 L_sun
    let l = luminosity::luminosity_solar(sun.radius_solar, sun.temperature_k);
    assert!((l - 1.0).abs() < 0.02, "Solar luminosity: {l}");

    // Absolute bolometric magnitude
    let m_bol = luminosity::absolute_bolometric_magnitude(l);
    assert!((m_bol - 4.74).abs() < 0.1, "Solar M_bol: {m_bol}");

    // Evolutionary phase at solar age
    let phase = evolution::evolutionary_phase(sun.mass_solar, sun.age_years);
    assert_eq!(phase, evolution::EvolutionaryPhase::MainSequence);

    // Remnant type
    assert_eq!(
        evolution::remnant_type(sun.mass_solar),
        evolution::RemnantType::WhiteDwarf
    );
}

/// Build a massive O-star and verify cross-module consistency.
#[test]
fn massive_star_cross_module() {
    let star = Star::builder(30.0, 10.0, 35_000.0, 200_000.0, 3e6)
        .build()
        .unwrap();

    assert_eq!(star.spectral_class, SpectralClass::O);

    // Short-lived
    let lifetime = evolution::main_sequence_lifetime(star.mass_solar);
    assert!(
        lifetime < 1e8,
        "O-star lifetime should be < 100 Myr: {lifetime:.2e}"
    );

    // Should become a black hole
    assert_eq!(
        evolution::remnant_type(star.mass_solar),
        evolution::RemnantType::BlackHole
    );

    // HR region
    let region = classification::hr_region(star.temperature_k, star.luminosity_solar);
    assert_eq!(region, classification::HrRegion::Supergiant);

    // Nucleosynthesis: CNO dominant at ~30 MK core
    assert_eq!(
        nucleosynthesis::dominant_process(30e6),
        nucleosynthesis::FusionProcess::CnoCycle
    );
}

/// Bridge module: round-trip through falak/prakash/tanmatra bridges.
#[test]
fn bridge_round_trips() {
    // Falak: mass → gravitational parameter
    let mu = bridge::stellar_mass_solar_to_mu(1.0);
    assert!(mu > 1e20, "Solar mu: {mu}");

    // Prakash: temperature → peak wavelength
    let peak = bridge::temperature_to_peak_wavelength_nm(5772.0);
    assert!((peak - 502.0).abs() < 5.0, "Solar peak: {peak} nm");

    // Limb darkening from log(g)
    let log_g = atmosphere::log_surface_gravity(1.0, 1.0);
    let u = bridge::surface_gravity_to_limb_darkening(log_g);
    assert!(u > 0.4 && u < 0.7, "Solar limb darkening: {u}");

    // Tanmatra: composition → lifetime scale
    let scale = bridge::composition_to_lifetime_scale(0.74);
    assert!(
        (scale - 1.0).abs() < 0.01,
        "Solar composition scale: {scale}"
    );

    // Tanmatra: core temp → CNO fraction
    let cno = bridge::core_temp_to_cno_fraction(15.7e6);
    assert!(cno < 0.5, "Solar core should be pp-dominated: {cno}");

    // Bridge: luminosity class → abs magnitude
    let m = bridge::luminosity_class_to_abs_magnitude(LuminosityClass::V);
    assert!((m - 5.0).abs() < 0.1);
}

/// Spectral: Planck → Doppler → broadening pipeline.
#[test]
fn spectral_pipeline() {
    // Planck radiance at H-alpha for solar temperature
    let h_alpha_nm = 656.3;
    let radiance = spectral::planck_radiance_nm(h_alpha_nm, 5772.0);
    assert!(radiance > 0.0, "H-alpha radiance: {radiance}");

    // Doppler shift at 100 km/s recession
    let shifted = spectral::doppler_shift(h_alpha_nm, 1e5);
    assert!(shifted > h_alpha_nm, "Redshifted: {shifted}");

    // Thermal broadening of H-alpha
    let sigma = spectral::thermal_broadening(656.3e-9, 5772.0, tara::constants::AMU);
    assert!(sigma > 0.0 && sigma < 1e-9, "Broadening sigma: {sigma}");

    // Voigt profile should produce non-zero value at line center
    let v = spectral::pseudo_voigt_profile(h_alpha_nm, h_alpha_nm, 0.01, 0.005);
    assert!(v > 0.0, "Voigt at center: {v}");
}

/// Atmosphere: T_eff → log(g) → grey atmosphere → opacity chain.
#[test]
fn atmosphere_chain() {
    let t_eff = atmosphere::effective_temperature(tara::constants::L_SUN, tara::constants::R_SUN);
    assert!((t_eff - 5772.0).abs() < 10.0, "Effective temp: {t_eff}");

    let log_g = atmosphere::log_surface_gravity(1.0, 1.0);
    assert!((log_g - 4.44).abs() < 0.02, "log(g): {log_g}");

    // Grey atmosphere: at τ=2/3 should recover ~T_eff
    let t_23 = atmosphere::grey_atmosphere_temperature(t_eff, 2.0 / 3.0);
    // T(τ=2/3) = (3/4 × T_eff⁴ × 4/3)^(1/4) = T_eff
    assert!((t_23 - t_eff).abs() < 1.0, "Grey T at τ=2/3: {t_23}");

    // Electron scattering for solar composition
    let kappa_es = atmosphere::electron_scattering_opacity(0.74);
    assert!((kappa_es - 0.348).abs() < 0.01, "kappa_es: {kappa_es}");

    // Scale height should be order ~100 km
    let g_cgs = 10.0_f64.powf(log_g); // cm/s²
    let g_si = g_cgs / 100.0; // m/s²
    let h = atmosphere::scale_height(t_eff, 1.3, g_si);
    assert!(h > 1e5 && h < 3e5, "Scale height: {h} m");
}

/// Luminosity: distance modulus ↔ apparent/absolute magnitude consistency.
#[test]
fn magnitude_distance_consistency() {
    // Sun at 10 pc: m = M (distance modulus = 0)
    let m_abs = luminosity::absolute_bolometric_magnitude(1.0);
    let m_app = luminosity::apparent_magnitude(m_abs, 10.0);
    assert!(
        (m_app - m_abs).abs() < f64::EPSILON,
        "At 10 pc, m should equal M"
    );

    // At 100 pc, apparent should be 5 mag dimmer
    let m_100 = luminosity::apparent_magnitude(m_abs, 100.0);
    assert!(
        (m_100 - m_abs - 5.0).abs() < 1e-10,
        "100 pc distance modulus: {}",
        m_100 - m_abs
    );

    // Round-trip: magnitude → luminosity → magnitude
    let l = luminosity::luminosity_from_abs_bol_mag(m_abs);
    let m_back = luminosity::absolute_bolometric_magnitude(l);
    assert!(
        (m_back - m_abs).abs() < 1e-10,
        "Magnitude roundtrip: {m_abs} → L={l} → {m_back}"
    );
}
