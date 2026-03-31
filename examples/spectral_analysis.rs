//! Spectral analysis: blackbody curves, Doppler shifts, and line profiles.

use tara::constants;
use tara::spectral;
use tara::star::Star;

fn main() {
    let sun = Star::sun().unwrap();

    // Blackbody spectrum sampling
    println!("=== Solar Blackbody Spectrum ===");
    println!("{:>10} {:>15}", "λ (nm)", "B(λ,T) (W/sr/m³)");
    for lam_nm in (200..=1200).step_by(100) {
        let b = spectral::planck_radiance_nm(lam_nm as f64, sun.temperature_k);
        println!("{lam_nm:>10} {b:>15.4e}");
    }

    // Wien peak
    let peak_nm = constants::WIEN_B_NM / sun.temperature_k;
    println!("\nWien peak wavelength: {peak_nm:.1} nm");

    // Doppler shifts
    println!("\n=== Doppler Shifts (H-alpha 656.3 nm) ===");
    let h_alpha = 656.3;
    let velocities = [-300_000.0, -100_000.0, 0.0, 100_000.0, 300_000.0, 1e8];
    for v in velocities {
        let shifted = spectral::doppler_shift(h_alpha, v);
        let rel_shifted = spectral::doppler_shift_relativistic(h_alpha, v);
        println!(
            "v = {:>+12.0} m/s → classical: {shifted:>8.3} nm, relativistic: {rel_shifted:>8.3} nm",
            v
        );
    }

    // Thermal broadening
    println!("\n=== Thermal Broadening ===");
    let sigma_h = spectral::thermal_broadening(656.3e-9, sun.temperature_k, constants::AMU);
    let sigma_fe =
        spectral::thermal_broadening(656.3e-9, sun.temperature_k, 55.845 * constants::AMU);
    println!(
        "H-alpha broadening (hydrogen): σ = {:.4e} m = {:.4} pm",
        sigma_h,
        sigma_h * 1e12
    );
    println!(
        "H-alpha broadening (iron):     σ = {:.4e} m = {:.4} pm",
        sigma_fe,
        sigma_fe * 1e12
    );
    println!("(Heavier atoms have narrower thermal lines)");

    // Line profiles
    println!("\n=== Line Profile Comparison (centered at 656.3 nm) ===");
    let center = 656.3;
    let sigma = 0.01; // nm
    let gamma = 0.005; // nm
    println!(
        "{:>10} {:>12} {:>12} {:>12}",
        "Δλ (nm)", "Gaussian", "Lorentzian", "Voigt"
    );
    for i in -10..=10 {
        let lam = center + i as f64 * 0.005;
        let g = spectral::gaussian_profile(lam, center, sigma);
        let l = spectral::lorentzian_profile(lam, center, gamma);
        let v = spectral::pseudo_voigt_profile(lam, center, sigma, gamma);
        println!(
            "{:>+10.4} {:>12.4} {:>12.4} {:>12.4}",
            lam - center,
            g,
            l,
            v
        );
    }
}
