//! Basic star creation, classification, and property computation.

use tara::classification::{self, LuminosityClass};
use tara::evolution;
use tara::luminosity;
use tara::star::Star;

fn main() {
    // Create the Sun using convenience constructor
    let sun = Star::sun().unwrap();
    println!("=== The Sun ===");
    println!(
        "Classification: {}",
        classification::format_classification(sun.temperature_k, sun.luminosity_class)
    );
    println!("Temperature: {} K", sun.temperature_k);
    println!("Luminosity: {:.2} L_sun", sun.luminosity_solar);
    println!(
        "Absolute magnitude: {:.2}",
        luminosity::absolute_bolometric_magnitude(sun.luminosity_solar)
    );
    println!(
        "Main-sequence lifetime: {:.2e} yr",
        evolution::main_sequence_lifetime(sun.mass_solar)
    );
    println!(
        "Evolutionary phase: {:?}",
        evolution::evolutionary_phase(sun.mass_solar, sun.age_years)
    );
    println!(
        "Remnant type: {:?}",
        evolution::remnant_type(sun.mass_solar)
    );

    // Create a massive O-type star with the builder
    println!("\n=== Rigel (approximation) ===");
    let rigel = Star::builder(21.0, 78.9, 12_100.0, 120_000.0, 8e6)
        .luminosity_class(LuminosityClass::Ia)
        .metallicity(-0.06)
        .build()
        .unwrap();
    println!(
        "Classification: {}",
        classification::format_classification(rigel.temperature_k, rigel.luminosity_class)
    );
    println!("Mass: {:.1} M_sun", rigel.mass_solar);
    println!(
        "Main-sequence lifetime: {:.2e} yr",
        evolution::main_sequence_lifetime(rigel.mass_solar)
    );
    println!(
        "Remnant type: {:?}",
        evolution::remnant_type(rigel.mass_solar)
    );

    // Create a red dwarf
    println!("\n=== Proxima Centauri (approximation) ===");
    let proxima = Star::builder(0.122, 0.154, 3042.0, 0.0017, 4.85e9)
        .luminosity_class(LuminosityClass::V)
        .build()
        .unwrap();
    println!(
        "Classification: {}",
        classification::format_classification(proxima.temperature_k, proxima.luminosity_class)
    );
    println!(
        "Main-sequence lifetime: {:.2e} yr",
        evolution::main_sequence_lifetime(proxima.mass_solar)
    );
    println!("(Much longer than the age of the universe!)");

    // Demonstrate spectral class auto-derivation
    println!("\n=== Auto-classification ===");
    for temp in [
        35_000.0, 15_000.0, 8_000.0, 6_500.0, 5_772.0, 4_500.0, 3_000.0,
    ] {
        let class = classification::spectral_class_from_temperature(temp);
        let sub = classification::spectral_subclass(temp);
        println!("{temp:>8.0} K → {class}{sub}");
    }

    // Distance and apparent magnitude
    println!("\n=== Distance ===");
    let abs_mag = luminosity::absolute_bolometric_magnitude(sun.luminosity_solar);
    for dist in [1.0, 10.0, 100.0, 1000.0] {
        let app = luminosity::apparent_magnitude(abs_mag, dist);
        println!("Sun at {dist:>6.0} pc: apparent mag = {app:+.2}");
    }
}
