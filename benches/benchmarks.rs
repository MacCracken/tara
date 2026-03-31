use criterion::{Criterion, black_box, criterion_group, criterion_main};
use tara::bridge;
use tara::classification::LuminosityClass;
use tara::star::{SpectralClass, Star};
use tara::{atmosphere, evolution, luminosity, nucleosynthesis, spectral, sse};

fn bench_star_construction(c: &mut Criterion) {
    c.bench_function("star_new", |b| {
        b.iter(|| {
            Star::new(
                black_box(1.0),
                black_box(1.0),
                black_box(5772.0),
                black_box(1.0),
                black_box(4.6e9),
                black_box(SpectralClass::G),
            )
        })
    });

    c.bench_function("star_builder", |b| {
        b.iter(|| {
            Star::builder(
                black_box(1.0),
                black_box(1.0),
                black_box(5772.0),
                black_box(1.0),
                black_box(4.6e9),
            )
            .luminosity_class(black_box(LuminosityClass::V))
            .metallicity(black_box(0.0))
            .build()
        })
    });
}

fn bench_bridge_conversions(c: &mut Criterion) {
    c.bench_function("stellar_mass_to_mu", |b| {
        b.iter(|| bridge::stellar_mass_solar_to_mu(black_box(1.0)))
    });

    c.bench_function("temperature_to_peak_wavelength", |b| {
        b.iter(|| bridge::temperature_to_peak_wavelength_nm(black_box(5772.0)))
    });

    c.bench_function("surface_gravity_to_limb_darkening", |b| {
        b.iter(|| bridge::surface_gravity_to_limb_darkening(black_box(4.44)))
    });

    c.bench_function("composition_to_lifetime_scale", |b| {
        b.iter(|| bridge::composition_to_lifetime_scale(black_box(0.74)))
    });

    c.bench_function("core_temp_to_cno_fraction", |b| {
        b.iter(|| bridge::core_temp_to_cno_fraction(black_box(17e6)))
    });

    c.bench_function("luminosity_class_to_abs_magnitude", |b| {
        b.iter(|| bridge::luminosity_class_to_abs_magnitude(black_box(LuminosityClass::V)))
    });
}

fn bench_luminosity(c: &mut Criterion) {
    c.bench_function("luminosity_from_radius_temp", |b| {
        b.iter(|| luminosity::luminosity_from_radius_temp(black_box(6.957e8), black_box(5772.0)))
    });

    c.bench_function("absolute_bolometric_magnitude", |b| {
        b.iter(|| luminosity::absolute_bolometric_magnitude(black_box(1.0)))
    });

    c.bench_function("distance_modulus", |b| {
        b.iter(|| luminosity::distance_modulus(black_box(100.0)))
    });

    c.bench_function("bolometric_correction", |b| {
        b.iter(|| luminosity::bolometric_correction(black_box(5772.0)))
    });
}

fn bench_atmosphere(c: &mut Criterion) {
    c.bench_function("effective_temperature", |b| {
        b.iter(|| atmosphere::effective_temperature(black_box(3.828e26), black_box(6.957e8)))
    });

    c.bench_function("log_surface_gravity", |b| {
        b.iter(|| atmosphere::log_surface_gravity(black_box(1.0), black_box(1.0)))
    });

    c.bench_function("grey_atmosphere_temperature", |b| {
        b.iter(|| atmosphere::grey_atmosphere_temperature(black_box(5772.0), black_box(0.667)))
    });

    c.bench_function("kramers_opacity", |b| {
        b.iter(|| {
            atmosphere::kramers_opacity(
                black_box(1.0),
                black_box(1e6),
                black_box(0.74),
                black_box(0.02),
            )
        })
    });
}

fn bench_evolution(c: &mut Criterion) {
    c.bench_function("mass_luminosity", |b| {
        b.iter(|| evolution::mass_luminosity(black_box(1.0)))
    });

    c.bench_function("main_sequence_lifetime", |b| {
        b.iter(|| evolution::main_sequence_lifetime(black_box(1.0)))
    });

    c.bench_function("evolutionary_phase", |b| {
        b.iter(|| evolution::evolutionary_phase(black_box(1.0), black_box(4.6e9)))
    });
}

fn bench_spectral(c: &mut Criterion) {
    c.bench_function("planck_radiance_nm", |b| {
        b.iter(|| spectral::planck_radiance_nm(black_box(500.0), black_box(5772.0)))
    });

    c.bench_function("doppler_shift", |b| {
        b.iter(|| spectral::doppler_shift(black_box(500.0), black_box(1e6)))
    });

    c.bench_function("gaussian_profile", |b| {
        b.iter(|| spectral::gaussian_profile(black_box(500.0), black_box(500.5), black_box(0.1)))
    });

    c.bench_function("pseudo_voigt_profile", |b| {
        b.iter(|| {
            spectral::pseudo_voigt_profile(
                black_box(500.0),
                black_box(500.5),
                black_box(0.1),
                black_box(0.05),
            )
        })
    });
}

fn bench_nucleosynthesis(c: &mut Criterion) {
    c.bench_function("pp_chain_rate", |b| {
        b.iter(|| {
            nucleosynthesis::pp_chain_rate(black_box(150.0), black_box(0.34), black_box(15.7e6))
        })
    });

    c.bench_function("cno_cycle_rate", |b| {
        b.iter(|| {
            nucleosynthesis::cno_cycle_rate(
                black_box(100.0),
                black_box(0.70),
                black_box(0.01),
                black_box(25e6),
            )
        })
    });

    c.bench_function("triple_alpha_rate", |b| {
        b.iter(|| {
            nucleosynthesis::triple_alpha_rate(black_box(1e4), black_box(0.98), black_box(200e6))
        })
    });
}

fn bench_sse(c: &mut Criterion) {
    c.bench_function("sse_zams_luminosity", |b| {
        b.iter(|| sse::zams_luminosity(black_box(1.0), black_box(0.02)))
    });

    c.bench_function("sse_zams_radius", |b| {
        b.iter(|| sse::zams_radius(black_box(1.0), black_box(0.02)))
    });

    c.bench_function("sse_zams_properties", |b| {
        b.iter(|| sse::zams_properties(black_box(1.0), black_box(0.02)))
    });

    c.bench_function("sse_t_bgb", |b| {
        b.iter(|| sse::t_bgb(black_box(1.0), black_box(0.02)))
    });

    c.bench_function("sse_ms_lifetime", |b| {
        b.iter(|| sse::ms_lifetime(black_box(1.0), black_box(0.02)))
    });

    c.bench_function("sse_tms_luminosity", |b| {
        b.iter(|| sse::tms_luminosity(black_box(1.0), black_box(0.02)))
    });

    c.bench_function("sse_tms_radius", |b| {
        b.iter(|| sse::tms_radius(black_box(1.0), black_box(0.02)))
    });

    c.bench_function("sse_ms_properties", |b| {
        b.iter(|| sse::ms_properties(black_box(1.0), black_box(0.02), black_box(4.6e9)))
    });

    c.bench_function("sse_lbgb", |b| {
        b.iter(|| sse::lbgb(black_box(1.0), black_box(0.02)))
    });

    c.bench_function("sse_hg_luminosity", |b| {
        b.iter(|| sse::hg_luminosity(black_box(1.0), black_box(0.02), black_box(0.5)))
    });

    c.bench_function("sse_hg_radius", |b| {
        b.iter(|| sse::hg_radius(black_box(1.0), black_box(0.02), black_box(0.5)))
    });

    c.bench_function("sse_rgb_luminosity", |b| {
        b.iter(|| sse::rgb_luminosity(black_box(1.0), black_box(0.02), black_box(0.5)))
    });

    c.bench_function("sse_rgb_radius", |b| {
        b.iter(|| sse::rgb_radius(black_box(1.0), black_box(100.0), black_box(0.02)))
    });

    c.bench_function("sse_evolve", |b| {
        b.iter(|| sse::evolve(black_box(1.0), black_box(0.02), black_box(4.6e9)))
    });
}

criterion_group!(
    benches,
    bench_star_construction,
    bench_bridge_conversions,
    bench_luminosity,
    bench_atmosphere,
    bench_evolution,
    bench_spectral,
    bench_nucleosynthesis,
    bench_sse,
);
criterion_main!(benches);
