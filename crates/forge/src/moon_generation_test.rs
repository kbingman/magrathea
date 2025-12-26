use rand::SeedableRng;
use rand_chacha::ChaChaRng;

use crate::moon_generation::{
    GAS_GIANT_MOON_DIST, hill_radius, planet_category, roche_limit, sample_moon_count,
    tidal_heat_flux,
};

#[test]
fn test_hill_radius() {
    // Earth-Sun system: Hill radius should be ~1.5 million km
    let earth_mass_kg = 5.972e24;
    let sun_mass_kg = 1.989e30;
    let au_in_m = 1.496e11;

    let r_hill = hill_radius(au_in_m, earth_mass_kg, sun_mass_kg);

    // Expected: ~1.5e9 m (1.5 million km)
    assert!(r_hill > 1.0e9 && r_hill < 2.0e9, "Hill radius: {}", r_hill);
}

#[test]
fn test_roche_limit() {
    // Earth: Roche limit for Moon-density body
    // R_Roche = 2.44 × R_planet × (ρ_planet / ρ_moon)^(1/3)
    let earth_radius_m = 6.371e6;
    let earth_density = 5515.0; // kg/m³
    let moon_density = 3340.0; // kg/m³

    let r_roche = roche_limit(earth_radius_m, earth_density, moon_density);

    // Expected: ~18,000 km for Moon-density body
    // (Moon orbits at 384,000 km, well outside this)
    assert!(
        r_roche > 15.0e6 && r_roche < 22.0e6,
        "Roche limit: {}",
        r_roche
    );
}

#[test]
fn test_moon_count_sampling() {
    let mut rng = ChaChaRng::seed_from_u64(42);

    // Gas giants should usually have moons
    let mut has_moons = 0;
    for _ in 0..100 {
        let count = sample_moon_count(&mut rng, GAS_GIANT_MOON_DIST);
        if count > 0 {
            has_moons += 1;
        }
    }

    assert!(
        has_moons > 90,
        "Gas giants should usually have moons: {}/100",
        has_moons
    );
}

#[test]
fn test_planet_category() {
    assert_eq!(planet_category(0.5), 0); // Terrestrial
    assert_eq!(planet_category(1.0), 0); // Terrestrial (Earth)
    assert_eq!(planet_category(5.0), 1); // Super-Earth
    assert_eq!(planet_category(17.0), 2); // Ice Giant (Neptune-mass)
    assert_eq!(planet_category(318.0), 3); // Gas Giant (Jupiter-mass)
}

#[test]
fn test_tidal_heat_flux_io() {
    // Io parameters - should give ~2 W/m²
    let jupiter_mass = 1.898e27; // kg
    let io_radius = 1.8216e6; // m
    let io_sma = 4.217e8; // m
    let io_eccentricity = 0.0041;

    let flux = tidal_heat_flux(jupiter_mass, io_radius, io_sma, io_eccentricity);

    assert!(
        (flux - 2.0).abs() < 0.1,
        "Io should have ~2 W/m² heat flux, got {}",
        flux
    );
}

#[test]
fn test_tidal_heat_flux_europa() {
    // Europa parameters - should give ~0.05-0.1 W/m²
    let jupiter_mass = 1.898e27; // kg
    let europa_radius = 1.5608e6; // m
    let europa_sma = 6.709e8; // m
    let europa_eccentricity = 0.009; // Forced by resonance

    let flux = tidal_heat_flux(jupiter_mass, europa_radius, europa_sma, europa_eccentricity);

    // Europa's heat flux is estimated at ~0.05-0.2 W/m²
    assert!(
        flux > 0.01 && flux < 0.5,
        "Europa should have ~0.05-0.2 W/m² heat flux, got {}",
        flux
    );
}

#[test]
fn test_tidal_heat_flux_scales_with_eccentricity() {
    let jupiter_mass = 1.898e27;
    let moon_radius = 1.5e6;
    let sma = 5.0e8;

    let flux_low_e = tidal_heat_flux(jupiter_mass, moon_radius, sma, 0.001);
    let flux_high_e = tidal_heat_flux(jupiter_mass, moon_radius, sma, 0.01);

    // 10x eccentricity should give 100x heat flux (e² scaling)
    let ratio = flux_high_e / flux_low_e;
    assert!(
        (ratio - 100.0).abs() < 1.0,
        "Heat flux should scale as e², ratio = {}",
        ratio
    );
}
