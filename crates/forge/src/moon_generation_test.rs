use rand::SeedableRng;
use rand_chacha::ChaChaRng;

use crate::moon_generation::{
    GAS_GIANT_MOON_DIST, hill_radius, planet_category, roche_limit, sample_moon_count,
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
