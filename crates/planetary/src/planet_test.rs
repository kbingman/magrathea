use rand_chacha::ChaChaRng;

use crate::{
    planet::{
        calculate_equilibrium_temp, earth_analog, jupiter_analog, neptune_analog, planet_letter,
        radius_with_envelope,
    },
    planet_class::PlanetClass,
    planet_type::PlanetType,
};

use rand::SeedableRng;

#[test]
fn test_earth() {
    let e = earth_analog();
    assert_eq!(e.class, PlanetClass::Compact);
    assert!(matches!(
        e.planet_type,
        PlanetType::Terran { .. } | PlanetType::Desert
    ));

    // Density ~5.5 g/cm³
    let d = e.density();
    assert!(d > 5.0 && d < 6.0, "Earth density: {}", d);
}

#[test]
fn test_jupiter() {
    let j = jupiter_analog();
    assert_eq!(j.class, PlanetClass::Giant);
    assert!(matches!(j.planet_type, PlanetType::GasGiant { .. }));

    // Density ~1.3 g/cm³
    let d = j.density();
    assert!(d > 1.0 && d < 2.0, "Jupiter density: {}", d);
}

#[test]
fn test_neptune() {
    let n = neptune_analog();
    assert_eq!(n.class, PlanetClass::Volatile);
    assert!(matches!(n.planet_type, PlanetType::IceGiant { .. }));
}

#[test]
fn test_envelope_radius() {
    let mut rng = ChaChaRng::seed_from_u64(42);

    // Sub-Neptune: 5 M⊕ with 10% envelope should be ~2-3 R⊕
    let r = radius_with_envelope(5.0, 0.10, &mut rng);
    assert!(r > 1.5 && r < 4.0, "Sub-Neptune radius: {}", r);

    // Same core mass without envelope should be smaller
    let r_bare = PlanetClass::Transitional.radius_from_mass(4.5, false, &mut rng);
    assert!(r_bare < r, "Bare core should be smaller than with envelope");
}

#[test]
fn test_equilibrium_temp() {
    let t_earth = calculate_equilibrium_temp(1.0, 1.0);
    assert!(t_earth > 250.0 && t_earth < 290.0);

    let t_mercury = calculate_equilibrium_temp(0.39, 1.0);
    assert!(t_mercury > t_earth);

    let t_jupiter = calculate_equilibrium_temp(5.2, 1.0);
    assert!(t_jupiter < t_earth);
}

#[test]
fn test_planet_letter() {
    // Standard letters: b through z
    assert_eq!(planet_letter(0), "b"); // innermost
    assert_eq!(planet_letter(1), "c");
    assert_eq!(planet_letter(2), "d");
    assert_eq!(planet_letter(24), "z");

    // Extended: aa, ab, ac... for systems with >25 planets
    assert_eq!(planet_letter(25), "aa");
    assert_eq!(planet_letter(26), "ab");
    assert_eq!(planet_letter(50), "az");
    assert_eq!(planet_letter(51), "ba");
}
