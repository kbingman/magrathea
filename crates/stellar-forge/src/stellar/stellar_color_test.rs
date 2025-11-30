use rand::SeedableRng;
use rand_chacha::ChaChaRng;
use units::Mass;

use crate::stellar::generation::{neutron_star, solar_analog};
use crate::stellar::stellar_color::StellarColor;
use crate::stellar::stellar_objects::{BlackHole, StellarObject};

#[test]
fn stellar_color_from_temperature() {
    // Hot O-type star (~30000K) - should be blue-white
    let hot = StellarColor::from_temperature(30000.0);
    assert!(hot.b > hot.r, "Hot stars should have more blue than red");

    // Sun-like G-type star (~5800K) - should be yellowish-white
    let solar = StellarColor::from_temperature(5800.0);
    assert!(
        solar.r > solar.b,
        "Solar-type should have more red than blue"
    );
    assert!(
        solar.g > solar.b,
        "Solar-type should have more green than blue"
    );

    // Cool M-type star (~3000K) - should be reddish
    let cool = StellarColor::from_temperature(3000.0);
    assert!(cool.r > cool.b, "Cool stars should have more red than blue");
    assert!(
        cool.r > cool.g,
        "Cool stars should have more red than green"
    );

    // Temperature gradient: hotter stars should have higher blue/red ratio
    assert!(
        (hot.b as f64 / hot.r as f64) > (solar.b as f64 / solar.r as f64),
        "Blue/red ratio should increase with temperature"
    );
}

#[test]
fn stellar_color_temperature_clamping() {
    // Extreme temperatures should be clamped, not panic
    let very_cold = StellarColor::from_temperature(100.0);
    let very_hot = StellarColor::from_temperature(100000.0);

    // Very cold should be reddish (clamped to 1000K)
    assert!(very_cold.r > very_cold.b);

    // Very hot should be bluish-white (clamped to 40000K)
    assert!(very_hot.b >= very_hot.r);
}

#[test]
fn stellar_color_hex_roundtrip() {
    let color = StellarColor::new(255, 128, 64);
    assert_eq!(color.to_hex(), "#FF8040");

    let black = StellarColor::new(0, 0, 0);
    assert_eq!(black.to_hex(), "#000000");

    let white = StellarColor::new(255, 255, 255);
    assert_eq!(white.to_hex(), "#FFFFFF");

    // from_hex should parse back
    assert_eq!(StellarColor::from_hex("#FF8040").unwrap(), color);
    assert_eq!(StellarColor::from_hex("FF8040").unwrap(), color); // without #
    assert_eq!(StellarColor::from_hex("#000000").unwrap(), black);
}

#[test]
fn stellar_object_color() {
    let mut rng = ChaChaRng::seed_from_u64(42);

    // Main sequence star should use temperature
    let ms = StellarObject::MainSequence(solar_analog());
    let ms_color = ms.color();
    let direct_color = StellarColor::from_temperature(5800.0);
    assert_eq!(ms_color, direct_color);

    // Black hole without accretion should be nearly black
    let bh = StellarObject::BlackHole(BlackHole {
        mass: Mass::from_solar_masses(10.0),
        spin: 0.5,
        has_accretion: false,
    });
    let bh_color = bh.color();
    assert!(bh_color.r < 50 && bh_color.g < 50 && bh_color.b < 50);

    // Black hole with accretion should be brighter (orange-white)
    let bh_acc = StellarObject::BlackHole(BlackHole {
        mass: Mass::from_solar_masses(10.0),
        spin: 0.5,
        has_accretion: true,
    });
    let bh_acc_color = bh_acc.color();
    assert!(bh_acc_color.r > 200);

    // Neutron star should be blue-ish
    let ns = StellarObject::NeutronStar(neutron_star(&mut rng));
    let ns_color = ns.color();
    assert!(ns_color.b > ns_color.r);
}
