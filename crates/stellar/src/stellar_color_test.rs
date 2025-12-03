use rand::SeedableRng;
use rand_chacha::ChaChaRng;

use crate::generation::{black_hole, neutron_star, solar_analog};
use crate::stellar_color::StellarColor;

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

    // Main sequence star should have color based on temperature
    let ms = solar_analog();
    let expected_color = StellarColor::from_temperature(5800.0);
    assert_eq!(ms.color, expected_color);

    // Black hole color depends on accretion state
    let bh = black_hole(&mut rng, 25.0, 0.0);
    if bh.has_accretion {
        // Orange-white accretion glow
        assert_eq!(bh.color, StellarColor::new(255, 200, 150));
    } else {
        // Nearly black
        assert_eq!(bh.color, StellarColor::new(20, 20, 30));
    }

    // Neutron star should be blue-ish
    let ns = neutron_star(&mut rng);
    assert!(ns.color.b > ns.color.r);
}
