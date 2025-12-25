use crate::{
    composition::Composition,
    planet_class::PlanetClass,
    planet_type::{PlanetType, evaporation_threshold},
};

#[test]
fn test_serialization_includes_name_and_description() {
    let gas_giant = PlanetType::GasGiant { has_rings: true };
    let json = serde_json::to_string(&gas_giant).unwrap();

    assert!(json.contains(r#""type":"gasGiant""#));
    assert!(json.contains(r#""name":"Gas Giant""#));
    assert!(json.contains(r#""description":"Gas giant: Jupiter-like, banded clouds""#));
    assert!(json.contains(r#""hasRings":true"#));

    let terran = PlanetType::Terran {
        ocean_fraction: 0.7,
        tectonically_active: true,
    };
    let json = serde_json::to_string(&terran).unwrap();

    assert!(json.contains(r#""type":"terran""#));
    assert!(json.contains(r#""name":"Terran""#));
    assert!(json.contains(r#""description":"Terran: Earth-like with oceans and continents""#));
    assert!(json.contains(r#""oceanFraction""#));

    let desert = PlanetType::Desert;
    let json = serde_json::to_string(&desert).unwrap();

    assert!(json.contains(r#""type":"desert""#));
    assert!(json.contains(r#""name":"Desert""#));
    assert!(json.contains(r#""description":"Desert: thin atmosphere, no surface liquid""#));
}

#[test]
fn test_class_mapping() {
    assert_eq!(
        PlanetType::Terran {
            ocean_fraction: 0.7,
            tectonically_active: true
        }
        .class(),
        PlanetClass::Compact
    );
    assert_eq!(
        PlanetType::MiniNeptune {
            envelope_fraction: 0.05
        }
        .class(),
        PlanetClass::Transitional
    );
    assert_eq!(
        PlanetType::IceGiant { has_rings: true }.class(),
        PlanetClass::Volatile
    );
    assert_eq!(
        PlanetType::GasGiant { has_rings: true }.class(),
        PlanetClass::Giant
    );
}

#[test]
fn test_habitability() {
    assert!(
        PlanetType::Terran {
            ocean_fraction: 0.7,
            tectonically_active: true
        }
        .potentially_habitable()
    );
    assert!(PlanetType::Hycean.potentially_habitable());
    assert!(
        !PlanetType::Lava {
            tidally_locked: true,
            surface_temp_k: 2000.0
        }
        .potentially_habitable()
    );
    assert!(!PlanetType::GasGiant { has_rings: true }.potentially_habitable());
}

// =============================================================================
// Photoevaporation tests
// =============================================================================

#[test]
fn test_evaporation_threshold_scales_with_mass() {
    // More massive planets can resist more flux
    let threshold_2m = evaporation_threshold(2.0);
    let threshold_5m = evaporation_threshold(5.0);
    let threshold_10m = evaporation_threshold(10.0);

    assert!(
        threshold_5m > threshold_2m,
        "5 M⊕ should resist more flux than 2 M⊕"
    );
    assert!(
        threshold_10m > threshold_5m,
        "10 M⊕ should resist more flux than 5 M⊕"
    );

    // Check approximate values (base is 20 F⊕ at 2 M⊕, exponent 1.8)
    // At 2 M⊕: 20 * (2/2)^1.8 = 20 F⊕
    assert!(
        (threshold_2m - 20.0).abs() < 0.1,
        "Threshold at 2 M⊕ should be ~20 F⊕"
    );

    // At 5 M⊕: 20 * (5/2)^1.8 ≈ 20 * 4.9 ≈ 98 F⊕
    assert!(
        threshold_5m > 80.0 && threshold_5m < 120.0,
        "Threshold at 5 M⊕ should be ~98 F⊕, got {}",
        threshold_5m
    );
}

#[test]
fn test_transitional_stripped_at_high_flux() {
    // Composition with H/He envelope
    let with_envelope = Composition::new(0.3, 0.5, 0.1, 0.1); // 10% H/He

    // Low-mass transitional planet (2 M⊕) at high flux should be stripped
    let high_flux = 100.0; // Well above threshold for 2 M⊕
    let low_mass = 2.0;
    let planet_type = PlanetType::from_environment(
        PlanetClass::Transitional,
        &with_envelope,
        300.0, // temperature
        high_flux,
        low_mass,
        1.0, // stellar mass
        0.1, // close-in
    );

    assert!(
        matches!(planet_type, PlanetType::SuperTerran),
        "Close-in low-mass transitional should be stripped to SuperTerran, got {:?}",
        planet_type
    );
}

#[test]
fn test_transitional_retains_envelope_at_low_flux() {
    // Composition with H/He envelope
    let with_envelope = Composition::new(0.3, 0.5, 0.1, 0.1); // 10% H/He

    // Same mass but at low flux should retain envelope
    let low_flux = 5.0; // Below threshold
    let planet_type = PlanetType::from_environment(
        PlanetClass::Transitional,
        &with_envelope,
        200.0,
        low_flux,
        2.0,
        1.0,
        1.0,
    );

    assert!(
        matches!(planet_type, PlanetType::MiniNeptune { .. }),
        "Low-flux transitional should retain envelope as MiniNeptune, got {:?}",
        planet_type
    );
}

#[test]
fn test_massive_transitional_resists_stripping() {
    // Composition with H/He envelope
    let with_envelope = Composition::new(0.3, 0.5, 0.1, 0.1); // 10% H/He

    // 5 M⊕ at moderate flux - should retain envelope because mass is higher
    let moderate_flux = 50.0; // Above 2 M⊕ threshold but below 5 M⊕ threshold
    let planet_type = PlanetType::from_environment(
        PlanetClass::Transitional,
        &with_envelope,
        250.0,
        moderate_flux,
        5.0,
        1.0,
        0.5,
    );

    assert!(
        matches!(planet_type, PlanetType::MiniNeptune { .. }),
        "Massive transitional at moderate flux should retain envelope, got {:?}",
        planet_type
    );
}

#[test]
fn test_volatile_becomes_chthonian_under_extreme_flux() {
    // Low-mass volatile (just above transitional threshold) under extreme irradiation
    // At 6 M⊕: threshold = 20 * (6/2)^1.8 ≈ 128 F⊕
    // Need > 128 * 3 = 384 F⊕ for Chthonian, use 500 F⊕
    let planet_type = PlanetType::from_environment(
        PlanetClass::Volatile,
        &Composition::ice_giant(),
        1500.0, // Hot
        500.0,  // Extreme flux
        6.0,    // Low-mass volatile (just above 5 M⊕ transitional boundary)
        1.0,
        0.05, // Very close-in
    );

    assert!(
        matches!(planet_type, PlanetType::Chthonian),
        "Low-mass volatile under extreme flux should become Chthonian, got {:?}",
        planet_type
    );
}
