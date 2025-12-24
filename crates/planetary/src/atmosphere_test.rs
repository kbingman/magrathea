//! Tests for atmospheric classification

use crate::atmosphere::AtmosphereType;
use crate::planet_class::PlanetClass;

// ========== Basic Classification Tests ==========

#[test]
fn test_earth_atmosphere() {
    // Earth: 1 M⊕, 288K, 11.2 km/s escape velocity
    let atm = AtmosphereType::classify(1.0, 288.0, &PlanetClass::Compact, 11.2);
    assert_eq!(atm, AtmosphereType::NitrogenOxygen);
}

#[test]
fn test_mars_atmosphere() {
    // Mars: 0.107 M⊕, 210K, 5.0 km/s escape velocity
    let atm = AtmosphereType::classify(0.107, 210.0, &PlanetClass::Compact, 5.0);
    assert_eq!(atm, AtmosphereType::ThinCO2);
}

#[test]
fn test_venus_atmosphere() {
    // Venus: 0.815 M⊕, 735K (surface), 10.4 km/s escape velocity
    let atm = AtmosphereType::classify(0.815, 500.0, &PlanetClass::Compact, 10.4);
    assert_eq!(atm, AtmosphereType::ThickCO2);
}

#[test]
fn test_mercury_atmosphere() {
    // Mercury: 0.055 M⊕, 440K, 4.3 km/s escape velocity
    // Too hot and low escape velocity for significant atmosphere
    let atm = AtmosphereType::classify(0.055, 440.0, &PlanetClass::Compact, 4.3);
    assert_eq!(atm, AtmosphereType::Trace);
}

#[test]
fn test_moon_atmosphere() {
    // Moon: 0.012 M⊕, 250K, 2.4 km/s escape velocity
    let atm = AtmosphereType::classify(0.012, 250.0, &PlanetClass::Compact, 2.4);
    assert_eq!(atm, AtmosphereType::Trace);
}

// ========== Giant Planet Tests ==========

#[test]
fn test_jupiter_atmosphere() {
    // Jupiter: 318 M⊕, 165K, cold H/He atmosphere
    let atm = AtmosphereType::classify(318.0, 165.0, &PlanetClass::Giant, 59.5);
    assert_eq!(atm, AtmosphereType::HydrogenHelium);
}

#[test]
fn test_hot_jupiter_atmosphere() {
    // Hot Jupiter: ~318 M⊕, 1500K
    let atm = AtmosphereType::classify(318.0, 1500.0, &PlanetClass::Giant, 60.0);
    assert_eq!(atm, AtmosphereType::HydrogenMetalOxides);
}

#[test]
fn test_ultra_hot_jupiter_atmosphere() {
    // Ultra-hot Jupiter: ~318 M⊕, 2500K
    let atm = AtmosphereType::classify(318.0, 2500.0, &PlanetClass::Giant, 60.0);
    assert_eq!(atm, AtmosphereType::HydrogenIonizedMetal);
}

// ========== Ice Giant Tests ==========

#[test]
fn test_neptune_atmosphere() {
    // Neptune: 17.1 M⊕, 72K, H/CH4 atmosphere
    let atm = AtmosphereType::classify(17.1, 72.0, &PlanetClass::Volatile, 23.5);
    assert_eq!(atm, AtmosphereType::HydrogenMethane);
}

#[test]
fn test_warm_neptune_atmosphere() {
    // Warm Neptune: ~17 M⊕, 700K
    let atm = AtmosphereType::classify(17.0, 700.0, &PlanetClass::Volatile, 23.0);
    assert_eq!(atm, AtmosphereType::HydrogenSulfide);
}

// ========== Small Body Tests ==========

#[test]
fn test_titan_like_atmosphere() {
    // Titan: 0.023 M⊕, 94K, 2.6 km/s - thick methane/nitrogen
    let atm = AtmosphereType::classify(0.023, 94.0, &PlanetClass::Compact, 2.6);
    assert_eq!(atm, AtmosphereType::MethaneRich);
}

#[test]
fn test_pluto_like_atmosphere() {
    // Pluto: 0.0022 M⊕, 44K, 1.2 km/s - thin nitrogen
    let atm = AtmosphereType::classify(0.0022, 44.0, &PlanetClass::Compact, 1.2);
    assert_eq!(atm, AtmosphereType::ThinNitrogen);
}

// ========== Transitional Planet Tests ==========

#[test]
fn test_mini_neptune_atmosphere() {
    // Mini-Neptune: 4 M⊕, 300K, good retention
    let atm = AtmosphereType::classify(4.0, 300.0, &PlanetClass::Transitional, 12.0);
    // Could be Hycean or H/He depending on conditions
    assert!(atm.is_primordial());
}

#[test]
fn test_super_venus_atmosphere() {
    // Hot super-terrestrial: 3 M⊕, 700K
    let atm = AtmosphereType::classify(3.0, 700.0, &PlanetClass::Transitional, 15.0);
    assert_eq!(atm, AtmosphereType::WaterVapor);
}

// ========== Geological Context Tests ==========

#[test]
fn test_io_like_atmosphere() {
    // Io: 0.015 M⊕, volcanically active, thin SO2
    let atm = AtmosphereType::classify_with_geology(
        0.015,
        130.0,
        &PlanetClass::Compact,
        2.6,
        true, // volcanically active
        None,
    );
    assert_eq!(atm, AtmosphereType::ThinSulfur);
}

#[test]
fn test_europa_like_atmosphere() {
    // Europa: 0.008 M⊕, icy with radiolysis O2
    let atm = AtmosphereType::classify_with_geology(
        0.008,
        100.0,
        &PlanetClass::Compact,
        2.0,
        false,
        Some(0.6), // high water fraction
    );
    assert_eq!(atm, AtmosphereType::ThinOxygen);
}

// ========== Greenhouse Effect Tests ==========

#[test]
fn test_no_greenhouse_effect() {
    assert_eq!(AtmosphereType::None.greenhouse_effect(), 0.0);
    assert_eq!(AtmosphereType::Trace.greenhouse_effect(), 0.0);
}

#[test]
fn test_earth_greenhouse_effect() {
    assert!((AtmosphereType::NitrogenOxygen.greenhouse_effect() - 33.0).abs() < 1.0);
}

#[test]
fn test_venus_greenhouse_effect() {
    assert!(AtmosphereType::ThickCO2.greenhouse_effect() > 400.0);
}

#[test]
fn test_mars_greenhouse_effect() {
    assert!((AtmosphereType::ThinCO2.greenhouse_effect() - 5.0).abs() < 1.0);
}

// ========== Property Tests ==========

#[test]
fn test_primordial_atmospheres() {
    assert!(AtmosphereType::HydrogenHelium.is_primordial());
    assert!(AtmosphereType::HydrogenMethane.is_primordial());
    assert!(AtmosphereType::Hydrogen.is_primordial());
    assert!(!AtmosphereType::NitrogenOxygen.is_primordial());
    assert!(!AtmosphereType::ThickCO2.is_primordial());
}

#[test]
fn test_potentially_habitable() {
    assert!(AtmosphereType::NitrogenOxygen.potentially_habitable());
    assert!(!AtmosphereType::None.potentially_habitable());
    assert!(!AtmosphereType::ThickCO2.potentially_habitable());
}

#[test]
fn test_display() {
    assert_eq!(format!("{}", AtmosphereType::NitrogenOxygen), "N₂/O₂");
    assert_eq!(format!("{}", AtmosphereType::HydrogenHelium), "H₂/He");
    assert_eq!(format!("{}", AtmosphereType::ThickCO2), "Thick CO₂");
}
