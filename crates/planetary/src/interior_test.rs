//! Tests for interior thermal state calculations

use crate::interior::{DifferentiationState, HeatSource, ThermalState};

// ========== ThermalState Tests ==========

#[test]
fn test_earth_thermal_state() {
    let thermal = ThermalState::calculate(1.0, 4.5, 0.01);

    // Earth's core ~5200 K
    assert!(
        thermal.core_temperature > 4000.0 && thermal.core_temperature < 7000.0,
        "Earth core temp should be ~5000K, got {}",
        thermal.core_temperature
    );

    // Earth heat flow ~0.087 W/m²
    assert!(
        thermal.heat_flow > 0.05 && thermal.heat_flow < 0.15,
        "Earth heat flow should be ~0.087 W/m², got {}",
        thermal.heat_flow
    );

    // Radiogenic should dominate
    assert!(
        thermal.radiogenic_power > thermal.tidal_power,
        "Earth should be radiogenically heated"
    );
}

#[test]
fn test_mars_thermal_state() {
    let thermal = ThermalState::calculate(0.107, 4.5, 0.01);

    // Mars is cooler
    assert!(
        thermal.core_temperature < 4000.0,
        "Mars core should be cooler than Earth"
    );

    // Lower heat flow
    assert!(
        thermal.heat_flow < 0.05,
        "Mars heat flow should be low, got {}",
        thermal.heat_flow
    );
}

#[test]
fn test_io_like_tidal_heating() {
    // Io: 0.015 M⊕, strong tidal heating (factor ~100-150)
    let thermal = ThermalState::calculate(0.015, 4.5, 150.0);

    // Tidal should dominate
    assert!(
        thermal.tidal_power > thermal.radiogenic_power,
        "Io should be tidally heated: tidal={}, radiogenic={}",
        thermal.tidal_power,
        thermal.radiogenic_power
    );

    // High heat flow despite small mass
    assert!(
        thermal.heat_flow > 1.0,
        "Io should have high heat flow, got {}",
        thermal.heat_flow
    );
}

#[test]
fn test_young_planet_hotter() {
    let young = ThermalState::calculate(1.0, 0.5, 0.01);
    let old = ThermalState::calculate(1.0, 4.5, 0.01);

    // Young planets retain more heat
    assert!(
        young.core_temperature > old.core_temperature,
        "Young planet should be hotter"
    );

    assert!(
        young.heat_flow > old.heat_flow,
        "Young planet should have higher heat flow"
    );
}

#[test]
fn test_massive_planet_hotter() {
    let small = ThermalState::calculate(0.5, 4.5, 0.01);
    let large = ThermalState::calculate(2.0, 4.5, 0.01);

    // Larger planets retain more heat
    assert!(
        large.core_temperature > small.core_temperature,
        "Larger planet should be hotter"
    );
}

#[test]
fn test_heat_sources_ordering() {
    // Earth-like: radiogenic should come first
    let earth = ThermalState::calculate(1.0, 4.5, 0.01);
    assert!(
        matches!(earth.heat_sources.first(), Some(HeatSource::Radiogenic)),
        "Earth should have radiogenic as primary"
    );

    // Io-like: tidal should come first
    let io = ThermalState::calculate(0.015, 4.5, 150.0);
    assert!(
        matches!(io.heat_sources.first(), Some(HeatSource::Tidal)),
        "Io should have tidal as primary"
    );

    // Young planet: should include accretion
    let young = ThermalState::calculate(1.0, 0.05, 0.01);
    assert!(
        young.heat_sources.contains(&HeatSource::Accretion),
        "Young planet should include accretion heat"
    );
}

#[test]
fn test_is_tidally_heated() {
    let earth = ThermalState::calculate(1.0, 4.5, 0.01);
    assert!(!earth.is_tidally_heated());

    let io = ThermalState::calculate(0.015, 4.5, 150.0);
    assert!(io.is_tidally_heated());
}

#[test]
fn test_is_geologically_active() {
    let earth = ThermalState::calculate(1.0, 4.5, 0.01);
    assert!(earth.is_geologically_active());

    // Small old body should be inactive
    let moon = ThermalState::calculate(0.012, 4.5, 0.01);
    assert!(!moon.is_geologically_active());
}

// ========== DifferentiationState Tests ==========

#[test]
fn test_differentiation_undifferentiated() {
    let state = DifferentiationState::classify(0.05, 4.5);
    assert_eq!(state, DifferentiationState::Undifferentiated);
}

#[test]
fn test_differentiation_partial() {
    let state = DifferentiationState::classify(1.0, 0.05);
    assert_eq!(state, DifferentiationState::PartiallyDifferentiated);
}

#[test]
fn test_differentiation_full() {
    let state = DifferentiationState::classify(1.0, 4.5);
    assert_eq!(state, DifferentiationState::FullyDifferentiated);
}

#[test]
fn test_radiogenic_power_scales_with_mass() {
    let small = ThermalState::calculate(0.5, 4.5, 0.0);
    let large = ThermalState::calculate(2.0, 4.5, 0.0);

    assert!(
        large.radiogenic_power > small.radiogenic_power,
        "Larger planet should have more radiogenic power"
    );
}

#[test]
fn test_radiogenic_power_decays_with_age() {
    let young = ThermalState::calculate(1.0, 1.0, 0.0);
    let old = ThermalState::calculate(1.0, 8.0, 0.0);

    assert!(
        young.radiogenic_power > old.radiogenic_power,
        "Young planet should have more radiogenic power"
    );
}
