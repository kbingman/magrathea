//! Tests for temperature calculations

use crate::atmosphere::AtmosphereType;
use crate::planet_class::PlanetClass;
use crate::temperature::{Temperature, TemperatureClass};

// ========== Equilibrium Temperature Tests ==========

#[test]
fn test_equilibrium_temperature_earth() {
    // Earth: 1 AU, solar luminosity, albedo ~0.3
    let temp = Temperature::equilibrium_temperature(1.0, 1.0, 0.3);
    // Earth's equilibrium temperature is about 255 K
    assert!(
        (temp - 255.0).abs() < 5.0,
        "Earth equilibrium temp should be ~255K, got {}",
        temp
    );
}

#[test]
fn test_equilibrium_temperature_mars() {
    // Mars: 1.52 AU, albedo ~0.25
    let temp = Temperature::equilibrium_temperature(1.52, 1.0, 0.25);
    // Mars equilibrium temperature is about 210 K
    assert!(
        (temp - 210.0).abs() < 10.0,
        "Mars equilibrium temp should be ~210K, got {}",
        temp
    );
}

#[test]
fn test_equilibrium_temperature_venus() {
    // Venus: 0.72 AU, albedo ~0.75 (very reflective clouds)
    let temp = Temperature::equilibrium_temperature(0.72, 1.0, 0.75);
    // Venus equilibrium temperature is about 230 K (before greenhouse)
    assert!(
        (temp - 230.0).abs() < 15.0,
        "Venus equilibrium temp should be ~230K, got {}",
        temp
    );
}

#[test]
fn test_equilibrium_temperature_jupiter() {
    // Jupiter: 5.2 AU, albedo ~0.5
    let temp = Temperature::equilibrium_temperature(5.2, 1.0, 0.5);
    // Jupiter equilibrium temperature is about 110 K
    assert!(
        (temp - 110.0).abs() < 15.0,
        "Jupiter equilibrium temp should be ~110K, got {}",
        temp
    );
}

#[test]
fn test_equilibrium_temperature_hot_jupiter() {
    // Hot Jupiter: 0.05 AU, albedo ~0.1
    let temp = Temperature::equilibrium_temperature(0.05, 1.0, 0.1);
    // Should be very hot, over 1200 K
    assert!(temp > 1200.0, "Hot Jupiter should be >1200K, got {}", temp);
    assert!(temp < 1500.0, "Hot Jupiter should be <1500K, got {}", temp);
}

#[test]
fn test_equilibrium_temperature_scales_with_distance() {
    let near = Temperature::equilibrium_temperature(0.5, 1.0, 0.3);
    let far = Temperature::equilibrium_temperature(2.0, 1.0, 0.3);

    // Temperature decreases with distance (proportional to 1/sqrt(distance))
    assert!(near > far, "Closer planet should be hotter");
    assert!(near > far * 1.5, "Temperature should scale significantly");
}

#[test]
fn test_equilibrium_temperature_scales_with_luminosity() {
    let dim_star = Temperature::equilibrium_temperature(1.0, 0.5, 0.3);
    let bright_star = Temperature::equilibrium_temperature(1.0, 2.0, 0.3);

    // Temperature increases with stellar luminosity
    assert!(bright_star > dim_star, "Brighter star should be hotter");
}

#[test]
fn test_equilibrium_temperature_albedo_effect() {
    let dark = Temperature::equilibrium_temperature(1.0, 1.0, 0.1); // Low albedo (dark)
    let bright = Temperature::equilibrium_temperature(1.0, 1.0, 0.7); // High albedo (bright)

    // Dark surfaces absorb more, get hotter
    assert!(
        dark > bright,
        "Dark planet should be hotter than bright planet"
    );
}

// ========== Greenhouse Effect Tests ==========

#[test]
fn test_greenhouse_effect_none() {
    let temp = Temperature::calculate(1.0, 255.0, &PlanetClass::Compact, &AtmosphereType::None);

    // No atmosphere = no greenhouse effect (only minimal internal heating)
    assert!(
        (temp.effective - temp.equilibrium).abs() < 1.0,
        "No atmosphere should have minimal temperature increase"
    );
}

#[test]
fn test_greenhouse_effect_earth_like() {
    let temp = Temperature::calculate(
        1.0,
        255.0,
        &PlanetClass::Compact,
        &AtmosphereType::NitrogenOxygen,
    );

    // Earth's greenhouse effect is about 33 K
    assert!(
        (temp.effective - temp.equilibrium - 33.0).abs() < 2.0,
        "Earth-like atmosphere should have ~33K greenhouse effect, got {}K increase",
        temp.effective - temp.equilibrium
    );
}

#[test]
fn test_greenhouse_effect_venus_like() {
    let temp = Temperature::calculate(
        0.815,
        230.0,
        &PlanetClass::Compact,
        &AtmosphereType::ThickCO2,
    );

    // Venus has extreme greenhouse effect (~500 K)
    assert!(
        temp.effective > 700.0,
        "Venus-like atmosphere should have extreme greenhouse, got {}K",
        temp.effective
    );
}

#[test]
fn test_greenhouse_effect_mars_like() {
    let temp = Temperature::calculate(
        0.107,
        210.0,
        &PlanetClass::Compact,
        &AtmosphereType::ThinCO2,
    );

    // Mars has minimal greenhouse effect (~5 K)
    assert!(
        (temp.effective - temp.equilibrium - 5.0).abs() < 1.0,
        "Thin CO2 atmosphere should have ~5K greenhouse effect"
    );
}

#[test]
fn test_greenhouse_effect_water_vapor() {
    let temp = Temperature::calculate(
        1.0,
        300.0,
        &PlanetClass::Compact,
        &AtmosphereType::WaterVapor,
    );

    // Water vapor creates strong greenhouse effect
    assert!(
        temp.effective > temp.equilibrium + 50.0,
        "Water vapor should create strong greenhouse effect"
    );
}

// ========== Internal Heating Tests ==========

#[test]
fn test_internal_heating_terrestrial_minimal() {
    let temp = Temperature::calculate(1.0, 255.0, &PlanetClass::Compact, &AtmosphereType::None);

    // Terrestrial planets have negligible internal heating
    assert!(
        (temp.effective - temp.equilibrium).abs() < 1.0,
        "Terrestrial planets should have minimal internal heating"
    );
}

#[test]
fn test_internal_heating_gas_giant() {
    let temp = Temperature::calculate(
        318.0,
        110.0,
        &PlanetClass::Giant,
        &AtmosphereType::HydrogenHelium,
    );

    // Gas giants have significant internal heating
    assert!(
        temp.effective > temp.equilibrium + 50.0,
        "Gas giants should have significant internal heating, got effective {}K from equilibrium {}K",
        temp.effective,
        temp.equilibrium
    );
}

#[test]
fn test_internal_heating_ice_giant() {
    let temp = Temperature::calculate(
        17.1,
        72.0,
        &PlanetClass::Volatile,
        &AtmosphereType::HydrogenMethane,
    );

    // Ice giants have some internal heating but less than gas giants
    assert!(
        temp.effective > temp.equilibrium,
        "Ice giants should have some internal heating"
    );
}

#[test]
fn test_internal_heating_scales_with_mass() {
    let small = Temperature::calculate(
        50.0,
        200.0,
        &PlanetClass::Volatile,
        &AtmosphereType::HydrogenMethane,
    );

    let large = Temperature::calculate(
        150.0,
        200.0,
        &PlanetClass::Giant,
        &AtmosphereType::HydrogenHelium,
    );

    // More massive planets have more internal heating
    assert!(
        large.effective > small.effective,
        "More massive giants should have more internal heating"
    );
}

// ========== Combined Effects Tests ==========

#[test]
fn test_effective_temperature_earth() {
    let temp = Temperature::calculate(
        1.0,
        255.0,
        &PlanetClass::Compact,
        &AtmosphereType::NitrogenOxygen,
    );

    // Earth's effective temperature is about 288 K (255 + 33 greenhouse)
    assert!(
        (temp.effective - 288.0).abs() < 5.0,
        "Earth effective temp should be ~288K, got {}",
        temp.effective
    );
}

#[test]
fn test_effective_temperature_venus() {
    let temp = Temperature::calculate(
        0.815,
        230.0,
        &PlanetClass::Compact,
        &AtmosphereType::ThickCO2,
    );

    // Venus's effective temperature is about 735 K
    assert!(
        (temp.effective - 730.0).abs() < 20.0,
        "Venus effective temp should be ~730K, got {}",
        temp.effective
    );
}

#[test]
fn test_effective_temperature_jupiter() {
    let temp = Temperature::calculate(
        318.0,
        110.0,
        &PlanetClass::Giant,
        &AtmosphereType::HydrogenHelium,
    );

    // Jupiter's effective temperature includes internal heat
    assert!(
        temp.effective > 150.0,
        "Jupiter effective temp should be >150K (with internal heat), got {}",
        temp.effective
    );
}

// ========== Edge Cases and Boundaries ==========

#[test]
fn test_temperature_very_close_orbit() {
    let temp = Temperature::equilibrium_temperature(0.01, 1.0, 0.1);
    // Should be extremely hot
    assert!(temp > 2000.0, "Ultra-close orbit should be >2000K");
}

#[test]
fn test_temperature_very_distant_orbit() {
    let temp = Temperature::equilibrium_temperature(100.0, 1.0, 0.3);
    // Should be very cold
    assert!(temp < 30.0, "Very distant orbit should be <30K");
}

#[test]
fn test_temperature_bright_star() {
    let temp = Temperature::equilibrium_temperature(1.0, 10.0, 0.3);
    // 10x solar luminosity should produce much higher temperature
    assert!(
        temp > 400.0,
        "Bright star should produce high temperature, got {}K",
        temp
    );
}

#[test]
fn test_temperature_dim_star() {
    let temp = Temperature::equilibrium_temperature(1.0, 0.1, 0.3);
    // 0.1x solar luminosity should produce lower temperature
    assert!(
        temp < 150.0,
        "Dim star should produce low temperature, got {}K",
        temp
    );
}

#[test]
fn test_temperature_struct_fields() {
    let temp = Temperature::calculate(
        1.0,
        255.0,
        &PlanetClass::Compact,
        &AtmosphereType::NitrogenOxygen,
    );

    // Verify both fields are set correctly
    assert_eq!(temp.equilibrium, 255.0);
    assert!(temp.effective > temp.equilibrium);
}

#[test]
fn test_dwarf_planet_minimal_heating() {
    let temp = Temperature::calculate(0.01, 100.0, &PlanetClass::Compact, &AtmosphereType::None);

    // Dwarf planets have negligible internal heating
    assert!(
        (temp.effective - temp.equilibrium).abs() < 1.0,
        "Dwarf planets should have minimal internal heating"
    );
}

// ========== Temperature Classification Tests ==========

#[test]
fn test_temperature_class_cold() {
    assert_eq!(TemperatureClass::classify(50.0), TemperatureClass::Cold);
    assert_eq!(TemperatureClass::classify(149.0), TemperatureClass::Cold);
}

#[test]
fn test_temperature_class_temperate() {
    assert_eq!(
        TemperatureClass::classify(150.0),
        TemperatureClass::Temperate
    );
    assert_eq!(
        TemperatureClass::classify(288.0),
        TemperatureClass::Temperate
    );
    assert_eq!(
        TemperatureClass::classify(399.0),
        TemperatureClass::Temperate
    );
}

#[test]
fn test_temperature_class_warm() {
    assert_eq!(TemperatureClass::classify(400.0), TemperatureClass::Warm);
    assert_eq!(TemperatureClass::classify(700.0), TemperatureClass::Warm);
}

#[test]
fn test_temperature_class_hot() {
    assert_eq!(TemperatureClass::classify(1000.0), TemperatureClass::Hot);
    assert_eq!(TemperatureClass::classify(1500.0), TemperatureClass::Hot);
}

#[test]
fn test_temperature_class_ultrahot() {
    assert_eq!(
        TemperatureClass::classify(2000.0),
        TemperatureClass::UltraHot
    );
    assert_eq!(
        TemperatureClass::classify(3000.0),
        TemperatureClass::UltraHot
    );
}

#[test]
fn test_temperature_class_ranges() {
    assert_eq!(TemperatureClass::Cold.range(), (0.0, 150.0));
    assert_eq!(TemperatureClass::Temperate.range(), (150.0, 400.0));
    assert_eq!(TemperatureClass::Warm.range(), (400.0, 1000.0));
    assert_eq!(TemperatureClass::Hot.range(), (1000.0, 2000.0));
    assert_eq!(TemperatureClass::UltraHot.range().0, 2000.0);
}

#[test]
fn test_temperature_class_habitability() {
    assert!(!TemperatureClass::Cold.potentially_habitable());
    assert!(TemperatureClass::Temperate.potentially_habitable());
    assert!(!TemperatureClass::Warm.potentially_habitable());
    assert!(!TemperatureClass::Hot.potentially_habitable());
    assert!(!TemperatureClass::UltraHot.potentially_habitable());
}

#[test]
fn test_temperature_class_display() {
    assert_eq!(format!("{}", TemperatureClass::Cold), "Cold");
    assert_eq!(format!("{}", TemperatureClass::Temperate), "Temperate");
    assert_eq!(format!("{}", TemperatureClass::UltraHot), "Ultra-Hot");
}

// ========== Albedo Estimation Tests ==========

#[test]
fn test_estimate_albedo_gas_giant() {
    let cold = Temperature::estimate_albedo(&PlanetClass::Giant, 100.0);
    let hot = Temperature::estimate_albedo(&PlanetClass::Giant, 1500.0);

    assert!(cold > 0.4, "Cold gas giant should have high albedo");
    assert!(hot < 0.2, "Hot gas giant should have low albedo");
}

#[test]
fn test_estimate_albedo_rocky() {
    let icy = Temperature::estimate_albedo(&PlanetClass::Compact, 100.0);
    let lava = Temperature::estimate_albedo(&PlanetClass::Compact, 1500.0);

    assert!(icy > 0.5, "Icy world should have high albedo");
    assert!(lava < 0.2, "Lava world should have low albedo");
}
