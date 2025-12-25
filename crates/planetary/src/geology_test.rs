//! Tests for geological activity classification

use crate::geology::{GeologicalActivity, TectonicRegime, VolcanismLevel};
use crate::interior::ThermalState;
use crate::planet_class::PlanetClass;

#[test]
fn test_earth_like_geology() {
    let thermal = ThermalState::calculate(1.0, 4.5, 0.01);
    let geology = GeologicalActivity::calculate(
        &thermal,
        1.0,
        4.5,
        288.0,
        &PlanetClass::Compact,
        false, // not water-rich
    );

    // Earth should have active volcanism
    assert_eq!(geology.volcanism, VolcanismLevel::Active);

    // Earth should have plate tectonics
    assert_eq!(geology.tectonics, TectonicRegime::ActivePlates);

    // Surface should be relatively young
    assert!(geology.surface_age < 1.0);

    // Crustal thickness ~35 km
    assert!(
        geology.crustal_thickness > 30.0 && geology.crustal_thickness < 60.0,
        "Crustal thickness should be ~35-50 km, got {}",
        geology.crustal_thickness
    );
}

#[test]
fn test_mars_like_geology() {
    // Mars: 0.107 M⊕
    let thermal = ThermalState::calculate(0.107, 4.5, 0.01);
    let geology =
        GeologicalActivity::calculate(&thermal, 0.107, 4.5, 210.0, &PlanetClass::Compact, false);

    // Mars should have ancient or no volcanism
    assert!(
        matches!(
            geology.volcanism,
            VolcanismLevel::Ancient | VolcanismLevel::None
        ),
        "Mars volcanism should be Ancient or None, got {:?}",
        geology.volcanism
    );

    // Mars should have stagnant lid
    assert_eq!(geology.tectonics, TectonicRegime::StagnantLid);

    // Surface should be old
    assert!(geology.surface_age > 3.0);
}

#[test]
fn test_io_like_hyperactive() {
    // Io-like: 0.015 M⊕, tidally heated (factor ~100-200x)
    let thermal = ThermalState::calculate(0.015, 4.5, 150.0);
    let geology =
        GeologicalActivity::calculate(&thermal, 0.015, 4.5, 110.0, &PlanetClass::Compact, false);

    // Io should be hyperactive
    assert_eq!(geology.volcanism, VolcanismLevel::Hyperactive);

    // Surface should be very young
    assert!(geology.surface_age < 0.1);
}

#[test]
fn test_ocean_world_cryovolcanism() {
    // Cold water-rich world
    let thermal = ThermalState::calculate(3.0, 4.5, 0.01);
    let geology = GeologicalActivity::calculate(
        &thermal,
        3.0,
        4.5,
        200.0,
        &PlanetClass::Transitional,
        true, // water-rich
    );

    // Should have cryovolcanism
    assert_eq!(geology.volcanism, VolcanismLevel::Cryovolcanic);
}

#[test]
fn test_gas_giant_no_volcanism() {
    let thermal = ThermalState::calculate(318.0, 4.5, 0.01);
    let geology =
        GeologicalActivity::calculate(&thermal, 318.0, 4.5, 165.0, &PlanetClass::Giant, false);

    // Gas giants have no volcanism (no solid surface)
    assert_eq!(geology.volcanism, VolcanismLevel::None);

    // No tectonics
    assert_eq!(geology.tectonics, TectonicRegime::None);

    // No crust
    assert_eq!(geology.crustal_thickness, 0.0);
}

#[test]
fn test_ice_giant_no_volcanism() {
    let thermal = ThermalState::calculate(17.0, 4.5, 0.01);
    let geology =
        GeologicalActivity::calculate(&thermal, 17.0, 4.5, 72.0, &PlanetClass::Volatile, false);

    // Ice giants have no volcanism
    assert_eq!(geology.volcanism, VolcanismLevel::None);

    // No tectonics
    assert_eq!(geology.tectonics, TectonicRegime::None);
}

#[test]
fn test_young_planet_active() {
    // Young Earth-sized planet
    let thermal = ThermalState::calculate(1.0, 0.5, 0.01);
    let geology =
        GeologicalActivity::calculate(&thermal, 1.0, 0.5, 300.0, &PlanetClass::Compact, false);

    // Young planets should be active
    assert!(
        matches!(
            geology.volcanism,
            VolcanismLevel::Active | VolcanismLevel::Hyperactive
        ),
        "Young planet should be volcanically active"
    );
}

#[test]
fn test_small_body_dead() {
    // Moon-sized body
    let thermal = ThermalState::calculate(0.012, 4.5, 0.01);
    let geology =
        GeologicalActivity::calculate(&thermal, 0.012, 4.5, 250.0, &PlanetClass::Compact, false);

    // Small old bodies are dead
    assert_eq!(geology.volcanism, VolcanismLevel::None);

    // Stagnant lid
    assert_eq!(geology.tectonics, TectonicRegime::StagnantLid);
}

#[test]
fn test_venus_like_episodic() {
    // Venus: 0.815 M⊕, hot, no water
    let thermal = ThermalState::calculate(0.815, 4.5, 0.01);
    let geology = GeologicalActivity::calculate(
        &thermal,
        0.815,
        4.5,
        737.0, // Venus surface temp
        &PlanetClass::Compact,
        false,
    );

    // Should have stagnant lid or episodic overturn (no water for plate tectonics)
    assert!(
        matches!(
            geology.tectonics,
            TectonicRegime::StagnantLid | TectonicRegime::EpisodicOverturn
        ),
        "Venus should have stagnant lid or episodic overturn, got {:?}",
        geology.tectonics
    );
}

#[test]
fn test_surface_age_bounds() {
    let thermal = ThermalState::calculate(1.0, 4.5, 0.01);
    let geology =
        GeologicalActivity::calculate(&thermal, 1.0, 4.5, 288.0, &PlanetClass::Compact, false);

    // Surface age should not exceed planetary age
    assert!(geology.surface_age <= 4.5);
    assert!(geology.surface_age >= 0.0);
}

#[test]
fn test_crustal_thickness_scales_with_mass() {
    let thermal1 = ThermalState::calculate(0.1, 4.5, 0.01);
    let thermal2 = ThermalState::calculate(2.0, 4.5, 0.01);

    let small =
        GeologicalActivity::calculate(&thermal1, 0.1, 4.5, 250.0, &PlanetClass::Compact, false);
    let large =
        GeologicalActivity::calculate(&thermal2, 2.0, 4.5, 300.0, &PlanetClass::Compact, false);

    // Larger planets should have thicker crust
    assert!(
        large.crustal_thickness > small.crustal_thickness,
        "Large planet should have thicker crust"
    );
}

#[test]
fn test_cold_distant_world_dead() {
    // Cold outer system body
    let thermal = ThermalState::calculate(0.5, 4.5, 0.01);
    let geology = GeologicalActivity::calculate(
        &thermal,
        0.5,
        4.5,
        50.0, // Very cold
        &PlanetClass::Compact,
        false,
    );

    // Should be geologically dead
    assert_eq!(geology.volcanism, VolcanismLevel::None);
}

#[test]
fn test_is_active() {
    let thermal = ThermalState::calculate(1.0, 4.5, 0.01);

    let active =
        GeologicalActivity::calculate(&thermal, 1.0, 4.5, 288.0, &PlanetClass::Compact, false);
    assert!(active.is_active());

    let dead = GeologicalActivity::calculate(
        &ThermalState::calculate(0.012, 4.5, 0.01),
        0.012,
        4.5,
        250.0,
        &PlanetClass::Compact,
        false,
    );
    assert!(!dead.is_active());
}

#[test]
fn test_has_plate_tectonics() {
    let thermal = ThermalState::calculate(1.0, 4.5, 0.01);
    let earth =
        GeologicalActivity::calculate(&thermal, 1.0, 4.5, 288.0, &PlanetClass::Compact, false);
    assert!(earth.has_plate_tectonics());

    let venus =
        GeologicalActivity::calculate(&thermal, 0.815, 4.5, 737.0, &PlanetClass::Compact, false);
    assert!(!venus.has_plate_tectonics());
}

#[test]
fn test_volcanism_display() {
    assert_eq!(format!("{}", VolcanismLevel::Active), "Active");
    assert_eq!(format!("{}", VolcanismLevel::Cryovolcanic), "Cryovolcanic");
}

#[test]
fn test_tectonic_display() {
    assert_eq!(format!("{}", TectonicRegime::ActivePlates), "Active Plates");
    assert_eq!(format!("{}", TectonicRegime::StagnantLid), "Stagnant Lid");
}
