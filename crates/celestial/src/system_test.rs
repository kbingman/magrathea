use crate::{
    GenerationMethod, HabitableZone, PlanetarySystem, SystemArchitecture, SystemMetadata, snow_line,
};
use planetary::planet::earth_analog;
use stellar::{StellarObject, solar_analog};
use units::Length;

#[test]
fn test_habitable_zone() {
    // Sun-like star
    let hz = HabitableZone::from_luminosity(1.0);
    assert!(
        hz.inner_edge < 1.0 && hz.inner_edge > 0.8,
        "Inner edge {} should be between 0.8 and 1.0 AU",
        hz.inner_edge
    );
    assert!(
        hz.outer_edge > 1.5 && hz.outer_edge < 2.0,
        "Outer edge {} should be between 1.5 and 2.0 AU",
        hz.outer_edge
    );

    // M dwarf (L = 0.01 L☉)
    let hz_m = HabitableZone::from_luminosity(0.01);
    assert!(
        hz_m.inner_edge < 0.15,
        "M dwarf inner edge {} should be < 0.15 AU",
        hz_m.inner_edge
    );
    assert!(
        hz_m.outer_edge < 0.3,
        "M dwarf outer edge {} should be < 0.3 AU",
        hz_m.outer_edge
    );
}

#[test]
fn test_snow_line() {
    // Sun-like: snow line ~2.7 AU
    let sl = snow_line(1.0);
    assert!(sl > 2.0 && sl < 4.0, "Snow line {} should be 2-4 AU", sl);

    // M dwarf: much closer
    let sl_m = snow_line(0.01);
    assert!(sl_m < 0.5, "M dwarf snow line {} should be < 0.5 AU", sl_m);
}

#[test]
fn test_system_stability() {
    let earth = earth_analog();
    let mut earth2 = earth.clone();
    earth2.semi_major_axis = Length::from_au(1.01); // Way too close

    let star = StellarObject::MainSequence(solar_analog());
    let metadata = SystemMetadata::new_random(GenerationMethod::Manual, SystemArchitecture::Mixed);

    let system = PlanetarySystem::new(vec![star], vec![earth, earth2], metadata);

    assert!(
        !system.is_stable(),
        "Planets at 1.0 and 1.01 AU should be unstable"
    );
}

#[test]
fn test_system_accessors() {
    let star = StellarObject::MainSequence(solar_analog());
    let metadata = SystemMetadata::new_random(
        GenerationMethod::Statistical,
        SystemArchitecture::CompactMulti,
    );

    let system = PlanetarySystem::new(vec![star], vec![earth_analog()], metadata);

    // Test accessor methods
    assert!(!system.is_binary());
    assert!(system.effective_mass() > 0.9 && system.effective_mass() < 1.1);
    assert!(system.total_luminosity() > 0.5 && system.total_luminosity() < 2.0);
    assert_eq!(system.architecture(), SystemArchitecture::CompactMulti);

    // HZ should contain Earth
    let hz = system.habitable_zone();
    assert!(hz.inner_edge < 1.0 && hz.outer_edge > 1.0);

    // Snow line should be reasonable
    let sl = system.snow_line();
    assert!(sl > 2.0 && sl < 4.0);
}

#[test]
fn test_planets_sorted_by_sma() {
    let star = StellarObject::MainSequence(solar_analog());
    let metadata = SystemMetadata::new_random(GenerationMethod::Manual, SystemArchitecture::Mixed);

    let mut inner = earth_analog();
    inner.semi_major_axis = Length::from_au(0.5);

    let mut outer = earth_analog();
    outer.semi_major_axis = Length::from_au(2.0);

    // Pass planets in wrong order
    let system = PlanetarySystem::new(vec![star], vec![outer, inner], metadata);

    // Should be sorted by SMA
    assert!(
        system.planets[0].semi_major_axis.to_au() < system.planets[1].semi_major_axis.to_au(),
        "Planets should be sorted by semi-major axis"
    );
}

#[test]
fn test_planet_naming() {
    let star = StellarObject::MainSequence(solar_analog());
    let metadata = SystemMetadata::new_random(GenerationMethod::Manual, SystemArchitecture::Mixed);
    let catalog_name = metadata.catalog_name.clone();
    let system_id = metadata.id.to_string();

    let mut inner = earth_analog();
    inner.semi_major_axis = Length::from_au(0.5);

    let mut middle = earth_analog();
    middle.semi_major_axis = Length::from_au(1.0);

    let mut outer = earth_analog();
    outer.semi_major_axis = Length::from_au(2.0);

    // Pass planets in scrambled order
    let system = PlanetarySystem::new(vec![star], vec![outer, inner, middle], metadata);

    // Check IDs are assigned correctly (innermost = b)
    assert_eq!(system.planets[0].id, format!("{}-b", system_id));
    assert_eq!(system.planets[1].id, format!("{}-c", system_id));
    assert_eq!(system.planets[2].id, format!("{}-d", system_id));

    // Check names follow IAU convention
    assert_eq!(system.planets[0].name, format!("{} b", catalog_name));
    assert_eq!(system.planets[1].name, format!("{} c", catalog_name));
    assert_eq!(system.planets[2].name, format!("{} d", catalog_name));
}
