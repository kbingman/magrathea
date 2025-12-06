use crate::planet_class::PlanetClass;

#[test]
fn test_classification_boundaries() {
    assert_eq!(PlanetClass::from_earth_masses(1.0), PlanetClass::Compact);
    assert_eq!(
        PlanetClass::from_earth_masses(2.0),
        PlanetClass::Transitional
    );
    assert_eq!(PlanetClass::from_earth_masses(5.0), PlanetClass::Volatile);
    assert_eq!(PlanetClass::from_earth_masses(160.0), PlanetClass::Giant);
}

#[test]
fn test_solar_system() {
    assert_eq!(PlanetClass::from_earth_masses(1.0), PlanetClass::Compact); // Earth
    assert_eq!(PlanetClass::from_earth_masses(17.1), PlanetClass::Volatile); // Neptune
    assert_eq!(PlanetClass::from_earth_masses(317.8), PlanetClass::Giant); // Jupiter
}
