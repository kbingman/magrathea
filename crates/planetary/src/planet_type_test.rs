use crate::{planet_class::PlanetClass, planet_type::PlanetType};

#[test]
fn test_class_mapping() {
    assert_eq!(
        PlanetType::Terran {
            ocean_fraction: 0.7,
            tectonically_active: true
        }
        .class(),
        PlanetClass::Rocky
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
