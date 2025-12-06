use crate::{planet_class::PlanetClass, planet_type::PlanetType};

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
