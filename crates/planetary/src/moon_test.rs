use crate::moon::{MoonFormation, MoonSystem, moon_numeral};

#[test]
fn test_moon_numeral() {
    assert_eq!(moon_numeral(0), "I");
    assert_eq!(moon_numeral(1), "II");
    assert_eq!(moon_numeral(2), "III");
    assert_eq!(moon_numeral(3), "IV");
    assert_eq!(moon_numeral(4), "V");
    assert_eq!(moon_numeral(5), "VI");
    assert_eq!(moon_numeral(6), "VII");
    assert_eq!(moon_numeral(7), "VIII");
    assert_eq!(moon_numeral(8), "IX");
    assert_eq!(moon_numeral(9), "X");
    assert_eq!(moon_numeral(13), "XIV");
    assert_eq!(moon_numeral(19), "XX");
}

#[test]
fn test_moon_system_basics() {
    let system = MoonSystem::new();
    assert!(!system.has_moons());
    assert_eq!(system.moon_count(), 0);
    assert!(!system.has_rings);
}

#[test]
fn test_moon_formation_display() {
    assert_eq!(MoonFormation::GiantImpact.to_string(), "Giant Impact");
    assert_eq!(
        MoonFormation::Capture { retrograde: true }.to_string(),
        "Captured (retrograde)"
    );
    assert_eq!(
        MoonFormation::CoAccretion {
            resonance_order: Some((2, 1))
        }
        .to_string(),
        "Co-accretion (2:1 resonance)"
    );
}
