use crate::sampling::sample_planet_mass;
use planetary::planet_class::PlanetClass;

use rand::SeedableRng;
use rand_chacha::ChaChaRng;

#[test]
fn test_mass_sampling() {
    let mut rng = ChaChaRng::seed_from_u64(42);
    let mut class_counts = [0usize; 4];

    for _ in 0..1000 {
        let mass = sample_planet_mass(&mut rng, 0.0);
        let class = PlanetClass::from_earth_masses(mass);
        match class {
            PlanetClass::Compact => class_counts[0] += 1,
            PlanetClass::Transitional => class_counts[1] += 1,
            PlanetClass::Volatile => class_counts[2] += 1,
            PlanetClass::Giant => class_counts[3] += 1,
        }
    }

    // Rocky, Transitional, and Volatile should be well-represented
    // (sample_planet_mass is for inner system planets, so no Jupiter-mass giants)
    assert!(class_counts[0] > 100, "Rocky planets should be common");
    assert!(
        class_counts[1] > 100,
        "Transitional planets should be common"
    );
    assert!(class_counts[2] > 50, "Volatile planets should appear");

    // Giant class should be rare or absent (max mass is 160 M⊕, just at threshold)
    // The sub-Saturn bin (50-160 M⊕) is Volatile class, not Giant
    assert!(
        class_counts[3] < 50,
        "Giant planets should be rare from inner sampling"
    );
}
