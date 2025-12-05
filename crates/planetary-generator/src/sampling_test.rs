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
            PlanetClass::Rocky => class_counts[0] += 1,
            PlanetClass::Transitional => class_counts[1] += 1,
            PlanetClass::Volatile => class_counts[2] += 1,
            PlanetClass::Giant => class_counts[3] += 1,
        }
    }

    // All classes should be represented
    for count in &class_counts {
        assert!(*count > 0);
    }
}
