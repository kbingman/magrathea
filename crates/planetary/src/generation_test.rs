use rand::SeedableRng;
use rand_chacha::ChaChaRng;

use crate::generation::generate_planetary_system;

#[test]
fn test_generate_system() {
    let mut rng = ChaChaRng::seed_from_u64(42);
    let system = generate_planetary_system(&mut rng, 1.0, 1.0, 5778.0, 0.0, "G");

    assert!(system.is_stable());
}

#[test]
fn test_reproducibility() {
    let s1 = {
        let mut rng = ChaChaRng::seed_from_u64(12345);
        generate_planetary_system(&mut rng, 1.0, 1.0, 5778.0, 0.0, "G")
    };
    let s2 = {
        let mut rng = ChaChaRng::seed_from_u64(12345);
        generate_planetary_system(&mut rng, 1.0, 1.0, 5778.0, 0.0, "G")
    };

    assert_eq!(s1.planets.len(), s2.planets.len());
}
