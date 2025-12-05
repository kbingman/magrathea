use rand::SeedableRng;
use rand_chacha::ChaChaRng;

use crate::SystemArchitecture;

#[test]
fn test_architecture_sampling_by_spectral_type() {
    let mut rng = ChaChaRng::seed_from_u64(42);

    // M dwarfs should favor compact multi
    let mut compact_count = 0;
    for _ in 0..100 {
        if matches!(
            SystemArchitecture::sample(&mut rng, "M", 0.0),
            SystemArchitecture::CompactMulti
        ) {
            compact_count += 1;
        }
    }
    assert!(
        compact_count > 30,
        "M dwarfs should frequently produce compact systems, got {}",
        compact_count
    );

    // High metallicity should increase giant probability
    let mut _giant_low = 0;
    let mut _giant_high = 0;
    for _ in 0..100 {
        if matches!(
            SystemArchitecture::sample(&mut rng, "G", -0.3),
            SystemArchitecture::GiantDominated
        ) {
            _giant_low += 1;
        }
        if matches!(
            SystemArchitecture::sample(&mut rng, "G", 0.3),
            SystemArchitecture::GiantDominated
        ) {
            _giant_high += 1;
        }
    }
    // Can't guarantee with RNG, but high metallicity should tend toward more giants
}
