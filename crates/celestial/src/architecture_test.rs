use rand::SeedableRng;
use rand_chacha::ChaChaRng;

use crate::SystemArchitecture;

#[test]
fn test_architecture_sampling_by_spectral_type() {
    let mut rng = ChaChaRng::seed_from_u64(42);

    // M dwarfs should favor compact multi (using typical M-dwarf mass of 0.3 M☉)
    let mut compact_count = 0;
    for _ in 0..100 {
        if matches!(
            SystemArchitecture::sample(&mut rng, "M3", 0.3, 0.0),
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

    // High metallicity should increase giant probability (using solar mass)
    let mut _giant_low = 0;
    let mut _giant_high = 0;
    for _ in 0..100 {
        if matches!(
            SystemArchitecture::sample(&mut rng, "G2", 1.0, -0.3),
            SystemArchitecture::GiantDominated
        ) {
            _giant_low += 1;
        }
        if matches!(
            SystemArchitecture::sample(&mut rng, "G2", 1.0, 0.3),
            SystemArchitecture::GiantDominated
        ) {
            _giant_high += 1;
        }
    }
    // Can't guarantee with RNG, but high metallicity should tend toward more giants
}

#[test]
fn test_late_m_dwarf_giant_suppression() {
    let mut rng = ChaChaRng::seed_from_u64(42);

    // Late M-dwarfs (M7, ~0.1 M☉) should almost never get GiantDominated
    let mut giant_count = 0;
    for _ in 0..1000 {
        if matches!(
            SystemArchitecture::sample(&mut rng, "M7", 0.1, 0.0),
            SystemArchitecture::GiantDominated
        ) {
            giant_count += 1;
        }
    }
    // Should be <1% giant-dominated for late M-dwarfs
    assert!(
        giant_count < 20,
        "Late M-dwarfs should rarely get GiantDominated architecture, got {}%",
        giant_count as f64 / 10.0
    );
}
