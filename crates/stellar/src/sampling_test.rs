use rand::SeedableRng;
use rand_chacha::ChaChaRng;

use crate::sampling::{sample_gaussian, sample_mass_kroupa, sample_metallicity, sample_power_law};

#[test]
fn sample_gaussian_produces_reasonable_values() {
    let mut rng = ChaChaRng::seed_from_u64(42);

    // Sample many values and check they're roughly centered on mean
    let samples: Vec<f64> = (0..1000)
        .map(|_| sample_gaussian(&mut rng, 5.0, 1.0))
        .collect();
    let mean: f64 = samples.iter().sum::<f64>() / samples.len() as f64;

    // Mean should be close to 5.0
    assert!(
        (mean - 5.0).abs() < 0.2,
        "Mean {} should be close to 5.0",
        mean
    );

    // Check standard deviation is roughly 1.0
    let variance: f64 =
        samples.iter().map(|x| (x - mean).powi(2)).sum::<f64>() / samples.len() as f64;
    let std_dev = variance.sqrt();
    assert!(
        (std_dev - 1.0).abs() < 0.2,
        "Std dev {} should be close to 1.0",
        std_dev
    );
}

#[test]
fn sample_power_law_respects_bounds() {
    let mut rng = ChaChaRng::seed_from_u64(42);

    for _ in 0..100 {
        let sample = sample_power_law(0.5, 10.0, -2.3, &mut rng);
        assert!(sample >= 0.5, "Sample {} should be >= 0.5", sample);
        assert!(sample <= 10.0, "Sample {} should be <= 10.0", sample);
    }
}

#[test]
fn sample_power_law_negative_exponent_favors_lower_values() {
    let mut rng = ChaChaRng::seed_from_u64(42);

    // With negative exponent, lower values should be more common
    let samples: Vec<f64> = (0..1000)
        .map(|_| sample_power_law(1.0, 100.0, -2.3, &mut rng))
        .collect();
    let median = {
        let mut sorted = samples.clone();
        sorted.sort_by(|a, b| a.partial_cmp(b).unwrap());
        sorted[500]
    };

    // Median should be much closer to min than max for steep negative power law
    assert!(
        median < 20.0,
        "Median {} should be < 20 for steep power law",
        median
    );
}

#[test]
fn sample_metallicity_within_bounds() {
    let mut rng = ChaChaRng::seed_from_u64(42);

    for _ in 0..100 {
        let metallicity = sample_metallicity(&mut rng);
        assert!(
            metallicity >= -0.5,
            "Metallicity {} should be >= -0.5",
            metallicity
        );
        assert!(
            metallicity <= 0.4,
            "Metallicity {} should be <= 0.4",
            metallicity
        );
    }
}

#[test]
fn sample_metallicity_centered_on_solar() {
    let mut rng = ChaChaRng::seed_from_u64(42);

    let samples: Vec<f64> = (0..1000).map(|_| sample_metallicity(&mut rng)).collect();
    let mean: f64 = samples.iter().sum::<f64>() / samples.len() as f64;

    // Mean should be close to 0.0 (solar)
    assert!(
        mean.abs() < 0.1,
        "Mean metallicity {} should be close to 0.0",
        mean
    );
}

#[test]
fn sample_mass_kroupa_low_mass_segment() {
    // Use a seed that hits the low mass segment (rand < 0.80)
    let mut rng = ChaChaRng::seed_from_u64(100);

    // Sample many times - most should be in low mass range
    let mut low_mass_count = 0;
    for _ in 0..100 {
        let mass = sample_mass_kroupa(&mut rng, 150.0);
        if mass < 0.5 {
            low_mass_count += 1;
        }
    }

    // Should get a significant number of low mass stars
    assert!(
        low_mass_count > 50,
        "Expected >50 low mass stars, got {}",
        low_mass_count
    );
}

#[test]
fn sample_mass_kroupa_intermediate_segment() {
    let mut rng = ChaChaRng::seed_from_u64(42);

    // Sample many times to hit all segments
    let mut intermediate_count = 0;
    for _ in 0..1000 {
        let mass = sample_mass_kroupa(&mut rng, 150.0);
        if mass >= 0.5 && mass < 1.0 {
            intermediate_count += 1;
        }
    }

    // Should get some intermediate mass stars (~15% expected)
    assert!(
        intermediate_count > 50,
        "Expected >50 intermediate mass stars, got {}",
        intermediate_count
    );
}

#[test]
fn sample_mass_kroupa_high_mass_segment() {
    let mut rng = ChaChaRng::seed_from_u64(42);

    // Sample many times to hit high mass segment
    let mut high_mass_count = 0;
    for _ in 0..1000 {
        let mass = sample_mass_kroupa(&mut rng, 150.0);
        if mass >= 1.0 {
            high_mass_count += 1;
        }
    }

    // Should get some high mass stars (~5% expected)
    assert!(
        high_mass_count > 10,
        "Expected >10 high mass stars, got {}",
        high_mass_count
    );
}

#[test]
fn sample_mass_kroupa_respects_max_mass() {
    let mut rng = ChaChaRng::seed_from_u64(42);

    // Test that max_mass limits the high mass segment
    // When max_mass = 2.0, stars in the high mass segment should be <= 2.0
    for _ in 0..100 {
        let mass = sample_mass_kroupa(&mut rng, 2.0);
        assert!(
            mass <= 2.0,
            "Mass {} should be <= 2.0 when max_mass=2.0",
            mass
        );
    }
}
