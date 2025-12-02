//! Tests for size distributions.

use approx::assert_relative_eq;
use units::{Density, Length, Mass};

use crate::particles::SizeDistribution;

// Helper to create common test values
fn silicate_density() -> Density {
    Density::from_grams_per_cm3(3.0)
}

fn test_mass() -> Mass {
    Mass::from_grams(1000.0)
}

// =============================================================================
// Constructor tests
// =============================================================================

#[test]
fn power_law_stores_parameters() {
    let dist = SizeDistribution::power_law(
        Length::from_cm(1e-4),
        Length::from_cm(0.1),
        3.5,
        test_mass(),
        silicate_density(),
    );

    assert_relative_eq!(dist.min_size().to_cm(), 1e-4);
    assert_relative_eq!(dist.max_size().to_cm(), 0.1);
    assert_relative_eq!(dist.total_mass().to_grams(), 1000.0);
    assert_relative_eq!(dist.material_density().to_grams_per_cm3(), 3.0);
}

#[test]
fn mrn_has_correct_parameters() {
    let dist = SizeDistribution::mrn(test_mass(), silicate_density());

    // MRN: 0.1 μm to 1 μm, q = 3.5
    assert_relative_eq!(dist.min_size().to_cm(), 1e-5);
    assert_relative_eq!(dist.max_size().to_cm(), 1e-4);
    assert_relative_eq!(dist.total_mass().to_grams(), 1000.0);
}

#[test]
fn monodisperse_stores_parameters() {
    let dist = SizeDistribution::monodisperse(
        Length::from_cm(0.01),
        Mass::from_grams(500.0),
        Density::from_grams_per_cm3(2.5),
    );

    assert_relative_eq!(dist.min_size().to_cm(), 0.01);
    assert_relative_eq!(dist.max_size().to_cm(), 0.01);
    assert_relative_eq!(dist.total_mass().to_grams(), 500.0);
    assert_relative_eq!(dist.material_density().to_grams_per_cm3(), 2.5);
}

#[test]
fn binned_empty_has_zero_mass() {
    let dist = SizeDistribution::binned_empty(
        Length::from_cm(1e-4),
        Length::from_cm(1.0),
        20,
        silicate_density(),
    );

    assert_relative_eq!(dist.total_mass().to_grams(), 0.0);
    assert_relative_eq!(dist.min_size().to_cm(), 1e-4);
    assert_relative_eq!(dist.max_size().to_cm(), 1.0);
}

// =============================================================================
// Mass conservation tests
// =============================================================================

#[test]
fn power_law_mass_preserved() {
    let dist = SizeDistribution::power_law(
        Length::from_cm(1e-4),
        Length::from_cm(0.1),
        3.5,
        Mass::from_grams(1234.5),
        silicate_density(),
    );
    assert_relative_eq!(dist.total_mass().to_grams(), 1234.5);
}

#[test]
fn to_binned_preserves_mass_power_law() {
    let power_law = SizeDistribution::power_law(
        Length::from_cm(1e-4),
        Length::from_cm(0.1),
        3.5,
        test_mass(),
        silicate_density(),
    );
    let binned = power_law.to_binned(50);

    assert_relative_eq!(
        power_law.total_mass().to_grams(),
        binned.total_mass().to_grams(),
        epsilon = 1e-10
    );
}

#[test]
fn to_binned_preserves_mass_mrn() {
    let mrn = SizeDistribution::mrn(Mass::from_grams(5000.0), silicate_density());
    let binned = mrn.to_binned(100);

    assert_relative_eq!(
        mrn.total_mass().to_grams(),
        binned.total_mass().to_grams(),
        epsilon = 1e-10
    );
}

#[test]
fn to_binned_preserves_mass_monodisperse() {
    let mono = SizeDistribution::monodisperse(
        Length::from_cm(0.01),
        Mass::from_grams(750.0),
        silicate_density(),
    );
    let binned = mono.to_binned(30);

    assert_relative_eq!(
        mono.total_mass().to_grams(),
        binned.total_mass().to_grams(),
        epsilon = 1e-10
    );
}

#[test]
fn to_binned_different_exponents() {
    // Test several exponents including edge cases
    for q in [2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0] {
        let power_law = SizeDistribution::power_law(
            Length::from_cm(1e-4),
            Length::from_cm(0.1),
            q,
            test_mass(),
            silicate_density(),
        );
        let binned = power_law.to_binned(50);

        let m1 = power_law.total_mass().to_grams();
        let m2 = binned.total_mass().to_grams();
        assert!(
            (m1 - m2).abs() < 1e-10,
            "mass not preserved for q = {}: {} vs {}",
            q,
            m1,
            m2
        );
    }
}

// =============================================================================
// Mean size tests
// =============================================================================

#[test]
fn monodisperse_mean_equals_size() {
    let dist =
        SizeDistribution::monodisperse(Length::from_cm(0.05), test_mass(), silicate_density());
    assert_relative_eq!(dist.mean_size().to_cm(), 0.05);
}

#[test]
fn mean_size_between_bounds() {
    let dist = SizeDistribution::power_law(
        Length::from_cm(1e-4),
        Length::from_cm(0.1),
        3.5,
        test_mass(),
        silicate_density(),
    );
    let mean = dist.mean_size().to_cm();

    assert!(mean > 1e-4, "mean {} should be > s_min {}", mean, 1e-4);
    assert!(mean < 0.1, "mean {} should be < s_max {}", mean, 0.1);
}

#[test]
fn mean_size_shifts_with_exponent() {
    // Lower exponent -> more mass in large particles -> larger mean
    let steep = SizeDistribution::power_law(
        Length::from_cm(1e-4),
        Length::from_cm(0.1),
        4.5,
        test_mass(),
        silicate_density(),
    );
    let shallow = SizeDistribution::power_law(
        Length::from_cm(1e-4),
        Length::from_cm(0.1),
        2.5,
        test_mass(),
        silicate_density(),
    );

    assert!(
        shallow.mean_size().to_cm() > steep.mean_size().to_cm(),
        "shallow ({}) should have larger mean than steep ({})",
        shallow.mean_size().to_cm(),
        steep.mean_size().to_cm()
    );
}

#[test]
fn binned_mean_close_to_power_law() {
    let power_law = SizeDistribution::power_law(
        Length::from_cm(1e-4),
        Length::from_cm(0.1),
        3.5,
        test_mass(),
        silicate_density(),
    );
    let binned = power_law.to_binned(100);

    // Should be within a few percent
    assert_relative_eq!(
        power_law.mean_size().to_cm(),
        binned.mean_size().to_cm(),
        epsilon = 0.05
    );
}

// =============================================================================
// Total number tests
// =============================================================================

#[test]
fn monodisperse_number_correct() {
    let size = 0.01; // cm
    let rho = 3.0; // g/cm³
    let total_mass = 1000.0; // g

    let dist = SizeDistribution::monodisperse(
        Length::from_cm(size),
        Mass::from_grams(total_mass),
        Density::from_grams_per_cm3(rho),
    );

    // m_particle = (4π/3) ρ s³
    let m_particle = (4.0 / 3.0) * std::f64::consts::PI * rho * size.powi(3);
    let expected_n = total_mass / m_particle;

    assert_relative_eq!(dist.total_number(), expected_n, epsilon = 1e-10);
}

#[test]
fn power_law_number_positive() {
    let dist = SizeDistribution::power_law(
        Length::from_cm(1e-4),
        Length::from_cm(0.1),
        3.5,
        test_mass(),
        silicate_density(),
    );
    assert!(dist.total_number() > 0.0);
}

#[test]
fn mrn_dominated_by_small_particles() {
    // MRN has q = 3.5 > 1, so number is dominated by small particles
    let dist = SizeDistribution::mrn(test_mass(), silicate_density());

    // Number should be huge (small particles are numerous)
    let n = dist.total_number();
    assert!(n > 1e10, "MRN should have many particles, got {}", n);
}

// =============================================================================
// Edge cases
// =============================================================================

#[test]
fn zero_mass_distribution() {
    let dist = SizeDistribution::power_law(
        Length::from_cm(1e-4),
        Length::from_cm(0.1),
        3.5,
        Mass::from_grams(0.0),
        silicate_density(),
    );

    assert_relative_eq!(dist.total_mass().to_grams(), 0.0);
    assert_relative_eq!(dist.total_number(), 0.0);
}

#[test]
fn narrow_size_range() {
    // Very narrow range - should still work
    let dist = SizeDistribution::power_law(
        Length::from_cm(0.01),
        Length::from_cm(0.0100001),
        3.5,
        Mass::from_grams(100.0),
        silicate_density(),
    );

    assert_relative_eq!(dist.total_mass().to_grams(), 100.0);
    // Mean should be very close to both bounds
    let mean = dist.mean_size().to_cm();
    assert!(mean > 0.01 && mean < 0.0100001);
}

#[test]
fn exponent_exactly_four() {
    // q = 4 is a special case (logarithmic integrals)
    let dist = SizeDistribution::power_law(
        Length::from_cm(1e-4),
        Length::from_cm(0.1),
        4.0,
        test_mass(),
        silicate_density(),
    );

    assert_relative_eq!(dist.total_mass().to_grams(), 1000.0);

    let mean = dist.mean_size().to_cm();
    assert!(mean > 1e-4 && mean < 0.1);

    let binned = dist.to_binned(50);
    assert_relative_eq!(
        dist.total_mass().to_grams(),
        binned.total_mass().to_grams(),
        epsilon = 1e-10
    );
}

#[test]
fn to_binned_is_idempotent() {
    let power_law = SizeDistribution::power_law(
        Length::from_cm(1e-4),
        Length::from_cm(0.1),
        3.5,
        test_mass(),
        silicate_density(),
    );
    let binned1 = power_law.to_binned(50);
    let binned2 = binned1.to_binned(50);

    // Converting an already-binned distribution should be a no-op
    assert_relative_eq!(
        binned1.total_mass().to_grams(),
        binned2.total_mass().to_grams()
    );
}
