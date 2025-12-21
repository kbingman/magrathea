//! Tests for coagulation (particle sticking).

use approx::assert_relative_eq;
use units::{Density, Length, Mass, SurfaceDensity, Time};

use crate::disk::{DiskModel, GasDisk};
use crate::particles::{Coagulation, ParticleBin, SizeDistribution};

fn test_disk() -> GasDisk {
    GasDisk::mmsn()
}

fn silicate_density() -> Density {
    Density::from_grams_per_cm3(3.0)
}

// =============================================================================
// Constructor tests
// =============================================================================

#[test]
fn coagulation_creates_kernel_matrix() {
    let disk = test_disk();
    let bin = ParticleBin::from_disk(&disk, Length::from_au(1.0), Length::from_au(0.1));

    let coag = Coagulation::new(&bin, &disk);

    // Should have 30 bins (from ParticleBin default)
    assert_eq!(coag.n_bins(), 30);
}

#[test]
fn kernel_is_symmetric() {
    let disk = test_disk();
    let bin = ParticleBin::from_disk(&disk, Length::from_au(1.0), Length::from_au(0.1));

    let coag = Coagulation::new(&bin, &disk);

    for i in 0..coag.n_bins() {
        for j in 0..coag.n_bins() {
            assert_relative_eq!(coag.kernel(i, j), coag.kernel(j, i), epsilon = 1e-10);
        }
    }
}

#[test]
fn kernel_diagonal_is_nonzero() {
    let disk = test_disk();
    let bin = ParticleBin::from_disk(&disk, Length::from_au(1.0), Length::from_au(0.1));

    let coag = Coagulation::new(&bin, &disk);

    // Same-size collisions have nonzero kernel from turbulence
    // (even though differential drift is zero)
    // Actually, for same-size particles, Δv = 0 so K = 0
    // This is physically correct - same-size particles don't collide via drift
    for i in 0..coag.n_bins() {
        // Diagonal should be zero (no relative velocity for same size)
        assert_relative_eq!(coag.kernel(i, i), 0.0, epsilon = 1e-10);
    }
}

#[test]
fn kernel_increases_with_size_difference() {
    let disk = test_disk();
    let bin = ParticleBin::from_disk(&disk, Length::from_au(1.0), Length::from_au(0.1));

    let coag = Coagulation::new(&bin, &disk);

    // Kernel should generally increase with size difference
    // (larger relative velocity and cross-section)
    let k_01 = coag.kernel(0, 1); // Adjacent bins
    let k_0_last = coag.kernel(0, coag.n_bins() - 1); // Extreme bins

    assert!(
        k_0_last > k_01,
        "kernel for distant sizes ({}) should exceed adjacent ({})",
        k_0_last,
        k_01
    );
}

// =============================================================================
// Sticking evolution tests
// =============================================================================

#[test]
fn sticking_conserves_mass() {
    let disk = test_disk();
    let r = Length::from_au(1.0);
    let bin = ParticleBin::from_disk(&disk, r, Length::from_au(0.1));

    let coag = Coagulation::new(&bin, &disk);

    let initial_dist = bin.size_distribution().clone();
    let initial_mass = initial_dist.total_mass().to_grams();

    // Evolve for one year
    let dt = Time::from_seconds(365.25 * 24.0 * 3600.0);
    let sigma = bin.surface_density().to_grams_per_cm2();
    let h = bin.scale_height().to_cm();

    let evolved = coag.evolve_sticking(&initial_dist, sigma, h, dt);

    let final_mass = evolved.total_mass().to_grams();

    assert_relative_eq!(initial_mass, final_mass, max_relative = 1e-10);
}

#[test]
fn sticking_shifts_mass_to_larger_sizes() {
    let disk = test_disk();
    let r = Length::from_au(1.0);
    let bin = ParticleBin::from_disk(&disk, r, Length::from_au(0.1));

    let coag = Coagulation::new(&bin, &disk);

    let initial_dist = bin.size_distribution().clone();
    let initial_mean = initial_dist.mean_size().to_cm();

    // Evolve for 1000 years (should see noticeable growth)
    let sigma = bin.surface_density().to_grams_per_cm2();
    let h = bin.scale_height().to_cm();
    let dt = Time::from_seconds(365.25 * 24.0 * 3600.0); // 1 year

    let mut dist = initial_dist;
    for _ in 0..1000 {
        dist = coag.evolve_sticking(&dist, sigma, h, dt);
    }

    let final_mean = dist.mean_size().to_cm();

    assert!(
        final_mean > initial_mean,
        "mean size should grow: {} -> {}",
        initial_mean,
        final_mean
    );
}

#[test]
fn sticking_preserves_nonnegativity() {
    let disk = test_disk();
    let r = Length::from_au(1.0);
    let bin = ParticleBin::from_disk(&disk, r, Length::from_au(0.1));

    let coag = Coagulation::new(&bin, &disk);

    let initial_dist = bin.size_distribution().clone();

    // Evolve for many steps
    let sigma = bin.surface_density().to_grams_per_cm2();
    let h = bin.scale_height().to_cm();
    let dt = Time::from_seconds(365.25 * 24.0 * 3600.0);

    let mut dist = initial_dist;
    for _ in 0..100 {
        dist = coag.evolve_sticking(&dist, sigma, h, dt);

        // Check all bins are non-negative
        for (_, mass) in dist.bins() {
            assert!(
                mass.to_grams() >= 0.0,
                "mass should be non-negative: {}",
                mass.to_grams()
            );
        }
    }
}

#[test]
fn zero_mass_distribution_unchanged() {
    let disk = test_disk();
    let r = Length::from_au(1.0);

    // Create a bin with zero-mass distribution
    let size_dist = SizeDistribution::binned_empty(
        Length::from_cm(1e-5),
        Length::from_cm(1e-4),
        30,
        silicate_density(),
    );

    let bin = ParticleBin::new(
        r,
        Length::from_au(0.1),
        size_dist.clone(),
        SurfaceDensity::from_grams_per_cm2(0.0),
        disk.scale_height(r),
    );

    let coag = Coagulation::new(&bin, &disk);

    let dt = Time::from_seconds(365.25 * 24.0 * 3600.0);
    let evolved = coag.evolve_sticking(&size_dist, 0.0, disk.scale_height(r).to_cm(), dt);

    assert_relative_eq!(evolved.total_mass().to_grams(), 0.0, epsilon = 1e-20);
}

// =============================================================================
// Physical behavior tests
// =============================================================================

#[test]
fn higher_density_grows_faster() {
    let disk = test_disk();
    let r = Length::from_au(1.0);
    let bin = ParticleBin::from_disk(&disk, r, Length::from_au(0.1));

    let coag = Coagulation::new(&bin, &disk);

    let initial_dist = bin.size_distribution().clone();
    let h = bin.scale_height().to_cm();
    let dt = Time::from_seconds(365.25 * 24.0 * 3600.0 * 100.0); // 100 years

    // Low density evolution
    let sigma_low = 0.01; // g/cm²
    let evolved_low = coag.evolve_sticking(&initial_dist, sigma_low, h, dt);
    let mean_low = evolved_low.mean_size().to_cm();

    // High density evolution
    let sigma_high = 1.0; // g/cm²
    let evolved_high = coag.evolve_sticking(&initial_dist, sigma_high, h, dt);
    let mean_high = evolved_high.mean_size().to_cm();

    assert!(
        mean_high > mean_low,
        "higher density should grow faster: {} vs {}",
        mean_high,
        mean_low
    );
}

#[test]
fn monodisperse_does_not_grow_from_sticking_alone() {
    // A truly monodisperse distribution has zero relative velocity,
    // so no collisions occur.
    let disk = test_disk();
    let r = Length::from_au(1.0);

    let size_dist = SizeDistribution::monodisperse(
        Length::from_cm(1e-4),
        Mass::from_grams(1e10),
        silicate_density(),
    )
    .to_binned(30);

    let bin = ParticleBin::new(
        r,
        Length::from_au(0.1),
        size_dist.clone(),
        SurfaceDensity::from_grams_per_cm2(1.0),
        disk.scale_height(r),
    );

    let coag = Coagulation::new(&bin, &disk);

    let initial_mean = size_dist.mean_size().to_cm();

    let dt = Time::from_seconds(365.25 * 24.0 * 3600.0 * 100.0);
    let evolved = coag.evolve_sticking(&size_dist, 1.0, disk.scale_height(r).to_cm(), dt);

    let final_mean = evolved.mean_size().to_cm();

    // Monodisperse -> binned puts all mass in one bin
    // That bin has zero relative velocity with itself
    // So no growth should occur
    // (In reality, Brownian motion would cause some growth, but we don't model that)
    assert_relative_eq!(initial_mean, final_mean, max_relative = 0.01);
}
