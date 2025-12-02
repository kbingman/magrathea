//! Tests for particle bins.

use approx::assert_relative_eq;
use units::{Density, Length, Mass, SurfaceDensity, Time};

use crate::disk::{DiskModel, GasDisk};
use crate::particles::{ParticleBin, SizeDistribution};

fn test_disk() -> GasDisk {
    GasDisk::mmsn()
}

// =============================================================================
// Constructor tests
// =============================================================================

#[test]
fn from_disk_creates_bin_at_correct_radius() {
    let disk = test_disk();
    let r = Length::from_au(1.0);
    let width = Length::from_au(0.1);

    let bin = ParticleBin::from_disk(&disk, r, width);

    assert_relative_eq!(bin.radial_center().to_au(), 1.0, epsilon = 1e-10);
    assert_relative_eq!(bin.radial_width().to_au(), 0.1, epsilon = 1e-10);
}

#[test]
fn from_disk_has_correct_dust_to_gas_ratio() {
    let disk = test_disk();
    let r = Length::from_au(1.0);
    let width = Length::from_au(0.1);

    let bin = ParticleBin::from_disk(&disk, r, width);

    let sigma_gas = disk.surface_density(r).to_grams_per_cm2();
    let sigma_dust = bin.surface_density().to_grams_per_cm2();

    let ratio = sigma_dust / sigma_gas;
    assert_relative_eq!(ratio, 0.01, epsilon = 1e-10);
}

#[test]
fn from_disk_scale_height_matches_gas() {
    let disk = test_disk();
    let r = Length::from_au(1.0);
    let width = Length::from_au(0.1);

    let bin = ParticleBin::from_disk(&disk, r, width);

    let h_gas = disk.scale_height(r).to_cm();
    let h_dust = bin.scale_height().to_cm();

    // Initially well-mixed
    assert_relative_eq!(h_dust, h_gas, epsilon = 1e-10);
}

#[test]
fn new_constructor_stores_parameters() {
    let size_dist =
        SizeDistribution::mrn(Mass::from_grams(1000.0), Density::from_grams_per_cm3(3.0));

    let bin = ParticleBin::new(
        Length::from_au(2.0),
        Length::from_au(0.2),
        size_dist,
        SurfaceDensity::from_grams_per_cm2(10.0),
        Length::from_cm(1e12),
    );

    assert_relative_eq!(bin.radial_center().to_au(), 2.0, epsilon = 1e-10);
    assert_relative_eq!(bin.radial_width().to_au(), 0.2, epsilon = 1e-10);
    assert_relative_eq!(
        bin.surface_density().to_grams_per_cm2(),
        10.0,
        epsilon = 1e-10
    );
    assert_relative_eq!(bin.scale_height().to_cm(), 1e12, epsilon = 1e-10);
}

// =============================================================================
// Geometric queries
// =============================================================================

#[test]
fn inner_outer_radius_correct() {
    let disk = test_disk();
    let r = Length::from_au(1.0);
    let width = Length::from_au(0.2);

    let bin = ParticleBin::from_disk(&disk, r, width);

    assert_relative_eq!(bin.inner_radius().to_au(), 0.9, epsilon = 1e-10);
    assert_relative_eq!(bin.outer_radius().to_au(), 1.1, epsilon = 1e-10);
}

#[test]
fn area_is_annulus_area() {
    let disk = test_disk();
    let r = Length::from_au(1.0);
    let width = Length::from_au(0.2);

    let bin = ParticleBin::from_disk(&disk, r, width);

    let r_in = bin.inner_radius().to_cm();
    let r_out = bin.outer_radius().to_cm();
    let expected_area = std::f64::consts::PI * (r_out.powi(2) - r_in.powi(2));

    // Use relative epsilon for large numbers
    assert_relative_eq!(bin.area(), expected_area, max_relative = 1e-10);
}

// =============================================================================
// Mass queries
// =============================================================================

#[test]
fn total_mass_equals_sigma_times_area() {
    let disk = test_disk();
    let r = Length::from_au(1.0);
    let width = Length::from_au(0.1);

    let bin = ParticleBin::from_disk(&disk, r, width);

    let expected = bin.surface_density().to_grams_per_cm2() * bin.area();
    let actual = bin.total_mass().to_grams();

    assert_relative_eq!(actual, expected, epsilon = 1e-10);
}

#[test]
fn midplane_density_positive() {
    let disk = test_disk();
    let bin = ParticleBin::from_disk(&disk, Length::from_au(1.0), Length::from_au(0.1));

    let rho = bin.midplane_density().to_grams_per_cm3();
    assert!(rho > 0.0);
}

// =============================================================================
// Aerodynamic properties
// =============================================================================

#[test]
fn mean_stokes_number_small_for_dust() {
    let disk = test_disk();
    let bin = ParticleBin::from_disk(&disk, Length::from_au(1.0), Length::from_au(0.1));

    // MRN distribution has μm-sized grains, so τ_s << 1
    let tau = bin.mean_stokes_number(&disk);
    assert!(tau < 1e-3, "τ_s = {} too large for μm grains", tau);
}

#[test]
fn mean_drift_velocity_inward() {
    let disk = test_disk();
    let bin = ParticleBin::from_disk(&disk, Length::from_au(1.0), Length::from_au(0.1));

    let v_r = bin.mean_drift_velocity(&disk).to_cm_per_sec();
    // Drift should be inward (negative)
    assert!(v_r < 0.0, "drift velocity {} should be negative", v_r);
}

#[test]
fn equilibrium_scale_height_less_than_gas_for_large_particles() {
    let disk = test_disk();
    let r = Length::from_au(1.0);
    let width = Length::from_au(0.1);

    // Create a bin with larger particles (mm-sized)
    let size_dist = SizeDistribution::monodisperse(
        Length::from_cm(0.1), // 1 mm
        Mass::from_grams(1e20),
        Density::from_grams_per_cm3(3.0),
    );

    let bin = ParticleBin::new(
        r,
        width,
        size_dist,
        SurfaceDensity::from_grams_per_cm2(1.0),
        disk.scale_height(r),
    );

    let h_eq = bin.equilibrium_scale_height(&disk).to_cm();
    let h_gas = disk.scale_height(r).to_cm();

    // Larger particles should settle to a thinner layer
    assert!(
        h_eq < h_gas,
        "h_eq = {} should be < h_gas = {}",
        h_eq,
        h_gas
    );
}

// =============================================================================
// Scale height relaxation
// =============================================================================

#[test]
fn scale_height_relaxes_toward_equilibrium() {
    let disk = test_disk();
    let r = Length::from_au(1.0);
    let width = Length::from_au(0.1);

    // Create a bin with larger particles that should settle
    let size_dist = SizeDistribution::monodisperse(
        Length::from_cm(0.1), // 1 mm
        Mass::from_grams(1e20),
        Density::from_grams_per_cm3(3.0),
    );

    let mut bin = ParticleBin::new(
        r,
        width,
        size_dist,
        SurfaceDensity::from_grams_per_cm2(1.0),
        disk.scale_height(r), // Start at gas scale height
    );

    let h_initial = bin.scale_height().to_cm();
    let h_eq = bin.equilibrium_scale_height(&disk).to_cm();

    // Evolve for many orbits
    let orbital_period = disk.orbital_period(r).to_seconds();
    let dt = Time::from_seconds(orbital_period * 100.0);
    bin.relax_scale_height(&disk, dt);

    let h_final = bin.scale_height().to_cm();

    // Should have moved toward equilibrium
    assert!(
        (h_final - h_eq).abs() < (h_initial - h_eq).abs(),
        "h should relax toward h_eq: {} -> {} (eq = {})",
        h_initial,
        h_final,
        h_eq
    );
}

#[test]
fn small_particles_stay_well_mixed() {
    let disk = test_disk();
    let r = Length::from_au(1.0);
    let width = Length::from_au(0.1);

    // MRN has μm grains - should stay well-mixed
    let mut bin = ParticleBin::from_disk(&disk, r, width);

    // Evolve for many orbits
    let orbital_period = disk.orbital_period(r).to_seconds();
    let dt = Time::from_seconds(orbital_period * 100.0);
    bin.relax_scale_height(&disk, dt);

    let h_final = bin.scale_height().to_cm();

    // Small particles should stay near gas scale height
    let h_gas = disk.scale_height(r).to_cm();
    assert_relative_eq!(h_final / h_gas, 1.0, epsilon = 0.1);
}

// =============================================================================
// Size distribution access
// =============================================================================

#[test]
fn size_distribution_accessible() {
    let disk = test_disk();
    let bin = ParticleBin::from_disk(&disk, Length::from_au(1.0), Length::from_au(0.1));

    // Should be able to query size distribution
    let mean_s = bin.mean_size().to_cm();
    assert!(mean_s > 0.0);

    // Should be MRN range
    let s_min = bin.size_distribution().min_size().to_cm();
    let s_max = bin.size_distribution().max_size().to_cm();
    assert_relative_eq!(s_min, 1e-5, epsilon = 1e-10); // 0.1 μm
    assert_relative_eq!(s_max, 1e-4, epsilon = 1e-10); // 1 μm
}

#[test]
fn material_density_is_silicate() {
    let disk = test_disk();
    let bin = ParticleBin::from_disk(&disk, Length::from_au(1.0), Length::from_au(0.1));

    assert_relative_eq!(
        bin.material_density().to_grams_per_cm3(),
        3.0,
        epsilon = 1e-10
    );
}

// =============================================================================
// Size-resolved queries
// =============================================================================

#[test]
fn stokes_numbers_returns_one_per_bin() {
    let disk = test_disk();
    let bin = ParticleBin::from_disk(&disk, Length::from_au(1.0), Length::from_au(0.1));

    let stokes = bin.stokes_numbers(&disk);

    // Should have 30 bins (N_SIZE_BINS constant)
    assert_eq!(stokes.len(), 30);

    // All Stokes numbers should be positive
    for (size, tau) in &stokes {
        assert!(tau > &0.0, "τ_s should be positive for size {:?}", size);
    }
}

#[test]
fn stokes_numbers_increase_with_size() {
    let disk = test_disk();
    let bin = ParticleBin::from_disk(&disk, Length::from_au(1.0), Length::from_au(0.1));

    let stokes = bin.stokes_numbers(&disk);

    // Stokes number should increase monotonically with size
    for i in 1..stokes.len() {
        assert!(
            stokes[i].1 > stokes[i - 1].1,
            "τ_s should increase with size: bin {} ({:?}) has τ={}, bin {} ({:?}) has τ={}",
            i - 1,
            stokes[i - 1].0,
            stokes[i - 1].1,
            i,
            stokes[i].0,
            stokes[i].1
        );
    }
}

#[test]
fn drift_velocities_returns_one_per_bin() {
    let disk = test_disk();
    let bin = ParticleBin::from_disk(&disk, Length::from_au(1.0), Length::from_au(0.1));

    let velocities = bin.drift_velocities(&disk);

    assert_eq!(velocities.len(), 30);
}

#[test]
fn drift_velocities_all_inward() {
    let disk = test_disk();
    let bin = ParticleBin::from_disk(&disk, Length::from_au(1.0), Length::from_au(0.1));

    let velocities = bin.drift_velocities(&disk);

    // All drift should be inward (negative)
    for (size, v) in &velocities {
        assert!(
            v.to_cm_per_sec() < 0.0,
            "drift should be inward for size {:?}",
            size
        );
    }
}

#[test]
fn relative_velocity_zero_for_same_size() {
    let disk = test_disk();
    let bin = ParticleBin::from_disk(&disk, Length::from_au(1.0), Length::from_au(0.1));

    let s = Length::from_cm(1e-5);
    let dv = bin.relative_velocity(&disk, s, s);

    // Same-size particles have zero differential drift
    // (turbulent component is also zero for equal τ_s)
    assert_relative_eq!(dv.to_cm_per_sec(), 0.0, epsilon = 1e-10);
}

#[test]
fn relative_velocity_nonzero_for_different_sizes() {
    let disk = test_disk();
    let bin = ParticleBin::from_disk(&disk, Length::from_au(1.0), Length::from_au(0.1));

    let s1 = Length::from_cm(1e-5); // 0.1 μm
    let s2 = Length::from_cm(1e-4); // 1 μm
    let dv = bin.relative_velocity(&disk, s1, s2);

    // Different sizes should have nonzero relative velocity
    assert!(
        dv.to_cm_per_sec() > 0.0,
        "relative velocity should be > 0 for different sizes"
    );
}

#[test]
fn relative_velocity_symmetric() {
    let disk = test_disk();
    let bin = ParticleBin::from_disk(&disk, Length::from_au(1.0), Length::from_au(0.1));

    let s1 = Length::from_cm(1e-5);
    let s2 = Length::from_cm(1e-4);

    let dv_12 = bin.relative_velocity(&disk, s1, s2);
    let dv_21 = bin.relative_velocity(&disk, s2, s1);

    // Should be symmetric: |v1 - v2| = |v2 - v1|
    assert_relative_eq!(
        dv_12.to_cm_per_sec(),
        dv_21.to_cm_per_sec(),
        epsilon = 1e-10
    );
}

#[test]
fn relative_velocity_larger_for_larger_size_difference() {
    let disk = test_disk();
    let bin = ParticleBin::from_disk(&disk, Length::from_au(1.0), Length::from_au(0.1));

    let s_small = Length::from_cm(1e-5);
    let s_medium = Length::from_cm(3e-5);
    let s_large = Length::from_cm(1e-4);

    let dv_small_diff = bin.relative_velocity(&disk, s_small, s_medium);
    let dv_large_diff = bin.relative_velocity(&disk, s_small, s_large);

    // Larger size difference should give larger relative velocity
    assert!(
        dv_large_diff.to_cm_per_sec() > dv_small_diff.to_cm_per_sec(),
        "larger size difference should give larger Δv: {} vs {}",
        dv_large_diff.to_cm_per_sec(),
        dv_small_diff.to_cm_per_sec()
    );
}
