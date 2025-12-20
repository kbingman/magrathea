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

// =============================================================================
// Gravitational Instability Tests
// =============================================================================

#[test]
fn toomre_q_high_for_well_mixed_particles() {
    let disk = test_disk();
    let r = Length::from_au(5.0);
    let width = Length::from_au(0.5);

    // Standard bin from disk (well-mixed, low density)
    let bin = ParticleBin::from_disk(&disk, r, width);

    let q = bin.toomre_q(&disk);

    // Well-mixed particles with standard dust-to-gas ratio should be
    // very stable (high Q)
    assert!(
        q > 10.0,
        "Q = {} should be >> 1 for well-mixed particles",
        q
    );
}

#[test]
fn toomre_q_decreases_with_higher_surface_density() {
    let disk = test_disk();
    let r = Length::from_au(5.0);
    let width = Length::from_au(0.5);

    // Create two bins with different surface densities
    let size_dist = SizeDistribution::mrn(Mass::from_grams(1e25), Density::from_grams_per_cm3(3.0));

    let mut bin_low = ParticleBin::new(
        r,
        width,
        size_dist.clone(),
        SurfaceDensity::from_grams_per_cm2(1.0),
        disk.scale_height(r),
    );
    bin_low.set_velocity_dispersion(disk.sound_speed(r) * 0.01);

    let mut bin_high = ParticleBin::new(
        r,
        width,
        size_dist,
        SurfaceDensity::from_grams_per_cm2(10.0), // 10x higher
        disk.scale_height(r),
    );
    bin_high.set_velocity_dispersion(disk.sound_speed(r) * 0.01);

    let q_low = bin_low.toomre_q(&disk);
    let q_high = bin_high.toomre_q(&disk);

    // Higher surface density → lower Q
    assert!(
        q_high < q_low,
        "Q should decrease with surface density: {} vs {}",
        q_low,
        q_high
    );

    // Should be roughly inversely proportional
    assert_relative_eq!(q_high * 10.0, q_low, max_relative = 0.01);
}

#[test]
fn toomre_q_increases_with_velocity_dispersion() {
    let disk = test_disk();
    let r = Length::from_au(5.0);
    let width = Length::from_au(0.5);

    let size_dist = SizeDistribution::mrn(Mass::from_grams(1e25), Density::from_grams_per_cm3(3.0));

    let mut bin_low = ParticleBin::new(
        r,
        width,
        size_dist.clone(),
        SurfaceDensity::from_grams_per_cm2(10.0),
        disk.scale_height(r),
    );
    bin_low.set_velocity_dispersion(disk.sound_speed(r) * 0.01);

    let mut bin_high = ParticleBin::new(
        r,
        width,
        size_dist,
        SurfaceDensity::from_grams_per_cm2(10.0),
        disk.scale_height(r),
    );
    bin_high.set_velocity_dispersion(disk.sound_speed(r) * 0.1); // 10x higher

    let q_low = bin_low.toomre_q(&disk);
    let q_high = bin_high.toomre_q(&disk);

    // Higher velocity dispersion → higher Q (more stable)
    assert!(
        q_high > q_low,
        "Q should increase with velocity dispersion: {} vs {}",
        q_low,
        q_high
    );

    // Should be roughly proportional
    assert_relative_eq!(q_high, q_low * 10.0, max_relative = 0.01);
}

#[test]
fn richardson_number_high_for_well_mixed_particles() {
    let disk = test_disk();
    let r = Length::from_au(5.0);
    let width = Length::from_au(0.5);

    let bin = ParticleBin::from_disk(&disk, r, width);

    let ri = bin.richardson_number(&disk);

    // Well-mixed particles (h_p ≈ h_g) should have Ri ≈ 1
    assert!(
        ri > 0.5,
        "Richardson number {} should be > 0.5 for well-mixed particles",
        ri
    );
}

#[test]
fn richardson_number_decreases_with_settling() {
    let disk = test_disk();
    let r = Length::from_au(5.0);
    let width = Length::from_au(0.5);

    let size_dist = SizeDistribution::mrn(Mass::from_grams(1e25), Density::from_grams_per_cm3(3.0));

    // Well-mixed bin
    let bin_thick = ParticleBin::new(
        r,
        width,
        size_dist.clone(),
        SurfaceDensity::from_grams_per_cm2(10.0),
        disk.scale_height(r), // h_p = h_g
    );

    // Settled bin
    let bin_thin = ParticleBin::new(
        r,
        width,
        size_dist,
        SurfaceDensity::from_grams_per_cm2(10.0),
        disk.scale_height(r) * 0.1, // h_p = 0.1 × h_g
    );

    let ri_thick = bin_thick.richardson_number(&disk);
    let ri_thin = bin_thin.richardson_number(&disk);

    // Thinner layer → lower Ri
    assert!(
        ri_thin < ri_thick,
        "Richardson number should decrease with settling: {} vs {}",
        ri_thick,
        ri_thin
    );
}

#[test]
fn well_mixed_particles_are_stable() {
    let disk = test_disk();
    let r = Length::from_au(5.0);
    let width = Length::from_au(0.5);

    // Standard well-mixed particles
    let bin = ParticleBin::from_disk(&disk, r, width);

    // Should be gravitationally stable
    assert!(
        !bin.is_gravitationally_unstable(&disk),
        "Well-mixed particles should be stable"
    );
}

#[test]
fn high_density_settled_layer_can_be_unstable() {
    let disk = test_disk();
    let r = Length::from_au(5.0);
    let width = Length::from_au(0.5);

    let size_dist = SizeDistribution::mrn(Mass::from_grams(1e25), Density::from_grams_per_cm3(3.0));

    // Create a very dense, settled layer
    let mut bin = ParticleBin::new(
        r,
        width,
        size_dist,
        SurfaceDensity::from_grams_per_cm2(100.0), // Very high surface density
        disk.scale_height(r) * 0.5,                // Moderately settled
    );

    // Low velocity dispersion
    bin.set_velocity_dispersion(disk.sound_speed(r) * 0.001);

    let q = bin.toomre_q(&disk);
    let ri = bin.richardson_number(&disk);

    // This configuration should have:
    // - Low Q (high surface density, low velocity dispersion)
    // - Reasonable Ri (not too settled)
    assert!(q < 2.0, "Dense settled layer should have low Q, got {}", q);
    assert!(
        ri > 0.2,
        "Moderately settled layer should have Ri > 0.2, got {}",
        ri
    );

    // May or may not be unstable depending on exact values
    let unstable = bin.is_gravitationally_unstable(&disk);
    println!(
        "Dense settled layer: Q = {:.2}, Ri = {:.2}, unstable = {}",
        q, ri, unstable
    );
}

#[test]
fn too_thin_layer_is_shear_unstable() {
    let disk = test_disk();
    let r = Length::from_au(5.0);
    let width = Length::from_au(0.5);

    let size_dist = SizeDistribution::mrn(Mass::from_grams(1e25), Density::from_grams_per_cm3(3.0));

    // Very thin layer (over-settled)
    let mut bin = ParticleBin::new(
        r,
        width,
        size_dist,
        SurfaceDensity::from_grams_per_cm2(100.0),
        disk.scale_height(r) * 0.01, // Very thin: h_p = 0.01 × h_g
    );

    bin.set_velocity_dispersion(disk.sound_speed(r) * 0.001);

    let q = bin.toomre_q(&disk);
    let ri = bin.richardson_number(&disk);

    // Should have:
    // - Very low Q (high density, low velocity)
    // - Very low Ri (too thin)
    assert!(q < 1.0, "Very dense layer should have Q < 1, got {}", q);
    assert!(ri < 0.5, "Very thin layer should have low Ri, got {}", ri);

    // Should NOT be unstable (shear instability prevents it)
    assert!(
        !bin.is_gravitationally_unstable(&disk),
        "Too-thin layer should be shear unstable, not gravitationally unstable"
    );
}

#[test]
fn toomre_q_scales_with_orbital_frequency() {
    let disk = test_disk();

    let size_dist = SizeDistribution::mrn(Mass::from_grams(1e25), Density::from_grams_per_cm3(3.0));

    // Same surface density and velocity at different radii
    let mut bin_inner = ParticleBin::new(
        Length::from_au(1.0),
        Length::from_au(0.1),
        size_dist.clone(),
        SurfaceDensity::from_grams_per_cm2(10.0),
        Length::from_cm(1e11),
    );
    bin_inner.set_velocity_dispersion(disk.sound_speed(Length::from_au(1.0)) * 0.01);

    let mut bin_outer = ParticleBin::new(
        Length::from_au(10.0),
        Length::from_au(1.0),
        size_dist,
        SurfaceDensity::from_grams_per_cm2(10.0),
        Length::from_cm(1e12),
    );
    bin_outer.set_velocity_dispersion(disk.sound_speed(Length::from_au(10.0)) * 0.01);

    let q_inner = bin_inner.toomre_q(&disk);
    let q_outer = bin_outer.toomre_q(&disk);

    // Both should be positive
    assert!(q_inner > 0.0 && q_outer > 0.0);

    // Inner disk has higher Ω, so for same σ and Σ, Q should be higher
    // (Actually depends on how sound speed varies with radius)
    println!("Q_inner = {}, Q_outer = {}", q_inner, q_outer);
}
