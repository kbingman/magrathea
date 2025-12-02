//! Tests for particle aerodynamics.

use approx::assert_relative_eq;
use units::{Density, Length};

use crate::disk::{DiskModel, GasDisk};
use crate::particles::{DragRegime, Particle};

/// Create a standard MMSN disk for testing.
fn test_disk() -> GasDisk {
    GasDisk::mmsn()
}

// =============================================================================
// Drag regime tests
// =============================================================================

#[test]
fn small_particles_are_epstein() {
    let disk = test_disk();
    let r = Length::from_au(1.0);

    // 1 μm grain - should be Epstein
    let dust = Particle::new(Length::from_cm(1e-4), Density::from_grams_per_cm3(3.0));

    assert_eq!(dust.drag_regime(&disk, r), DragRegime::Epstein);
}

#[test]
fn large_particles_are_stokes() {
    let disk = test_disk();
    let r = Length::from_au(1.0);

    // Mean free path at 1 AU in MMSN is ~10 cm
    // Particles > 9λ/4 ≈ 22 cm are in Stokes regime
    let boulder = Particle::new(Length::from_cm(100.0), Density::from_grams_per_cm3(3.0));

    assert_eq!(boulder.drag_regime(&disk, r), DragRegime::Stokes);
}

#[test]
fn regime_transition_near_mean_free_path() {
    let disk = test_disk();
    let r = Length::from_au(1.0);

    let mfp = disk.mean_free_path(r).to_cm();
    let transition_size = 9.0 * mfp / 4.0;

    // Just below transition
    let small = Particle::new(
        Length::from_cm(transition_size * 0.9),
        Density::from_grams_per_cm3(3.0),
    );
    assert_eq!(small.drag_regime(&disk, r), DragRegime::Epstein);

    // Just above transition
    let large = Particle::new(
        Length::from_cm(transition_size * 1.1),
        Density::from_grams_per_cm3(3.0),
    );
    assert_eq!(large.drag_regime(&disk, r), DragRegime::Stokes);
}

// =============================================================================
// Stopping time tests
// =============================================================================

#[test]
fn stopping_time_scales_with_size_epstein() {
    let disk = test_disk();
    let r = Length::from_au(1.0);

    // In Epstein regime, t_s ∝ s
    let p1 = Particle::new(Length::from_cm(0.01), Density::from_grams_per_cm3(3.0));
    let p2 = Particle::new(Length::from_cm(0.1), Density::from_grams_per_cm3(3.0));

    let t1 = p1.stopping_time(&disk, r).to_seconds();
    let t2 = p2.stopping_time(&disk, r).to_seconds();

    // 10x larger particle should have 10x longer stopping time
    assert_relative_eq!(t2 / t1, 10.0, epsilon = 0.01);
}

#[test]
fn stopping_time_scales_with_density() {
    let disk = test_disk();
    let r = Length::from_au(1.0);

    // In Epstein regime, t_s ∝ ρ_m
    let ice = Particle::new(Length::from_cm(0.1), Density::from_grams_per_cm3(1.0));
    let rock = Particle::new(Length::from_cm(0.1), Density::from_grams_per_cm3(3.0));

    let t_ice = ice.stopping_time(&disk, r).to_seconds();
    let t_rock = rock.stopping_time(&disk, r).to_seconds();

    // 3x denser particle should have 3x longer stopping time
    assert_relative_eq!(t_rock / t_ice, 3.0, epsilon = 0.01);
}

#[test]
fn stopping_time_reasonable_values() {
    let disk = test_disk();
    let r = Length::from_au(1.0);

    // 1 mm silicate grain at 1 AU
    let pebble = Particle::new(Length::from_cm(0.1), Density::from_grams_per_cm3(3.0));
    let t_s = pebble.stopping_time(&disk, r).to_seconds();

    // Should be on order of minutes to days
    // Actual value ~1400s (~23 min) for MMSN parameters
    let minute = 60.0;
    let day = 86400.0;
    assert!(t_s > minute, "stopping time {} s too short", t_s);
    assert!(t_s < 100.0 * day, "stopping time {} s too long", t_s);
}

// =============================================================================
// Stokes number tests
// =============================================================================

#[test]
fn stokes_number_unity_exists() {
    let disk = test_disk();
    let r = Length::from_au(1.0);
    let rho_m = Density::from_grams_per_cm3(3.0);

    // Search for particle size where τ_s ≈ 1 by scanning sizes
    // This validates the physics produces reasonable τ_s = 1 sizes
    let sizes = [0.001, 0.01, 0.1, 1.0, 10.0, 100.0, 1000.0];

    let stokes_numbers: Vec<(f64, f64)> = sizes
        .iter()
        .map(|&s| {
            let p = Particle::new(Length::from_cm(s), rho_m);
            (s, p.stokes_number(&disk, r))
        })
        .collect();

    // Find sizes that bracket τ_s = 1
    let below_one = stokes_numbers.iter().filter(|(_, tau)| *tau < 1.0).count();
    let above_one = stokes_numbers.iter().filter(|(_, tau)| *tau > 1.0).count();

    // Should have some particles with τ_s < 1 and some with τ_s > 1
    assert!(
        below_one > 0,
        "no particles with τ_s < 1: {:?}",
        stokes_numbers
    );
    assert!(
        above_one > 0,
        "no particles with τ_s > 1: {:?}",
        stokes_numbers
    );
}

#[test]
fn small_particles_have_small_stokes_number() {
    let disk = test_disk();
    let r = Length::from_au(1.0);

    // 1 μm grain
    let dust = Particle::new(Length::from_cm(1e-4), Density::from_grams_per_cm3(3.0));
    let tau = dust.stokes_number(&disk, r);

    assert!(tau < 1e-3, "τ_s = {} too large for μm grain", tau);
}

#[test]
fn stokes_number_proportional_to_size() {
    let disk = test_disk();
    let r = Length::from_au(1.0);

    // In Epstein regime, τ_s ∝ s
    let p1 = Particle::new(Length::from_cm(0.01), Density::from_grams_per_cm3(3.0));
    let p2 = Particle::new(Length::from_cm(0.1), Density::from_grams_per_cm3(3.0));

    let tau1 = p1.stokes_number(&disk, r);
    let tau2 = p2.stokes_number(&disk, r);

    assert_relative_eq!(tau2 / tau1, 10.0, epsilon = 0.01);
}

// =============================================================================
// Radial drift tests
// =============================================================================

#[test]
fn drift_peaks_at_stokes_one() {
    let disk = test_disk();
    let r = Length::from_au(1.0);

    // Create particles spanning τ_s << 1 to τ_s >> 1
    let sizes = [1e-4, 1e-3, 1e-2, 0.1, 1.0, 10.0, 100.0];
    let rho_m = Density::from_grams_per_cm3(3.0);

    let drifts: Vec<(f64, f64)> = sizes
        .iter()
        .map(|&s| {
            let p = Particle::new(Length::from_cm(s), rho_m);
            let tau = p.stokes_number(&disk, r);
            let v_r = p.radial_drift_velocity(&disk, r).to_cm_per_sec().abs();
            (tau, v_r)
        })
        .collect();

    // Find the maximum drift velocity
    let (tau_max, v_max) = drifts
        .iter()
        .max_by(|a, b| a.1.partial_cmp(&b.1).unwrap())
        .unwrap();

    // Peak should occur near τ_s = 1
    assert!(
        *tau_max > 0.1 && *tau_max < 10.0,
        "peak drift at τ_s = {}, expected ~1",
        tau_max
    );

    // Peak drift should be significant (tens of m/s)
    let v_max_m_s = v_max / 100.0; // cm/s to m/s
    assert!(v_max_m_s > 10.0, "peak drift {} m/s too slow", v_max_m_s);
    assert!(v_max_m_s < 200.0, "peak drift {} m/s too fast", v_max_m_s);
}

#[test]
fn drift_is_inward() {
    let disk = test_disk();
    let r = Length::from_au(1.0);

    let pebble = Particle::new(Length::from_cm(0.1), Density::from_grams_per_cm3(3.0));
    let v_r = pebble.radial_drift_velocity(&disk, r).to_cm_per_sec();

    // Drift should be inward (negative)
    assert!(v_r < 0.0, "drift velocity {} should be negative", v_r);
}

#[test]
fn small_particles_drift_slowly() {
    let disk = test_disk();
    let r = Length::from_au(1.0);

    // μm grain - should be tightly coupled
    let dust = Particle::new(Length::from_cm(1e-4), Density::from_grams_per_cm3(3.0));
    let v_r = dust.radial_drift_velocity(&disk, r).to_cm_per_sec().abs();

    // Should drift much slower than peak
    let v_r_m_s = v_r / 100.0;
    assert!(v_r_m_s < 1.0, "μm grain drifts {} m/s, too fast", v_r_m_s);
}

// =============================================================================
// Settling tests
// =============================================================================

#[test]
fn settling_toward_midplane() {
    let disk = test_disk();
    let r = Length::from_au(1.0);

    let pebble = Particle::new(Length::from_cm(0.1), Density::from_grams_per_cm3(3.0));

    // Above midplane
    let z = disk.scale_height(r);
    let v_z = pebble.settling_velocity(&disk, r, z).to_cm_per_sec();

    // Should settle downward (negative velocity for positive z)
    assert!(v_z < 0.0, "settling velocity {} should be negative", v_z);
}

#[test]
fn settling_scales_with_height() {
    let disk = test_disk();
    let r = Length::from_au(1.0);

    let pebble = Particle::new(Length::from_cm(0.1), Density::from_grams_per_cm3(3.0));

    let h = disk.scale_height(r);
    let v_1h = pebble.settling_velocity(&disk, r, h).to_cm_per_sec().abs();
    let v_2h = pebble
        .settling_velocity(&disk, r, Length::from_cm(2.0 * h.to_cm()))
        .to_cm_per_sec()
        .abs();

    // v_z ∝ z, so 2x height should give 2x settling velocity
    assert_relative_eq!(v_2h / v_1h, 2.0, epsilon = 0.01);
}

// =============================================================================
// Equilibrium scale height tests
// =============================================================================

#[test]
fn small_particles_well_mixed() {
    let disk = test_disk();
    let r = Length::from_au(1.0);

    // μm grain with τ_s << α
    let dust = Particle::new(Length::from_cm(1e-4), Density::from_grams_per_cm3(3.0));
    let h_p = dust.equilibrium_scale_height(&disk, r).to_cm();
    let h_g = disk.scale_height(r).to_cm();

    // Should be close to gas scale height
    assert_relative_eq!(h_p / h_g, 1.0, epsilon = 0.1);
}

#[test]
fn large_particles_settle_thin() {
    let disk = test_disk();
    let r = Length::from_au(1.0);

    // Large particle with τ_s >> α
    // α = 1e-3, so need τ_s >> 0.001
    let boulder = Particle::new(Length::from_cm(10.0), Density::from_grams_per_cm3(3.0));
    let h_p = boulder.equilibrium_scale_height(&disk, r).to_cm();
    let h_g = disk.scale_height(r).to_cm();

    // Should be much thinner than gas
    assert!(h_p / h_g < 0.5, "h_p/h_g = {} should be thin", h_p / h_g);
}
