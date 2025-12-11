use nalgebra::{Point2, Vector2};
use stellar::generation::main_sequence_star;

use crate::forces::gravity::{DirectGravity, hill_radius};
use crate::forces::{ForceModel, G};
use crate::state::SystemState;

#[test]
fn test_acceleration_toward_star() {
    let star = main_sequence_star(1.0, 0.0, 4_600.0);
    let mut system = SystemState::new(star);
    system.add_body(0.001, 0.01, Point2::new(1.0, 0.0), Vector2::new(0.0, 6.28));

    let gravity = DirectGravity::new();
    let accel = gravity.acceleration(0, &system);

    // Should point toward star (negative x)
    assert!(accel.x < 0.0);
    assert!(accel.y.abs() < 1e-10);
}

#[test]
fn test_acceleration_magnitude() {
    let star = main_sequence_star(1.0, 0.0, 4_600.0);
    let mut system = SystemState::new(star);
    system.add_body(0.001, 0.01, Point2::new(1.0, 0.0), Vector2::new(0.0, 0.0));

    let gravity = DirectGravity::new();
    let accel = gravity.acceleration(0, &system);

    // For circular orbit: a = v²/r = GM/r²
    // At 1 AU: a = G * 1.0 / 1.0² = G
    let expected = G;
    assert!((accel.magnitude() - expected).abs() / expected < 1e-10);
}

#[test]
fn test_two_body_acceleration() {
    let star = main_sequence_star(1.0, 0.0, 4_600.0);
    let mut system = SystemState::new(star);
    system.add_body(0.001, 0.01, Point2::new(1.0, 0.0), Vector2::new(0.0, 0.0));
    system.add_body(0.001, 0.01, Point2::new(2.0, 0.0), Vector2::new(0.0, 0.0));

    let gravity = DirectGravity::new();

    // Body 0 should feel pull from both star and body 1
    // Body 1 is farther out, so it pulls body 0 away from star (positive x)
    let accel0 = gravity.acceleration(0, &system);
    let accel_star_only = -G * 1.0 / 1.0_f64.powi(2);

    // Should be less negative (weaker toward star) due to body 1 pulling outward
    assert!(accel0.x > accel_star_only);

    // But still negative (net pull toward star)
    assert!(accel0.x < 0.0);
}

#[test]
fn test_softening_reduces_force() {
    let star = main_sequence_star(1.0, 0.0, 4_600.0);
    let mut system = SystemState::new(star);
    system.add_body(0.001, 0.01, Point2::new(0.01, 0.0), Vector2::new(0.0, 0.0));

    let gravity_hard = DirectGravity::new();
    let gravity_soft = DirectGravity::with_softening(0.01);

    let accel_hard = gravity_hard.acceleration(0, &system);
    let accel_soft = gravity_soft.acceleration(0, &system);

    // Softening should reduce force magnitude
    assert!(accel_soft.magnitude() < accel_hard.magnitude());
}

#[test]
fn test_hill_radius_earth() {
    use units::Mass;

    let earth_mass = Mass::from_earth_masses(1.0).to_solar_masses();
    let r_hill = hill_radius(earth_mass, 1.0, 1.0);

    // Earth's Hill radius is about 0.01 AU
    assert!((r_hill - 0.01).abs() < 0.001);
}

#[test]
fn test_potential_energy_negative() {
    let star = main_sequence_star(1.0, 0.0, 4_600.0);
    let mut system = SystemState::new(star);
    system.add_body(0.001, 0.01, Point2::new(1.0, 0.0), Vector2::new(0.0, 6.28));

    let gravity = DirectGravity::new();
    let pe = gravity.potential_energy(&system);

    // Gravitational potential is negative
    assert!(pe < 0.0);
}

#[test]
fn test_potential_energy_two_bodies() {
    let star = main_sequence_star(1.0, 0.0, 4_600.0);
    let mut system = SystemState::new(star);
    system.add_body(0.001, 0.01, Point2::new(1.0, 0.0), Vector2::new(0.0, 6.28));

    let gravity = DirectGravity::new();
    let pe_one = gravity.potential_energy(&system);

    // Add second body
    system.add_body(0.001, 0.01, Point2::new(2.0, 0.0), Vector2::new(0.0, 4.44));
    let pe_two = gravity.potential_energy(&system);

    // More bodies = more negative potential
    assert!(pe_two < pe_one);
}
