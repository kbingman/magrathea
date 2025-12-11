use nalgebra::{Point2, Vector2};
use stellar::generation::main_sequence_star;

use crate::forces::tree_gravity::TreeGravity;
use crate::forces::{DirectGravity, ForceModel};
use crate::state::SystemState;

#[test]
fn test_tree_acceleration_toward_star() {
    let star = main_sequence_star(1.0, 0.0, 4_600.0);
    let mut system = SystemState::new(star);
    system.add_body(0.001, 0.01, Point2::new(1.0, 0.0), Vector2::new(0.0, 6.28));

    let gravity = TreeGravity::new();
    let accel = gravity.acceleration(0, &system);

    // Should point toward star
    assert!(accel.x < 0.0);
    assert!(accel.y.abs() < 1e-10);
}

#[test]
fn test_tree_vs_direct_single_body() {
    let star = main_sequence_star(1.0, 0.0, 4_600.0);
    let mut system = SystemState::new(star);
    system.add_body(0.001, 0.01, Point2::new(1.0, 0.0), Vector2::new(0.0, 0.0));

    let tree_gravity = TreeGravity::new();
    let direct_gravity = DirectGravity::new();

    let accel_tree = tree_gravity.acceleration(0, &system);
    let accel_direct = direct_gravity.acceleration(0, &system);

    // For single body, should be identical (only star force)
    let error = (accel_tree - accel_direct).magnitude() / accel_direct.magnitude();
    assert!(error < 1e-10);
}

#[test]
fn test_tree_vs_direct_two_bodies() {
    let star = main_sequence_star(1.0, 0.0, 4_600.0);
    let mut system = SystemState::new(star);
    system.add_body(0.001, 0.01, Point2::new(1.0, 0.0), Vector2::new(0.0, 6.28));
    system.add_body(0.001, 0.01, Point2::new(2.0, 0.0), Vector2::new(0.0, 4.44));

    let tree_gravity = TreeGravity::with_theta(0.1); // Very accurate
    let direct_gravity = DirectGravity::new();

    let accel_tree = tree_gravity.acceleration(0, &system);
    let accel_direct = direct_gravity.acceleration(0, &system);

    // Should be very close with small theta
    let error = (accel_tree - accel_direct).magnitude() / accel_direct.magnitude();
    assert!(error < 0.01); // Within 1%
}

#[test]
fn test_theta_affects_accuracy() {
    let star = main_sequence_star(1.0, 0.0, 4_600.0);
    let mut system = SystemState::new(star);

    // Add several bodies
    for i in 0..5 {
        let r = 1.0 + i as f64 * 0.5;
        system.add_body(0.001, 0.01, Point2::new(r, 0.0), Vector2::new(0.0, 0.0));
    }

    let direct_gravity = DirectGravity::new();
    let accel_exact = direct_gravity.acceleration(0, &system);

    // Smaller theta should be more accurate
    let tree_accurate = TreeGravity::with_theta(0.1);
    let tree_fast = TreeGravity::with_theta(1.0);

    let accel_accurate = tree_accurate.acceleration(0, &system);
    let accel_fast = tree_fast.acceleration(0, &system);

    let error_accurate = (accel_accurate - accel_exact).magnitude() / accel_exact.magnitude();
    let error_fast = (accel_fast - accel_exact).magnitude() / accel_exact.magnitude();

    // Smaller theta = smaller error
    assert!(error_accurate < error_fast);
}

#[test]
fn test_potential_energy_negative() {
    let star = main_sequence_star(1.0, 0.0, 4_600.0);
    let mut system = SystemState::new(star);
    system.add_body(0.001, 0.01, Point2::new(1.0, 0.0), Vector2::new(0.0, 6.28));

    let gravity = TreeGravity::new();
    let pe = gravity.potential_energy(&system);

    assert!(pe < 0.0);
}

#[test]
fn test_tree_potential_matches_direct() {
    let star = main_sequence_star(1.0, 0.0, 4_600.0);
    let mut system = SystemState::new(star);
    system.add_body(0.001, 0.01, Point2::new(1.0, 0.0), Vector2::new(0.0, 6.28));
    system.add_body(0.001, 0.01, Point2::new(2.0, 0.0), Vector2::new(0.0, 4.44));

    let tree_gravity = TreeGravity::new();
    let direct_gravity = DirectGravity::new();

    let pe_tree = tree_gravity.potential_energy(&system);
    let pe_direct = direct_gravity.potential_energy(&system);

    // Potential uses direct calculation, should match exactly
    let error = (pe_tree - pe_direct).abs() / pe_direct.abs();
    assert!(error < 1e-10);
}

#[test]
fn test_many_bodies() {
    let star = main_sequence_star(1.0, 0.0, 4_600.0);
    let mut system = SystemState::new(star);

    // Add many bodies in a disk
    for i in 0..50 {
        let angle = (i as f64) * 2.0 * std::f64::consts::PI / 50.0;
        let r = 1.0 + (i as f64) * 0.1;
        let x = r * angle.cos();
        let y = r * angle.sin();
        system.add_body(0.0001, 0.001, Point2::new(x, y), Vector2::new(-y, x));
    }

    let gravity = TreeGravity::new();

    // Should compute without panic
    let accel = gravity.acceleration(0, &system);
    assert!(accel.magnitude() > 0.0);
}
