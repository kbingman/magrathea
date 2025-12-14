use nalgebra::{Point2, Vector2};
use stellar::generation::main_sequence_star;

use crate::forces::{DirectGravity, ForceModel, G};
use crate::integrator::{Euler, Integrator, Leapfrog};
use crate::state::SystemState;

fn make_test_system() -> SystemState {
    let star = main_sequence_star(1.0, 0.0, 4_600.0);
    let mut system = SystemState::new(star);

    // Add Earth at 1 AU in circular orbit
    // v_circular = sqrt(GM/r) = sqrt(G*1.0/1.0) ≈ 6.28 AU/year
    let v_circular = (G * 1.0 / 1.0).sqrt();
    system.add_body(
        3.0e-6, // Earth mass in solar masses
        0.01,
        Point2::new(1.0, 0.0),
        Vector2::new(0.0, v_circular),
    );

    system
}

#[test]
fn test_leapfrog_advances_time() {
    let mut system = make_test_system();
    let integrator = Leapfrog::new();
    let force = DirectGravity::new();

    assert_eq!(system.time, 0.0);

    integrator.step(&mut system, 0.01, &force);

    assert!((system.time - 0.01).abs() < 1e-15);
}

#[test]
fn test_leapfrog_changes_position() {
    let mut system = make_test_system();
    let initial_pos = system.bodies[0].position;

    let integrator = Leapfrog::new();
    let force = DirectGravity::new();

    integrator.step(&mut system, 0.01, &force);

    let final_pos = system.bodies[0].position;
    assert!((final_pos - initial_pos).magnitude() > 0.0);
}

#[test]
fn test_leapfrog_circular_orbit() {
    let mut system = make_test_system();
    let initial_r = system.bodies[0].position.coords.magnitude();

    let integrator = Leapfrog::new();
    let force = DirectGravity::new();

    // Integrate for one orbit (2π years for 1 AU orbit)
    let period = 2.0 * std::f64::consts::PI;
    let dt = period / 1000.0; // 1000 steps per orbit for better accuracy
    let n_steps = 1000;

    integrator.integrate(&mut system, dt, n_steps, &force);

    let final_r = system.bodies[0].position.coords.magnitude();

    // Radius should be preserved (within ~0.1% with 1000 steps)
    let error = (final_r - initial_r).abs() / initial_r;
    assert!(error < 0.001, "Radius error: {:.2e}", error);
}

#[test]
fn test_leapfrog_conserves_energy() {
    let mut system = make_test_system();
    let force = DirectGravity::new();

    let initial_ke = system.kinetic_energy();
    let initial_pe = force.potential_energy(&system);
    let initial_energy = initial_ke + initial_pe;

    let integrator = Leapfrog::new();
    let dt = 0.01;

    // Integrate for 100 steps
    integrator.integrate(&mut system, dt, 100, &force);

    let final_ke = system.kinetic_energy();
    let final_pe = force.potential_energy(&system);
    let final_energy = final_ke + final_pe;

    // Energy should be conserved to within 0.1%
    let energy_error = (final_energy - initial_energy).abs() / initial_energy.abs();
    assert!(energy_error < 1e-3, "Energy error: {:.2e}", energy_error);
}

#[test]
fn test_leapfrog_conserves_angular_momentum() {
    let mut system = make_test_system();
    let initial_l = system.total_angular_momentum();

    let integrator = Leapfrog::new();
    let force = DirectGravity::new();
    let dt = 0.01;

    integrator.integrate(&mut system, dt, 100, &force);

    let final_l = system.total_angular_momentum();

    // Angular momentum should be conserved to machine precision
    let error = (final_l - initial_l).abs() / initial_l.abs();
    assert!(error < 1e-10, "Angular momentum error: {:.2e}", error);
}

#[test]
fn test_dkd_vs_kdkd_equivalent_for_gravity() {
    let mut system_kdk = make_test_system();
    let mut system_dkd = make_test_system();

    let integrator_kdk = Leapfrog::new();
    let integrator_dkd = Leapfrog::new_dkd();
    let force = DirectGravity::new();

    let dt = 0.01;
    let n_steps = 10;

    integrator_kdk.integrate(&mut system_kdk, dt, n_steps, &force);
    integrator_dkd.integrate(&mut system_dkd, dt, n_steps, &force);

    // KDK and DKD should give similar results (both are 2nd order accurate)
    // but they differ by O(dt²) due to different ordering
    let pos_diff = (system_kdk.bodies[0].position - system_dkd.bodies[0].position).magnitude();
    let vel_diff = (system_kdk.bodies[0].velocity - system_dkd.bodies[0].velocity).magnitude();

    // Differences should be small (O(dt²))
    assert!(pos_diff < 1e-2, "Position diff: {:.2e}", pos_diff);
    assert!(vel_diff < 1e-2, "Velocity diff: {:.2e}", vel_diff);
}

#[test]
fn test_euler_advances_time() {
    let mut system = make_test_system();
    let force = DirectGravity::new();

    Euler.step(&mut system, 0.01, &force);

    assert!((system.time - 0.01).abs() < 1e-15);
}

#[test]
fn test_euler_worse_than_leapfrog() {
    let mut system_euler = make_test_system();
    let mut system_leapfrog = make_test_system();

    let force = DirectGravity::new();

    let initial_energy = system_euler.kinetic_energy() + force.potential_energy(&system_euler);

    let dt = 0.1; // Large timestep to exaggerate errors
    let n_steps = 10;

    Euler.integrate(&mut system_euler, dt, n_steps, &force);
    Leapfrog::new().integrate(&mut system_leapfrog, dt, n_steps, &force);

    let euler_energy = system_euler.kinetic_energy() + force.potential_energy(&system_euler);
    let leapfrog_energy =
        system_leapfrog.kinetic_energy() + force.potential_energy(&system_leapfrog);

    let euler_error = (euler_energy - initial_energy).abs() / initial_energy.abs();
    let leapfrog_error = (leapfrog_energy - initial_energy).abs() / initial_energy.abs();

    // Euler should have worse energy conservation
    assert!(euler_error > leapfrog_error);
    assert!(euler_error > 0.01); // Euler has significant error
    assert!(leapfrog_error < 0.01); // Leapfrog much better
}

#[test]
fn test_multi_step_integration() {
    let mut system = make_test_system();
    let integrator = Leapfrog::new();
    let force = DirectGravity::new();

    let dt = 0.01;
    let n_steps = 50;

    let final_time = integrator.integrate(&mut system, dt, n_steps, &force);

    let expected_time = dt * (n_steps as f64);
    assert!((final_time - expected_time).abs() < 1e-10);
    assert!((system.time - expected_time).abs() < 1e-10);
}

#[test]
fn test_two_body_system() {
    let star = main_sequence_star(1.0, 0.0, 4_600.0);
    let mut system = SystemState::new(star);

    // Add two planets that interact
    let v1 = (G * 1.0 / 1.0).sqrt();
    let v2 = (G * 1.0 / 2.0).sqrt();

    system.add_body(3.0e-6, 0.01, Point2::new(1.0, 0.0), Vector2::new(0.0, v1));
    system.add_body(3.0e-6, 0.01, Point2::new(2.0, 0.0), Vector2::new(0.0, v2));

    let integrator = Leapfrog::new();
    let force = DirectGravity::new();

    let initial_l = system.total_angular_momentum();

    integrator.integrate(&mut system, 0.01, 100, &force);

    let final_l = system.total_angular_momentum();

    // Angular momentum conserved
    let error = (final_l - initial_l).abs() / initial_l.abs();
    assert!(error < 1e-10);
}

#[test]
fn test_empty_system() {
    let star = main_sequence_star(1.0, 0.0, 4_600.0);
    let mut system = SystemState::new(star);

    let integrator = Leapfrog::new();
    let force = DirectGravity::new();

    // Should not panic with empty system
    integrator.step(&mut system, 0.01, &force);
    assert_eq!(system.bodies.len(), 0);
}

#[test]
fn test_small_timestep_more_accurate() {
    let mut system_small = make_test_system();
    let mut system_large = make_test_system();

    let force = DirectGravity::new();
    let integrator = Leapfrog::new();

    let initial_energy = system_small.kinetic_energy() + force.potential_energy(&system_small);

    let total_time = 1.0;

    // Small timestep
    let dt_small = 0.001;
    let n_small = (total_time / dt_small) as usize;
    integrator.integrate(&mut system_small, dt_small, n_small, &force);

    // Large timestep
    let dt_large = 0.1;
    let n_large = (total_time / dt_large) as usize;
    integrator.integrate(&mut system_large, dt_large, n_large, &force);

    let energy_small = system_small.kinetic_energy() + force.potential_energy(&system_small);
    let energy_large = system_large.kinetic_energy() + force.potential_energy(&system_large);

    let error_small = (energy_small - initial_energy).abs() / initial_energy.abs();
    let error_large = (energy_large - initial_energy).abs() / initial_energy.abs();

    // Smaller timestep should be more accurate
    assert!(error_small < error_large);
}
