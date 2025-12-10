use nalgebra::{Point2, Vector2};
use stellar::generation::main_sequence_star;

use crate::state::SystemState;

fn make_test_star() -> stellar::stellar_objects::MainSequenceStar {
    main_sequence_star(1.0, 0.0, 4_600.0) // 4.6 Gyr in Myr
}

#[test]
fn test_new_system() {
    let star = make_test_star();
    let system = SystemState::new(star);

    assert_eq!(system.time, 0.0);
    assert_eq!(system.star.mass.to_solar_masses(), 1.0);
    assert_eq!(system.body_count(), 0);
}

#[test]
fn test_add_body() {
    let star = make_test_star();
    let mut system = SystemState::new(star);

    let id1 = system.add_body(0.001, 0.01, Point2::new(1.0, 0.0), Vector2::new(0.0, 6.28));
    let id2 = system.add_body(0.002, 0.02, Point2::new(2.0, 0.0), Vector2::new(0.0, 4.44));

    assert_eq!(system.body_count(), 2);
    assert_eq!(id1.0, 0);
    assert_eq!(id2.0, 1);
}

#[test]
fn test_remove_body() {
    let star = make_test_star();
    let mut system = SystemState::new(star);

    let id = system.add_body(0.001, 0.01, Point2::new(1.0, 0.0), Vector2::new(0.0, 6.28));
    assert_eq!(system.body_count(), 1);

    let removed = system.remove_body(id);
    assert!(removed.is_some());
    assert_eq!(removed.unwrap().id, id);
    assert_eq!(system.body_count(), 0);
}

#[test]
fn test_remove_nonexistent_body() {
    use crate::body::BodyId;

    let star = make_test_star();
    let mut system = SystemState::new(star);

    let removed = system.remove_body(BodyId(999));
    assert!(removed.is_none());
}

#[test]
fn test_get_body() {
    let star = make_test_star();
    let mut system = SystemState::new(star);

    let id = system.add_body(0.001, 0.01, Point2::new(1.0, 0.0), Vector2::new(0.0, 6.28));

    let body = system.get_body(id);
    assert!(body.is_some());
    assert_eq!(body.unwrap().mass, 0.001);
}

#[test]
fn test_get_body_mut() {
    let star = make_test_star();
    let mut system = SystemState::new(star);

    let id = system.add_body(0.001, 0.01, Point2::new(1.0, 0.0), Vector2::new(0.0, 6.28));

    // Modify the body
    if let Some(body) = system.get_body_mut(id) {
        body.mass = 0.002;
    }

    assert_eq!(system.get_body(id).unwrap().mass, 0.002);
}

#[test]
fn test_total_planet_mass() {
    let star = make_test_star();
    let mut system = SystemState::new(star);

    system.add_body(0.001, 0.01, Point2::new(1.0, 0.0), Vector2::new(0.0, 6.28));
    system.add_body(0.002, 0.02, Point2::new(2.0, 0.0), Vector2::new(0.0, 4.44));
    system.add_body(0.003, 0.03, Point2::new(3.0, 0.0), Vector2::new(0.0, 3.63));

    assert_eq!(system.total_planet_mass(), 0.006);
}

#[test]
fn test_total_momentum_zero() {
    let star = make_test_star();
    let mut system = SystemState::new(star);

    // Add two bodies with equal and opposite momentum
    system.add_body(0.001, 0.01, Point2::new(1.0, 0.0), Vector2::new(0.0, 10.0));
    system.add_body(
        0.001,
        0.01,
        Point2::new(-1.0, 0.0),
        Vector2::new(0.0, -10.0),
    );

    let total_p = system.total_momentum();
    assert!(total_p.magnitude() < 1e-10);
}

#[test]
fn test_total_momentum_nonzero() {
    let star = make_test_star();
    let mut system = SystemState::new(star);

    system.add_body(0.001, 0.01, Point2::new(1.0, 0.0), Vector2::new(0.0, 10.0));

    let total_p = system.total_momentum();
    // momentum = mass * velocity = 0.001 * 10.0 = 0.01
    assert!((total_p.y - 0.01).abs() < 1e-10);
}

#[test]
fn test_total_angular_momentum() {
    let star = make_test_star();
    let mut system = SystemState::new(star);

    // Single body in circular orbit
    // L = m * r * v
    let mass = 0.001;
    let r = 1.0;
    let v = 6.28;
    system.add_body(mass, 0.01, Point2::new(r, 0.0), Vector2::new(0.0, v));

    let total_l = system.total_angular_momentum();
    let expected_l = mass * r * v;
    assert!((total_l - expected_l).abs() / expected_l < 1e-10);
}

#[test]
fn test_multiple_bodies_angular_momentum() {
    let star = make_test_star();
    let mut system = SystemState::new(star);

    // Two bodies, both with positive angular momentum
    system.add_body(0.001, 0.01, Point2::new(1.0, 0.0), Vector2::new(0.0, 6.28));
    system.add_body(0.002, 0.02, Point2::new(2.0, 0.0), Vector2::new(0.0, 4.44));

    let l1 = 0.001 * 1.0 * 6.28;
    let l2 = 0.002 * 2.0 * 4.44;
    let expected_total = l1 + l2;

    let total_l = system.total_angular_momentum();
    assert!((total_l - expected_total).abs() / expected_total < 1e-10);
}

#[test]
fn test_system_clone() {
    let star = make_test_star();
    let mut system1 = SystemState::new(star);
    system1.add_body(0.001, 0.01, Point2::new(1.0, 0.0), Vector2::new(0.0, 6.28));

    let system2 = system1.clone();

    assert_eq!(system1.body_count(), system2.body_count());
    assert_eq!(system1.time, system2.time);
    assert_eq!(
        system1.star.mass.to_solar_masses(),
        system2.star.mass.to_solar_masses()
    );
}

#[test]
fn test_time_evolution() {
    let star = make_test_star();
    let mut system = SystemState::new(star);

    assert_eq!(system.time, 0.0);

    system.time = 100.0;
    assert_eq!(system.time, 100.0);
}
