use nalgebra::{Point2, Vector2};
use stellar::generation::main_sequence_star;

use crate::body::{Body, BodyId};
use crate::collisions::CollisionEvent;
use crate::collisions::resolution::*;
use crate::state::SystemState;

#[test]
fn test_merge_bodies_mass_conservation() {
    let a = Body {
        id: BodyId(0),
        mass: 1.0,
        radius: 0.01,
        position: Point2::new(1.0, 0.0),
        velocity: Vector2::new(0.0, 5.0),
    };

    let b = Body {
        id: BodyId(1),
        mass: 2.0,
        radius: 0.01,
        position: Point2::new(1.1, 0.0),
        velocity: Vector2::new(0.0, 3.0),
    };

    let merged = merge_bodies(&a, &b, BodyId(2));

    // Mass should be conserved
    assert!((merged.mass - 3.0).abs() < 1e-10);
}

#[test]
fn test_merge_bodies_momentum_conservation() {
    let a = Body {
        id: BodyId(0),
        mass: 1.0,
        radius: 0.01,
        position: Point2::new(1.0, 0.0),
        velocity: Vector2::new(0.0, 5.0),
    };

    let b = Body {
        id: BodyId(1),
        mass: 2.0,
        radius: 0.01,
        position: Point2::new(1.1, 0.0),
        velocity: Vector2::new(0.0, 3.0),
    };

    let p_initial = a.momentum() + b.momentum();
    let merged = merge_bodies(&a, &b, BodyId(2));
    let p_final = merged.momentum();

    // Momentum should be conserved
    assert!((p_final.x - p_initial.x).abs() < 1e-10);
    assert!((p_final.y - p_initial.y).abs() < 1e-10);
}

#[test]
fn test_merge_bodies_center_of_mass() {
    let a = Body {
        id: BodyId(0),
        mass: 1.0,
        radius: 0.01,
        position: Point2::new(1.0, 0.0),
        velocity: Vector2::new(0.0, 5.0),
    };

    let b = Body {
        id: BodyId(1),
        mass: 1.0,
        radius: 0.01,
        position: Point2::new(2.0, 0.0),
        velocity: Vector2::new(0.0, 3.0),
    };

    let merged = merge_bodies(&a, &b, BodyId(2));

    // Equal masses: COM should be at midpoint
    assert!((merged.position.x - 1.5).abs() < 1e-10);
    assert!(merged.position.y.abs() < 1e-10);
}

#[test]
fn test_merge_bodies_radius() {
    let a = Body {
        id: BodyId(0),
        mass: 1.0,
        radius: 0.01,
        position: Point2::new(1.0, 0.0),
        velocity: Vector2::new(0.0, 5.0),
    };

    let b = Body {
        id: BodyId(1),
        mass: 1.0,
        radius: 0.01,
        position: Point2::new(1.1, 0.0),
        velocity: Vector2::new(0.0, 3.0),
    };

    let merged = merge_bodies(&a, &b, BodyId(2));

    // Volume conservation in 2D: r_new² = r_a² + r_b²
    let expected_radius = (0.01_f64.powi(2) + 0.01_f64.powi(2)).sqrt();
    assert!((merged.radius - expected_radius).abs() < 1e-10);
}

#[test]
fn test_resolve_single_collision() {
    let star = main_sequence_star(1.0, 0.0, 4_600.0);
    let mut system = SystemState::new(star);

    let id_a = system.add_body(0.001, 0.01, Point2::new(1.0, 0.0), Vector2::new(0.0, 6.28));
    let id_b = system.add_body(
        0.001,
        0.01,
        Point2::new(1.001, 0.0),
        Vector2::new(0.0, 6.28),
    );

    let initial_count = system.body_count();
    let initial_mass = system.total_planet_mass();

    let events = vec![CollisionEvent {
        body_a: id_a,
        body_b: id_b,
        separation: 0.001,
        collision_radius: 0.02,
    }];

    resolve_collisions(&mut system, events);

    // Should have one less body
    assert_eq!(system.body_count(), initial_count - 1);

    // Total mass should be conserved
    let final_mass = system.total_planet_mass();
    assert!((final_mass - initial_mass).abs() < 1e-10);
}

#[test]
fn test_resolve_multiple_collisions() {
    let star = main_sequence_star(1.0, 0.0, 4_600.0);
    let mut system = SystemState::new(star);

    let id_a = system.add_body(0.001, 0.01, Point2::new(1.0, 0.0), Vector2::new(0.0, 6.28));
    let id_b = system.add_body(
        0.001,
        0.01,
        Point2::new(1.001, 0.0),
        Vector2::new(0.0, 6.28),
    );
    let id_c = system.add_body(0.001, 0.01, Point2::new(2.0, 0.0), Vector2::new(0.0, 4.44));
    let id_d = system.add_body(
        0.001,
        0.01,
        Point2::new(2.001, 0.0),
        Vector2::new(0.0, 4.44),
    );

    let initial_count = system.body_count();
    let initial_mass = system.total_planet_mass();

    let events = vec![
        CollisionEvent {
            body_a: id_a,
            body_b: id_b,
            separation: 0.001,
            collision_radius: 0.02,
        },
        CollisionEvent {
            body_a: id_c,
            body_b: id_d,
            separation: 0.001,
            collision_radius: 0.02,
        },
    ];

    resolve_collisions(&mut system, events);

    // Should have two less bodies (two pairs merged)
    assert_eq!(system.body_count(), initial_count - 2);

    // Total mass should be conserved
    let final_mass = system.total_planet_mass();
    assert!((final_mass - initial_mass).abs() < 1e-10);
}

#[test]
fn test_resolve_cascade_collision() {
    let star = main_sequence_star(1.0, 0.0, 4_600.0);
    let mut system = SystemState::new(star);

    let id_a = system.add_body(0.001, 0.01, Point2::new(1.0, 0.0), Vector2::new(0.0, 6.28));
    let id_b = system.add_body(
        0.001,
        0.01,
        Point2::new(1.001, 0.0),
        Vector2::new(0.0, 6.28),
    );
    let id_c = system.add_body(
        0.001,
        0.01,
        Point2::new(1.002, 0.0),
        Vector2::new(0.0, 6.28),
    );

    let initial_mass = system.total_planet_mass();

    // Create events where all three bodies overlap
    // The sorting and consumed tracking should handle this correctly
    let events = vec![
        CollisionEvent {
            body_a: id_a,
            body_b: id_b,
            separation: 0.001,
            collision_radius: 0.02,
        },
        CollisionEvent {
            body_a: id_b,
            body_b: id_c,
            separation: 0.001,
            collision_radius: 0.02,
        },
        CollisionEvent {
            body_a: id_a,
            body_b: id_c,
            separation: 0.002,
            collision_radius: 0.02,
        },
    ];

    resolve_collisions(&mut system, events);

    // Should only process valid mergers (consumed bodies skipped)
    // Closest pair (a,b) merges first, then (merged,c) cannot happen
    // because b is already consumed
    assert!(system.body_count() <= 2);

    // Mass should still be conserved
    let final_mass = system.total_planet_mass();
    assert!((final_mass - initial_mass).abs() < 1e-10);
}

#[test]
fn test_remove_star_collisions() {
    let star = main_sequence_star(1.0, 0.0, 4_600.0);
    let mut system = SystemState::new(star);

    let id_close = system.add_body(
        0.001,
        0.01,
        Point2::new(0.001, 0.0),
        Vector2::new(0.0, 50.0),
    );
    let id_far = system.add_body(0.001, 0.01, Point2::new(1.0, 0.0), Vector2::new(0.0, 6.28));

    remove_star_collisions(&mut system, vec![id_close]);

    assert_eq!(system.body_count(), 1);
    assert!(system.get_body(id_far).is_some());
    assert!(system.get_body(id_close).is_none());
}

#[test]
fn test_remove_ejections() {
    let star = main_sequence_star(1.0, 0.0, 4_600.0);
    let mut system = SystemState::new(star);

    let id_ejected = system.add_body(
        0.001,
        0.01,
        Point2::new(100.0, 0.0),
        Vector2::new(0.0, 20.0),
    );
    let id_bound = system.add_body(0.001, 0.01, Point2::new(1.0, 0.0), Vector2::new(0.0, 6.28));

    remove_ejections(&mut system, vec![id_ejected]);

    assert_eq!(system.body_count(), 1);
    assert!(system.get_body(id_bound).is_some());
    assert!(system.get_body(id_ejected).is_none());
}

#[test]
fn test_empty_events_list() {
    let star = main_sequence_star(1.0, 0.0, 4_600.0);
    let mut system = SystemState::new(star);

    system.add_body(0.001, 0.01, Point2::new(1.0, 0.0), Vector2::new(0.0, 6.28));

    let initial_count = system.body_count();

    resolve_collisions(&mut system, vec![]);

    assert_eq!(system.body_count(), initial_count);
}

#[test]
fn test_momentum_conservation_in_resolution() {
    let star = main_sequence_star(1.0, 0.0, 4_600.0);
    let mut system = SystemState::new(star);

    let id_a = system.add_body(0.001, 0.01, Point2::new(1.0, 0.0), Vector2::new(1.0, 5.0));
    let id_b = system.add_body(
        0.002,
        0.01,
        Point2::new(1.001, 0.0),
        Vector2::new(-0.5, 3.0),
    );

    let initial_momentum = system.total_momentum();

    let events = vec![CollisionEvent {
        body_a: id_a,
        body_b: id_b,
        separation: 0.001,
        collision_radius: 0.02,
    }];

    resolve_collisions(&mut system, events);

    let final_momentum = system.total_momentum();

    // Momentum should be conserved
    assert!((final_momentum.x - initial_momentum.x).abs() < 1e-10);
    assert!((final_momentum.y - initial_momentum.y).abs() < 1e-10);
}
