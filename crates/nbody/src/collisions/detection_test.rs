use nalgebra::{Point2, Vector2};
use stellar::generation::main_sequence_star;
use units::Mass;

use crate::collisions::detection::*;
use crate::collisions::{CollisionCriteria, CollisionDetector, DirectDetector, TreeDetector};
use crate::state::SystemState;

fn make_test_system() -> SystemState {
    let star = main_sequence_star(1.0, 0.0, 4_600.0);
    SystemState::new(star)
}

#[test]
fn test_hill_radius_earth() {
    let earth_mass = Mass::from_earth_masses(1.0).to_solar_masses();
    let r_hill = hill_radius(earth_mass, 1.0, 1.0);

    // Earth's Hill radius is about 0.01 AU
    assert!((r_hill - 0.01).abs() < 0.001);
}

#[test]
fn test_mutual_hill_radius() {
    let m = Mass::from_earth_masses(1.0).to_solar_masses();

    let r_hill = mutual_hill_radius(m, 1.0, m, 1.5, 1.0);

    // Should be positive and reasonable
    assert!(r_hill > 0.0);
    assert!(r_hill < 0.1); // Less than 0.1 AU
}

#[test]
fn test_direct_detector_no_collisions() {
    let mut system = make_test_system();

    // Add two well-separated bodies
    system.add_body(0.001, 0.01, Point2::new(1.0, 0.0), Vector2::new(0.0, 6.28));
    system.add_body(0.001, 0.01, Point2::new(2.0, 0.0), Vector2::new(0.0, 4.44));

    let detector = DirectDetector;
    let criteria = CollisionCriteria::default();
    let collisions = detector.detect(&system, &criteria);

    assert_eq!(collisions.len(), 0);
}

#[test]
fn test_direct_detector_physical_collision() {
    let mut system = make_test_system();

    // Add two bodies very close together (physical collision)
    system.add_body(0.001, 0.01, Point2::new(1.0, 0.0), Vector2::new(0.0, 6.28));
    system.add_body(
        0.001,
        0.01,
        Point2::new(1.015, 0.0),
        Vector2::new(0.0, 6.28),
    ); // 0.015 AU apart, radii sum to 0.02 AU

    let detector = DirectDetector;
    let criteria = CollisionCriteria {
        hill_fraction: 0.0, // Only physical
        physical_collision: true,
    };
    let collisions = detector.detect(&system, &criteria);

    assert_eq!(collisions.len(), 1);
    assert!(collisions[0].separation < collisions[0].collision_radius);
}

#[test]
fn test_direct_detector_hill_collision() {
    let mut system = make_test_system();

    // Add two bodies within Hill sphere but not physically touching
    let m = Mass::from_earth_masses(1.0).to_solar_masses();
    system.add_body(m, 0.0001, Point2::new(1.0, 0.0), Vector2::new(0.0, 6.28));
    system.add_body(m, 0.0001, Point2::new(1.005, 0.0), Vector2::new(0.0, 6.28)); // Within Hill sphere

    let detector = DirectDetector;
    let criteria = CollisionCriteria {
        hill_fraction: 0.5,
        physical_collision: true,
    };
    let collisions = detector.detect(&system, &criteria);

    assert!(collisions.len() > 0);
}

#[test]
fn test_tree_detector_no_collisions() {
    let mut system = make_test_system();

    // Add several well-separated bodies
    system.add_body(0.001, 0.01, Point2::new(1.0, 0.0), Vector2::new(0.0, 6.28));
    system.add_body(0.001, 0.01, Point2::new(2.0, 0.0), Vector2::new(0.0, 4.44));
    system.add_body(0.001, 0.01, Point2::new(0.5, 0.0), Vector2::new(0.0, 8.87));

    let detector = TreeDetector;
    let criteria = CollisionCriteria::default();
    let collisions = detector.detect(&system, &criteria);

    assert_eq!(collisions.len(), 0);
}

#[test]
fn test_tree_detector_finds_collision() {
    let mut system = make_test_system();

    // Add two bodies very close together
    system.add_body(0.001, 0.01, Point2::new(1.0, 0.0), Vector2::new(0.0, 6.28));
    system.add_body(
        0.001,
        0.01,
        Point2::new(1.015, 0.0),
        Vector2::new(0.0, 6.28),
    );

    let detector = TreeDetector;
    let criteria = CollisionCriteria {
        hill_fraction: 0.0,
        physical_collision: true,
    };
    let collisions = detector.detect(&system, &criteria);

    assert_eq!(collisions.len(), 1);
}

#[test]
fn test_tree_vs_direct_same_results() {
    let mut system = make_test_system();

    // Add several bodies with some colliding
    system.add_body(0.001, 0.01, Point2::new(1.0, 0.0), Vector2::new(0.0, 6.28));
    system.add_body(
        0.001,
        0.01,
        Point2::new(1.015, 0.0),
        Vector2::new(0.0, 6.28),
    );
    system.add_body(0.001, 0.01, Point2::new(2.0, 0.0), Vector2::new(0.0, 4.44));
    system.add_body(
        0.001,
        0.01,
        Point2::new(2.015, 0.0),
        Vector2::new(0.0, 4.44),
    );

    let criteria = CollisionCriteria {
        hill_fraction: 0.0,
        physical_collision: true,
    };

    let collisions_direct = DirectDetector.detect(&system, &criteria);
    let collisions_tree = TreeDetector.detect(&system, &criteria);

    // Should find the same number of collisions
    assert_eq!(collisions_direct.len(), collisions_tree.len());
    assert_eq!(collisions_direct.len(), 2); // Two pairs colliding
}

#[test]
fn test_detect_star_collisions() {
    let mut system = make_test_system();

    // Add a body very close to the star
    system.add_body(
        0.001,
        0.01,
        Point2::new(0.001, 0.0),
        Vector2::new(0.0, 50.0),
    );

    // Add a normal body
    system.add_body(0.001, 0.01, Point2::new(1.0, 0.0), Vector2::new(0.0, 6.28));

    let star_collisions = detect_star_collisions(&system);

    assert_eq!(star_collisions.len(), 1);
}

#[test]
fn test_detect_ejections() {
    use crate::forces::G;

    let mut system = make_test_system();

    let r = 100.0; // Far out
    let star_mass = system.star.mass.to_solar_masses();
    let mu = G * star_mass;

    // Escape velocity at r
    let v_escape = (2.0 * mu / r).sqrt();

    // Add a body moving faster than escape velocity
    system.add_body(
        0.001,
        0.01,
        Point2::new(r, 0.0),
        Vector2::new(0.0, v_escape * 1.5),
    );

    // Add a bound body
    system.add_body(0.001, 0.01, Point2::new(1.0, 0.0), Vector2::new(0.0, 6.28));

    let ejections = detect_ejections(&system, 50.0);

    assert_eq!(ejections.len(), 1);
}

#[test]
fn test_no_ejections_below_escape_radius() {
    use crate::forces::G;

    let mut system = make_test_system();

    let r = 10.0; // Close to star
    let star_mass = system.star.mass.to_solar_masses();
    let mu = G * star_mass;
    let v_escape = (2.0 * mu / r).sqrt();

    // Fast but inside escape_radius
    system.add_body(
        0.001,
        0.01,
        Point2::new(r, 0.0),
        Vector2::new(0.0, v_escape * 1.5),
    );

    let ejections = detect_ejections(&system, 50.0); // escape_radius > r

    assert_eq!(ejections.len(), 0); // Not beyond escape radius
}

#[test]
fn test_empty_system() {
    let system = make_test_system();

    let collisions_direct = DirectDetector.detect(&system, &CollisionCriteria::default());
    let collisions_tree = TreeDetector.detect(&system, &CollisionCriteria::default());
    let star_collisions = detect_star_collisions(&system);
    let ejections = detect_ejections(&system, 100.0);

    assert_eq!(collisions_direct.len(), 0);
    assert_eq!(collisions_tree.len(), 0);
    assert_eq!(star_collisions.len(), 0);
    assert_eq!(ejections.len(), 0);
}

#[test]
fn test_single_body_no_collisions() {
    let mut system = make_test_system();
    system.add_body(0.001, 0.01, Point2::new(1.0, 0.0), Vector2::new(0.0, 6.28));

    let collisions_direct = DirectDetector.detect(&system, &CollisionCriteria::default());
    let collisions_tree = TreeDetector.detect(&system, &CollisionCriteria::default());

    assert_eq!(collisions_direct.len(), 0);
    assert_eq!(collisions_tree.len(), 0);
}
